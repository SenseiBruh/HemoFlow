"""A single owner of mutable physics state, separate from GUI drawing.

The viewer only receives immutable snapshots. Reset requests replace one
pending request, snapshots use an eight-entry ring, and injections are bounded.
"""
from collections import deque
from copy import deepcopy
from dataclasses import dataclass, replace
from threading import Condition, Thread
from time import perf_counter
import numpy as np
from particles import Tracers
from solver import FlowSolver


@dataclass(frozen=True)
class Snapshot:
    epoch: int
    sequence: int
    wall_time: float
    time: float
    config: object
    geometry: object
    fields: dict
    positions: np.ndarray
    particle_ids: np.ndarray
    trails: np.ndarray
    threads: int
    iteration: int
    dt: float
    particles_enabled: bool = True


class FlowWorker(Thread):
    def __init__(self, config, publish_interval=.04):
        super().__init__(name="HemoFlow physics",daemon=True)
        self.condition=Condition()
        self.frames=deque(maxlen=8)
        self.injections=deque(maxlen=8)
        self.epoch=0
        self.pending=(0,deepcopy(config))
        self.running=True
        self.stopping=False
        self.state="preparing"
        self.error=""
        self.publish_interval=publish_interval
        self.threads_cache={}
        self.total_steps=0
        self.particles_enabled=bool(getattr(config,"particles_visible",True))
        self.particle_update=None
        self.particle_restart=False
        self.recorder=None

    def reset(self, config):
        config.validate()
        with self.condition:
            self.epoch+=1
            self.pending=(self.epoch,deepcopy(config))
            self.frames.clear()
            self.injections.clear()
            self.particle_update=None
            self.particle_restart=False
            self.particles_enabled=bool(getattr(config,"particles_visible",True))
            self.state="preparing"
            self.error=""
            self.condition.notify_all()
            return self.epoch

    def set_running(self,running):
        with self.condition:
            self.running=bool(running)
            self.condition.notify_all()

    def set_particles_enabled(self, enabled):
        """Enable/disable tracer advection and rendering without resetting CFD."""
        enabled=bool(enabled)
        with self.condition:
            if enabled != self.particles_enabled:
                self.particles_enabled=enabled
                # Re-seed a bounded tracer set when particles are turned back
                # on; the fluid state and simulated time are left untouched.
                self.particle_restart=enabled
            self.condition.notify_all()

    def set_particle_parameters(self, model, diameter_um, density_kg_m3):
        """Queue one-way particle-property changes without rebuilding the flow."""
        model=str(model)
        diameter_um=float(diameter_um)
        density_kg_m3=float(density_kg_m3)
        with self.condition:
            base = self.pending[1] if self.pending is not None else None
            if base is None:
                # The worker owns the active solver config; this tuple is
                # consumed in its thread before the next batch.
                base = getattr(self, "active_config", None)
            if base is not None:
                updated=replace(base, particle_model=model,
                                particle_diameter_um=diameter_um,
                                particle_density_kg_m3=density_kg_m3).validate()
                self.particle_update=updated
            self.condition.notify_all()

    def set_recorder(self, recorder):
        """Attach/detach a recorder; recording never changes solver physics."""
        with self.condition:
            old=self.recorder
            self.recorder=recorder
            self.condition.notify_all()
        if old is not None and old is not recorder:
            old.close()

    def inject(self,x_mm,y_mm):
        with self.condition:
            self.injections.append((self.epoch,float(x_mm),float(y_mm)))
            self.condition.notify_all()

    def status(self):
        with self.condition:
            return self.state,self.error,tuple(self.frames)

    def latest(self):
        """Return the newest published snapshot without removing it."""
        with self.condition:
            return self.frames[-1] if self.frames else None

    def stop(self):
        with self.condition:
            self.stopping=True
            self.state="stopping"
            self.condition.notify_all()

    def _publish(self,solver,tracers,epoch,sequence):
        d=solver.diagnostics()
        for value in d.values():
            if isinstance(value,np.ndarray):
                value.setflags(write=False)
        if tracers is None:
            positions=np.empty((0,2),dtype=float)
            ids=np.empty((0,),dtype=np.uint32)
            trails=np.empty((0,0,2),dtype=float)
        else:
            positions=tracers.positions(solver)
            ids=tracers.ids.copy()
            trails=tracers.trails().copy()
        for value in (positions,ids,trails):value.setflags(write=False)
        snapshot=Snapshot(epoch,sequence,perf_counter(),solver.time,deepcopy(solver.config),
                          solver.geometry,d,positions,ids,trails,solver.threads_used,
                          solver.iteration,solver.dt,tracers is not None)
        with self.condition:
            if epoch == self.epoch and not self.stopping:
                self.frames.append(snapshot)
            recorder=self.recorder
        if recorder is not None:
            recorder.record(solver,d,epoch,tracers)

    def run(self):
        solver=tracers=None
        sequence=0
        active_epoch=-1
        last_publish=0
        # A modest block amortizes Python/Numba call overhead while keeping
        # pause, plaque edits, and reset requests responsive.
        batch_steps=96
        while True:
            with self.condition:
                if self.stopping:
                    self.state="stopped"
                    if self.recorder is not None:
                        self.recorder.close()
                    return
                pending=self.pending
                self.pending=None
                injections=tuple(self.injections)
                self.injections.clear()
                running=self.running
                particle_update=self.particle_update
                self.particle_update=None
                particle_restart=self.particle_restart
                self.particle_restart=False
                particles_enabled=self.particles_enabled
                if pending is None and particle_update is None and not particle_restart and not running and not injections:
                    self.state="error" if self.error else "paused"
                    self.condition.wait(.1)
                    continue
            try:
                if pending is not None:
                    active_epoch,c=pending
                    self.active_config=c
                    solver=FlowSolver(c)
                    key=(solver.nx,solver.ny,c.compute_threads)
                    if key in self.threads_cache:
                        solver.set_threads(self.threads_cache[key])
                    else:
                        self.threads_cache[key]=solver.tune_threads()
                    solver.step(1)
                    tracers=Tracers(solver,trail_length=28) if particles_enabled else None
                    if tracers is not None:
                        tracers.advance(solver,1)
                    sequence=0
                    self._publish(solver,tracers,active_epoch,sequence)
                    last_publish=perf_counter()
                if solver is not None and particle_update is not None:
                    solver.config=particle_update
                    self.active_config=particle_update
                    if particles_enabled:
                        tracers=Tracers(solver,trail_length=28)
                    else:
                        tracers=None
                if solver is not None and particle_restart:
                    tracers=Tracers(solver,trail_length=28) if particles_enabled else None
                if solver is not None and not particles_enabled:
                    tracers=None
                if solver is not None and solver.config.particles_visible != particles_enabled:
                    solver.config=replace(solver.config, particles_visible=particles_enabled).validate()
                    self.active_config=solver.config
                if solver is None:continue
                for epoch,x,y in injections:
                    if epoch==active_epoch and tracers is not None:tracers.inject(solver,x,y)
                if running:
                    solver.step(batch_steps)
                    solver.check_stability()
                    if tracers is not None:tracers.advance(solver,batch_steps)
                    self.total_steps+=batch_steps
                now=perf_counter()
                with self.condition:
                    self.state="running" if self.running else "paused"
                    # Publish the state reached just before a pause, then freeze.
                    just_paused=running and not self.running
                if now-last_publish>=self.publish_interval or injections or just_paused or not running:
                    if tracers is not None:tracers.record(solver)
                    sequence+=1
                    self._publish(solver,tracers,active_epoch,sequence)
                    last_publish=now
            except Exception as exc:
                with self.condition:
                    if active_epoch==self.epoch:
                        self.error=str(exc)
                        self.state="error"
                        self.running=False


def interpolate_frames(frames, wall_time, paused=False, delay=.065):
    """Interpolate computed points, never extrapolate or cross a particle recycle."""
    if not frames:return None
    target=wall_time-delay
    a=b=frames[-1]
    alpha=1.0
    if not paused:
        for left,right in zip(frames[:-1],frames[1:]):
            if left.epoch==right.epoch and left.wall_time<=target<=right.wall_time:
                a,b=left,right
                alpha=(target-left.wall_time)/max(right.wall_time-left.wall_time,1e-9)
                break
        else:
            if target<frames[0].wall_time:a=b=frames[0]
    positions=b.positions.copy()
    if a is not b:
        if a.particle_ids.shape == b.particle_ids.shape and a.positions.shape == b.positions.shape:
            same=a.particle_ids==b.particle_ids
            positions[same]=(1-alpha)*a.positions[same]+alpha*b.positions[same]
            positions[~same]=np.nan
            g=b.geometry
            finite=np.isfinite(positions).all(axis=1)
            ix=np.clip(np.rint(np.nan_to_num(positions[:,0])/1000/g.dx).astype(int),0,len(g.x)-1)
            iy=np.clip(np.rint((np.nan_to_num(positions[:,1])/1000-g.y[0])/g.dx).astype(int),0,len(g.y)-1)
            # A straight interpolation across a curved wall must not show a tracer in plaque.
            invalid=finite & g.solid[iy,ix]
            positions[invalid]=a.positions[invalid]
        # Particle visibility changes use a different array size.  Use the
        # newest state instead of attempting a broadcast across the toggle.
    return a,b,float(alpha),positions,(1-alpha)*a.time+alpha*b.time
