"""Controlled experimental inputs around the byte-preserved HemoFlow solver."""
from pathlib import Path
import sys
import math
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'source_frozen'))
from config import Config as DesktopConfig
from solver import FlowSolver
from audit_math import build_stencil

class Config(DesktopConfig):
    """Headless capacity override; all frozen physical validation stays active.

    The preserved validator's LAST check is a desktop-preview capacity ceiling.
    Only that exact error is eligible, with a bounded one-million-cell ceiling.
    """
    def validate(self):
        try:
            return super().validate()
        except ValueError as exc:
            if str(exc) != "This grid is too large for the desktop preview. Reduce length or cells_across.":
                raise
            nominal_cells = self.cells_across**2 * self.length / self.diameter
            if not 600_000 < nominal_cells <= 1_000_000:
                raise ValueError("Headless experiment exceeds the 1,000,000 nominal-cell capacity limit") from exc
            return self

VERSION='grid-refinement-2026-09-28.2' 
REFLECT=np.array([0,1,4,3,2,8,7,6,5])
STATIONS=np.arange(20.,55.0001,.5)


def plan(extended=False,smoke=False):
    cases=[]
    for n,dt in [(40,3e-6),(80,1.5e-6)]:
        cases.append(dict(key=f'healthy_steady_N{n:03}',kind='healthy',cells=n,dt_s=dt,
                          epsilon=0.,duration_s=.002 if smoke else 2.,cycles=None))
    for n,eps in [(80,0.),(80,.001),(80,-.001),(120,.001),(120,-.001)]:
        tag='zero' if eps==0 else ('plus' if eps>0 else 'minus')+'001'
        cases.append(dict(key=f'D75_N{n:03}_{tag}',kind='severe',cells=n,
            dt_s=.3e-6 if n==80 else .2e-6,epsilon=eps,
            duration_s=.002 if smoke else 10*60/72,cycles=10))
    if extended:
        for eps in [.0001,-.0001]:
            cases.append(dict(key='D75_N080_'+('plus' if eps>0 else 'minus')+'0001',kind='severe',
                cells=80,dt_s=.3e-6,epsilon=eps,duration_s=.002 if smoke else 10*60/72,cycles=10))
    return cases


def make_solver(spec,threads=1):
    healthy=spec['kind']=='healthy'
    c=Config(length_mm=spec.get('length_mm',60.),geometry_mode='healthy' if healthy else 'single',stenosis_percent=0 if healthy else 75,
        plaque_count=0 if healthy else 1,randomize_count=False,randomize_on_launch=False,
        cells_across=spec['cells'],pulsatility_percent=0 if healthy else 40,
        boundary_model='section_impedance',particles_visible=False,compute_threads=threads,
        time_scale=.3 if healthy else .1512)
    # Compute requested dt from the same scaling used by the frozen constructor.
    from geometry import make_geometry
    g=make_geometry(c)
    ul=min(.045,.1*(g.min_gap_cells/c.cells_across)/(1.5*(1+c.pulsatility_percent/100)))
    c.time_scale=spec['dt_s']*c.mean_velocity/(ul*g.dx)
    s=FlowSolver(c)
    # The frozen kernel overwrites every fluid entry before swapping buffers.
    # Initialize its scratch buffer solely to make t=0 diagnostic storage finite.
    s.other[:]=s.f
    s.set_threads(threads)
    if not math.isclose(s.dt,spec['dt_s'],rel_tol=1e-12):raise ValueError('Requested timestep was not realized')
    if not .501<=s.tau<=2:raise ValueError('Relaxation outside frozen solver constraints')
    if not np.array_equal(s.geometry.solid,s.geometry.solid[::-1]):raise ValueError('The geometry mask is not mirror symmetric')
    return s


def reflected_populations(f,ny,nx):
    return f.reshape(ny,nx,9)[::-1,:,REFLECT].reshape(-1,9).copy()


def initialize_pair(s,epsilon):
    """Negative branch starts from the reflected positive-branch initial state.

    This mirrors only roundoff-scale asymmetry already present in the nominally
    symmetric initialization. It is logged, not treated as an unchanged initial
    bit pattern. The zero-control path does not alter any state or inlet value.
    """
    original=s.profile.copy()
    if epsilon==0:
        return original,original.copy(),dict(epsilon=0.,profile_relative_flux_change=0.,initial_state_reflected=False)
    if not 0<abs(epsilon)<=.01:raise ValueError('Unsupported perturbation amplitude')
    plus=original*(1+abs(epsilon)*2*s.geometry.y/s.config.diameter)
    plus*=original.sum()/plus.sum()
    driven=plus if epsilon>0 else plus[::-1].copy()
    if epsilon<0:
        s.f=reflected_populations(s.f,s.ny,s.nx)
        shape=(s.ny,s.nx)
        s.rho=s.rho.reshape(shape)[::-1].copy().ravel()
        s.ux=s.ux.reshape(shape)[::-1].copy().ravel()
        s.uy=-s.uy.reshape(shape)[::-1].copy().ravel()
    error=float((driven.sum()-original.sum())/original.sum())
    if abs(error)>1e-14:raise ValueError('Perturbation changed total inlet flow')
    return original,driven,dict(epsilon=epsilon,profile_relative_flux_change=error,
        initial_state_reflected=epsilon<0,description='Axial inlet skew 1+epsilon*2y/D, discretely flux normalized, active for first heartbeat only; zero transverse inlet velocity.')


def advance_to(s,target_iteration,original,driven,cutoff):
    """Never step across the forcing cutoff without restoring the profile."""
    while s.iteration<target_iteration:
        forced=s.iteration<cutoff
        s.profile=driven if forced else original
        end=min(target_iteration,cutoff) if forced else target_iteration
        s.step(int(end-s.iteration))


class WallSampler:
    def __init__(self,s):
        self.s=s;self.entries=[]
        g=s.geometry
        for wall,surface,slope in [('lower',g.lower,s.lower_slope),('upper',g.upper,s.upper_slope)]:
            for x in STATIONS:
                j=int(np.argmin(abs(g.x-x/1000)))
                for radius in (4.,5.):
                    stencil=build_stencil(g.x,g.y,g.solid,g.x[j],surface[j],slope[j],wall,radius)
                    if stencil.reason:raise ValueError(f'Invalid wall stencil: {wall} {x}: {stencil.reason}')
                    self.entries.append((wall,x,j,radius,stencil))

    def sample(self,d):
        rows=[]
        for wall,x,j,radius,stencil in self.entries:
            r=stencil.evaluate(d['u'],d['v'],self.s.config.viscosity)
            if not r['valid']:raise ValueError(f'Invalid derivative at {wall} {x}: {r}')
            rows.append(dict(time_s=self.s.time,wall=wall,x_mm=float(x),radius_cells=radius,
                legacy_pa=float(d['wss_low' if wall=='lower' else 'wss_high'][j]),**r))
        return rows
