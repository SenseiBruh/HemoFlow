"""Interactive HemoFlow viewer.

The physics thread publishes snapshots while this Matplotlib thread only draws
and interpolates those snapshots. A plaque drag is a live geometry preview;
releasing the mouse commits it and starts a new fluid solution for that geometry.
"""
from collections import deque
from copy import deepcopy
from dataclasses import asdict, replace
from datetime import datetime
import json
from pathlib import Path
from time import perf_counter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.widgets import Button, Slider
from matplotlib.lines import Line2D
from engine import FlowWorker, interpolate_frames
from geometry import make_geometry, plaque_profile, stenosis_report
from physics import calculate_pressure_drop, calculate_wall_shear_stress, calculate_womersley_number, particle_relaxation_time
from recorder import RunRecorder
from rendering import BlitRenderer
from waveform import phase_label, pulse_curve

BG = "#0b1220"
PANEL = "#111d30"
TEXT = "#dfebf6"
MUTED = "#8fa5bc"
CYAN = "#55e5e2"
ORANGE = "#ffb96b"
SELECT = "#fff09a"
WALL = "#edaa8b"


class Dashboard:
    def __init__(self, config):
        config.validate()
        self.config = deepcopy(config)
        self.worker = FlowWorker(self.config, publish_interval=.035)
        self.worker.start()
        self.message = "Preparing the threaded solver; the window will update when its first snapshot is ready."
        self.paused = False
        self.lower_mode = "vorticity"
        self.animation = None
        self.renderer = None
        self.last_trail_key = None
        self.last_draw_time = None
        self.frame_interval = None
        self.history = deque(maxlen=600)
        self.last_history_time = -np.inf
        self.wall_start = perf_counter()
        self.frame = 0
        self.selected_id = None
        self.dragging = False
        self.drag_moved = False
        self.drag_mode = None
        self.drag_origin = None
        self.drag_base = None
        self.preview_plaque = None
        self.recording = None
        self.particles_enabled = bool(self.config.particles_visible)
        self.pending_heart_rate = float(self.config.heart_rate)
        self.snapshot = self._wait_for_first_snapshot()
        if self.snapshot is None:
            self.worker.stop()
            raise RuntimeError("The worker did not produce an initial snapshot.")
        self.geometry = self.snapshot.geometry
        self.fig = plt.figure(figsize=(15.2, 10.1), facecolor=BG)
        try:
            self.fig.canvas.manager.set_window_title("HemoFlow v3.1 — continuous flow laboratory")
        except AttributeError:
            pass
        self.fig.text(.065, .963, "HEMOFLOW", fontsize=24, weight="bold", color=TEXT)
        self.fig.text(.237, .964, "CONTINUOUS FLOW LABORATORY", fontsize=10, color=CYAN)
        self.fig.text(.92, .964, "THREADED SOLVER  /  2D  /  NEWTONIAN", fontsize=8, color=MUTED, ha="right")
        self.status = self.fig.text(.065, .922, "", fontsize=10.5, color=TEXT)
        self.flow_ax = self.fig.add_axes([.065, .653, .815, .216])
        self.lower_ax = self.fig.add_axes([.065, .411, .815, .162])
        self.flow_cax = self.fig.add_axes([.895, .653, .011, .216])
        self.lower_cax = self.fig.add_axes([.895, .411, .011, .162])
        self.pressure_ax = self.fig.add_axes([.065, .162, .247, .151])
        self.wss_ax = self.fig.add_axes([.37, .162, .247, .151])
        self.pulse_ax = self.fig.add_axes([.675, .162, .23, .151])
        for ax in (self.flow_ax, self.lower_ax, self.pressure_ax, self.wss_ax, self.pulse_ax):
            self._style(ax)
        self._setup_maps()
        self.pressure_line, = self.pressure_ax.plot([], [], color=CYAN, lw=1.5)
        self.wss_low, = self.wss_ax.plot([], [], color=CYAN, lw=1.1, label="Lower wall")
        self.wss_high, = self.wss_ax.plot([], [], color=ORANGE, lw=1.1, label="Upper wall")
        legend = self.wss_ax.legend(fontsize=7, loc="upper right", facecolor=PANEL, edgecolor=PANEL)
        for item in legend.get_texts():
            item.set_color(TEXT)
        self.pulse_line, = self.pulse_ax.plot([], [], color=ORANGE, lw=1.5, label="inlet speed")
        self.pulse_marker, = self.pulse_ax.plot([], [], marker="o", color=SELECT, ms=5, lw=0, zorder=8)
        self.pulse_phase = self.pulse_ax.text(.98, .90, "", transform=self.pulse_ax.transAxes,
                                             ha="right", va="top", color=SELECT, fontsize=8, weight="bold")
        for ax, title, ylabel in (
                (self.pressure_ax, "PRESSURE ALONG VESSEL", "Mean pressure (Pa)"),
                (self.wss_ax, "WALL SHEAR ESTIMATE", "|Shear| (Pa)"),
                (self.pulse_ax, "INLET WAVEFORM", "Mean speed (m/s)")):
            ax.set_title(title, fontsize=9, color=TEXT, loc="left", pad=5, weight="bold")
            ax.set_ylabel(ylabel, fontsize=8)
            ax.grid(alpha=.10)
        self.pressure_ax.set_xlabel("Distance (mm)", fontsize=8)
        self.wss_ax.set_xlabel("Distance (mm)", fontsize=8)
        self.pulse_ax.set_xlabel("Cardiac cycle (%)", fontsize=8)
        self.reference = self.fig.text(.065, .354, "", fontsize=8, color=MUTED)
        self.plaque_info = self.fig.text(.065, .334, "", fontsize=8.5, color=SELECT)
        self.footer = self.fig.text(.065, .013, "", fontsize=8, color=MUTED)
        self.buttons = []
        self.button_by_name = {}

        def add_button(name, x, width, label, callback, y=.067, primary=True):
            ax = self.fig.add_axes([x, y, width, .030])
            button = Button(ax, label, color="#1b2c43" if primary else "#16263a",
                            hovercolor="#304761" if primary else "#2d445e")
            button.label.set_color(TEXT if primary else MUTED)
            button.label.set_fontsize(8.0 if primary else 7.3)
            button.on_clicked(callback)
            self.buttons.append(button)
            self.button_by_name[name] = button

        add_button("pause", .065, .085, "Pause", self.toggle_pause)
        add_button("layout", .158, .095, "New layout", self.regenerate)
        add_button("add", .262, .070, "Add plaque", self.add_plaque)
        add_button("delete", .338, .100, "Delete selected", self.delete_selected)
        add_button("map", .445, .090, "Map: vorticity", self.toggle_map)
        add_button("particles", .541, .100, "Particles: ON", self.toggle_particles)
        add_button("record", .647, .090, "Record: OFF", self.toggle_recording)
        add_button("snapshot", .743, .095, "Save snapshot", self.save_snapshot)
        add_button("settings", .844, .095, "Save settings", self.save_settings)

        add_button("height_minus", .065, .042, "H −", lambda event: self.adjust_selected("height", -5), y=.032, primary=False)
        add_button("height_plus", .112, .042, "H +", lambda event: self.adjust_selected("height", 5), y=.032, primary=False)
        add_button("pulse", .160, .112, "Pulse: systolic", self.toggle_pulse_shape, y=.032, primary=False)
        add_button("bpm_minus", .278, .052, "BPM −", lambda event: self.adjust_heart_rate(-1), y=.032, primary=False)
        add_button("bpm_plus", .335, .052, "BPM +", lambda event: self.adjust_heart_rate(1), y=.032, primary=False)
        add_button("particle_minus", .395, .045, "P −", lambda event: self.adjust_particle_size(-2), y=.032, primary=False)
        add_button("particle_plus", .445, .045, "P +", lambda event: self.adjust_particle_size(2), y=.032, primary=False)
        add_button("particle_model", .495, .105, "Model: tracer", self.toggle_particle_model, y=.032, primary=False)

        hrax = self.fig.add_axes([.615, .041, .305, .012], facecolor=PANEL)
        self.heart_rate_slider = Slider(hrax, "BPM", 20, 240, valinit=self.config.heart_rate,
                                        valstep=1, color=ORANGE)
        self.heart_rate_slider.label.set_color(MUTED)
        self.heart_rate_slider.label.set_fontsize(7.2)
        self.heart_rate_slider.valtext.set_color(TEXT)
        self.heart_rate_slider.on_changed(self._preview_heart_rate)
        sax = self.fig.add_axes([.615, .020, .305, .012], facecolor=PANEL)
        # More latency gives the renderer a deeper interpolation buffer and
        # smoother particles; less latency makes the display respond sooner.
        self.playback = Slider(sax, "Smoothing", .25, 3.0, valinit=1.0, valstep=.25, color=CYAN)
        self.playback.label.set_color(MUTED)
        self.playback.label.set_fontsize(7.5)
        self.playback.valtext.set_color(TEXT)
        self.playback.on_changed(lambda _value: self.fig.canvas.draw_idle())
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)
        self.fig.canvas.mpl_connect("close_event", self.on_close)
        self.render()
        self.renderer = BlitRenderer(self.fig, self._animated_artists)

    def _animated_artists(self):
        return ([self.speed_map, self.lower_map, self.quiver, self.trail_lines, self.dots,
                 self.pressure_line, self.wss_low, self.wss_high, self.pulse_line,
                 self.pulse_marker, self.pulse_phase, self.status, self.reference,
                 self.plaque_info, self.footer] + self.wall_artists +
                [artist for group in self.plaque_artists.values() for artist in group])

    def _expand_limits(self, axis, values, nonnegative=False):
        values = np.asarray(values)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return
        lo, hi = float(finite.min()), float(finite.max())
        current_lo, current_hi = axis.get_ylim()
        if lo < current_lo or hi > current_hi or current_hi-current_lo <= 1.0:
            span = max(hi-lo, abs(hi)*.1, .1)
            limits = (0.0 if nonnegative else min(current_lo,lo-.15*span), max(current_hi,hi+.15*span))
            if limits != (current_lo,current_hi):
                axis.set_ylim(*limits)
                if self.renderer is not None:
                    self.renderer.invalidate()

    def _wait_for_first_snapshot(self, timeout=180):
        started = perf_counter()
        while perf_counter() - started < timeout:
            snap = self.worker.latest()
            state, error, _ = self.worker.status()
            if snap is not None:
                return snap
            if error:
                raise RuntimeError(error)
            # The solver's first JIT/thread benchmark is allowed to run without
            # freezing a GUI that has not been created yet.
            import time
            time.sleep(.02)
        raise TimeoutError("Timed out while preparing the numerical solver.")

    @staticmethod
    def _style(ax):
        ax.set_facecolor(PANEL)
        ax.tick_params(colors=MUTED, labelsize=8, length=3)
        ax.xaxis.label.set_color(MUTED)
        ax.yaxis.label.set_color(MUTED)
        for spine in ax.spines.values():
            spine.set_color("#293a50")

    def _remove_artist(self, artist):
        if artist is not None:
            try:
                artist.remove()
            except (ValueError, AttributeError):
                pass

    @staticmethod
    def _particle_marker_size(diameter_um):
        """Map micrometre diameter to readable screen area, not physical area."""
        return float(np.clip(2.0 + 0.18 * np.sqrt(max(float(diameter_um), .05)), 2.0, 22.0))

    def _setup_maps(self):
        s = self.snapshot
        g, c = s.geometry, s.config
        for artist in getattr(self, "wall_artists", []):
            self._remove_artist(artist)
        for artists in getattr(self, "plaque_artists", {}).values():
            for artist in artists:
                self._remove_artist(artist)
        for attr in ("speed_map", "lower_map", "trail_lines", "dots", "quiver"):
            self._remove_artist(getattr(self, attr, None))
        for attr in ("speed_bar", "lower_bar"):
            bar = getattr(self, attr, None)
            if bar is not None:
                try:
                    bar.remove()
                except (ValueError, AttributeError):
                    pass
        self.flow_ax.clear()
        self.lower_ax.clear()
        # Colorbar.remove() can detach a manually positioned cax from the
        # figure. Recreate it when that happens before clearing/reusing it.
        for attr, rect in (("flow_cax", [.895, .653, .011, .216]),
                           ("lower_cax", [.895, .411, .011, .162])):
            cax = getattr(self, attr, None)
            if cax is None or cax.get_figure() is None:
                cax = self.fig.add_axes(rect, facecolor=PANEL)
                setattr(self, attr, cax)
            else:
                cax.clear()
        self._style(self.flow_ax)
        self._style(self.lower_ax)
        self.flow_ax.set_xlim(0, g.x[-1]*1000)
        self.lower_ax.set_xlim(0, g.x[-1]*1000)
        self.flow_ax.set_ylim(-c.diameter_mm/2-.14, c.diameter_mm/2+.14)
        self.lower_ax.set_ylim(-c.diameter_mm/2-.14, c.diameter_mm/2+.14)
        for ax in (self.flow_ax, self.lower_ax):
            ax.set_xlabel("Distance along vessel (mm)", fontsize=8)
            ax.set_ylabel("y (mm)", fontsize=8)
        extent = [-g.dx*500, (g.x[-1]+g.dx/2)*1000,
                  (g.y[0]-g.dx/2)*1000, (g.y[-1]+g.dx/2)*1000]
        blank = np.ma.masked_where(g.solid, np.zeros(g.solid.shape))
        speed_cmap = plt.get_cmap("magma").copy()
        speed_cmap.set_bad(PANEL)
        self.speed_limit = max(.6, c.mean_velocity*(1+c.pulsatility_percent/100)*1.7*c.cells_across/max(g.min_gap_cells,1))
        self.speed_map = self.flow_ax.imshow(blank, extent=extent, origin="lower", aspect="auto",
                                              cmap=speed_cmap, vmin=0, vmax=self.speed_limit,
                                              interpolation="nearest", zorder=1)
        lower_cmap = plt.get_cmap("coolwarm").copy()
        lower_cmap.set_bad(PANEL)
        self.vort_limit = 8*c.mean_velocity/c.diameter
        self.lower_map = self.lower_ax.imshow(blank, extent=extent, origin="lower", aspect="auto",
                                               cmap=lower_cmap, vmin=-self.vort_limit, vmax=self.vort_limit,
                                               interpolation="nearest", zorder=1)
        self.wall_artists = []
        for ax in (self.flow_ax, self.lower_ax):
            self.wall_artists.extend([
                ax.fill_between(g.x*1000, -c.diameter_mm/2-.2, g.lower*1000,
                                facecolor="#794945", edgecolor="none", zorder=4),
                ax.fill_between(g.x*1000, g.upper*1000, c.diameter_mm/2+.2,
                                facecolor="#794945", edgecolor="none", zorder=4),
                ax.plot(g.x*1000, g.lower*1000, color=WALL, lw=1.0, zorder=5)[0],
                ax.plot(g.x*1000, g.upper*1000, color=WALL, lw=1.0, zorder=5)[0],
            ])
        self.flow_ax.set_title("VELOCITY  +  PARTICLE PATHS", loc="left", color=TEXT, fontsize=10, weight="bold", pad=10)
        self.flow_ax.set_title("Threaded physics  ·  click vessel to inject", loc="right", color=MUTED, fontsize=8, pad=10)
        self.lower_ax.set_title("VORTICITY  ·  signed local rotation (shear also contributes)", loc="left", color=TEXT, fontsize=9, pad=10)
        self.speed_bar = self.fig.colorbar(self.speed_map, cax=self.flow_cax)
        self.lower_bar = self.fig.colorbar(self.lower_map, cax=self.lower_cax)
        self.speed_bar.set_label("Speed (m/s)", color=MUTED, fontsize=8)
        self.lower_bar.set_label("Vorticity (1/s)", color=MUTED, fontsize=8)
        if self.lower_mode == "backflow":
            self.lower_map.set_clim(-self.speed_limit, self.speed_limit)
            self.lower_bar.set_label("Axial velocity (m/s)", color=MUTED, fontsize=8)
            self.lower_ax.set_title("AXIAL VELOCITY  ·  blue = backward, red = forward", loc="left",
                                    color=TEXT, fontsize=9, pad=10)
        else:
            self.lower_map.set_clim(-self.vort_limit, self.vort_limit)
        for bar in (self.speed_bar, self.lower_bar):
            bar.ax.tick_params(colors=MUTED, labelsize=7)
            bar.outline.set_edgecolor("#293a50")
        self.trail_lines = LineCollection([], linewidths=.6, colors=CYAN, alpha=.40, zorder=6)
        self.flow_ax.add_collection(self.trail_lines)
        self.dots = self.flow_ax.scatter([], [], s=3, color="#dbfffd", linewidths=0, zorder=7)
        self.plaque_artists = {}
        self._draw_plaque_artists()
        ys = np.arange(2, self.snapshot.geometry.solid.shape[0]-2, max(2,c.cells_across//9))
        xs = np.arange(3, self.snapshot.geometry.solid.shape[1]-3, max(4,self.snapshot.geometry.solid.shape[1]//44))
        self.qy, self.qx = np.meshgrid(ys,xs,indexing="ij")
        self.quiver = self.lower_ax.quiver(g.x[self.qx]*1000,g.y[self.qy]*1000,
                                           np.zeros_like(self.qx,dtype=float),np.zeros_like(self.qy,dtype=float),
                                           color="#17263a",angles="xy",scale_units="xy",scale=.8,
                                           width=.0019,headwidth=3.5,zorder=3)

    def _draw_plaque_artists(self, preview=None):
        for artists in getattr(self, "plaque_artists", {}).values():
            for artist in artists:
                self._remove_artist(artist)
        self.plaque_artists = {}
        g, c = self.geometry, self.snapshot.config
        plaques = [dict(p) for p in g.plaques]
        if preview is not None:
            plaques = [preview if p["id"] == preview["id"] else p for p in plaques]
        x_mm = g.x*1000
        for p in plaques:
            x0=p["center_mm"]-p["length_mm"]*p["asymmetry"]
            x1=p["center_mm"]+p["length_mm"]*(1-p["asymmetry"])
            xx=np.linspace(x0,x1,80)
            height=c.diameter_mm*p["severity_percent"]/100*plaque_profile(xx,p)
            color=SELECT if p["id"]==self.selected_id else ORANGE
            lw=2.3 if p["id"]==self.selected_id else 1.1
            artists=[]
            if p["wall"] in ("bottom","both"):
                yy=-c.diameter_mm/2+height/(2 if p["wall"]=="both" else 1)
                for ax in (self.flow_ax,self.lower_ax):artists.append(ax.plot(xx,yy,color=color,lw=lw,ls="--" if p["id"]==self.selected_id else "-",zorder=8)[0])
            if p["wall"] in ("top","both"):
                yy=c.diameter_mm/2-height/(2 if p["wall"]=="both" else 1)
                for ax in (self.flow_ax,self.lower_ax):artists.append(ax.plot(xx,yy,color=color,lw=lw,ls="--" if p["id"]==self.selected_id else "-",zorder=8)[0])
            # Selection handles: endpoints are width handles, center marker is height/move handle.
            if p["id"]==self.selected_id:
                for wall_sign in ((-1,) if p["wall"]=="bottom" else (1,) if p["wall"]=="top" else (-1,1)):
                    peak_y=wall_sign*(c.diameter_mm/2-(c.diameter_mm*p["severity_percent"]/100)/(2 if p["wall"]=="both" else 1))
                    x0h=p["center_mm"]-p["length_mm"]*p["asymmetry"]
                    x1h=p["center_mm"]+p["length_mm"]*(1-p["asymmetry"])
                    artists.extend(self.flow_ax.plot([x0h,x1h],[peak_y,peak_y],marker="s",ls="",ms=5,color=SELECT,zorder=9))
                    artists.extend(self.flow_ax.plot([p["center_mm"]],[peak_y],marker="o",ls="",ms=6,mec=TEXT,mfc=SELECT,zorder=9))
            self.plaque_artists[p["id"]]=artists

    def _install_snapshot(self, snapshot):
        changed = snapshot.epoch != self.snapshot.epoch or snapshot.geometry.solid.shape != self.snapshot.geometry.solid.shape
        self.snapshot = snapshot
        self.geometry = snapshot.geometry
        self.config = deepcopy(snapshot.config)
        self.particles_enabled = bool(getattr(snapshot, "particles_enabled", self.config.particles_visible))
        self.pending_heart_rate = float(self.config.heart_rate)
        if hasattr(self, "heart_rate_slider") and abs(self.heart_rate_slider.val - self.config.heart_rate) > 1e-6:
            self.heart_rate_slider.set_val(np.clip(self.config.heart_rate, self.heart_rate_slider.valmin,
                                                   self.heart_rate_slider.valmax))
        self.button_by_name.get("particles", None).label.set_text(
            "Particles: ON" if self.particles_enabled else "Particles: OFF")
        self.button_by_name.get("particle_model", None).label.set_text(
            f"Model: {self.config.particle_model}")
        self.button_by_name.get("pulse", None).label.set_text(
            f"Pulse: {self.config.pulse_shape}")
        if changed:
            if self.selected_id not in {p["id"] for p in self.geometry.plaques}:
                self.selected_id=None
            self._setup_maps()
            self.history.clear()
            self.last_history_time=-np.inf
            self.wall_start=perf_counter()
            self.pressure_ax.set_ylim(0,1)
            self.wss_ax.set_ylim(0,1)
            self.last_trail_key=None
        self._draw_plaque_artists(self.preview_plaque)
        if self.renderer is not None:
            self.renderer.sync()

    def _display(self):
        state,error,frames=self.worker.status()
        if error:
            self.message=error
        if not frames:
            return self.snapshot,self.snapshot,self.snapshot.time, self.snapshot.positions
        item=interpolate_frames(frames,perf_counter(),paused=self.paused,
                                delay=.065*float(self.playback.val))
        if item is None:return self.snapshot,self.snapshot,self.snapshot.time,self.snapshot.positions
        a,b,alpha,positions,display_time=item
        if b.epoch != self.snapshot.epoch:
            self._install_snapshot(b)
        return a,b,display_time,positions

    def render(self):
        a,b,display_time,positions=self._display()
        # A reset can make a/b differ from self.snapshot after installation; use the
        # newest geometry associated with the displayed endpoint.
        g=b.geometry
        fields={}
        alpha=float(np.clip((display_time-a.time)/max(b.time-a.time,1e-9),0,1))
        needed={"u","v","mean_p","wss_low","wss_high"}
        if self.lower_mode=="vorticity":
            needed.add("vorticity")
        for key,value in b.fields.items():
            if key in needed and a is not b and a.epoch==b.epoch and a.fields[key].shape==value.shape:
                fields[key]=(1-alpha)*a.fields[key]+alpha*value
            elif key in ("delta_p","backflow_percent") and a is not b and a.epoch==b.epoch:
                fields[key]=(1-alpha)*a.fields[key]+alpha*value
            else:
                fields[key]=value
        speed=np.hypot(fields["u"],fields["v"])
        self.speed_map.set_data(np.ma.masked_where(g.solid,speed))
        lower=fields["vorticity"] if self.lower_mode=="vorticity" else fields["u"]
        self.lower_map.set_data(np.ma.masked_where(g.solid,lower))
        c=b.config
        # Drawing every simulated tracer is unnecessary for a readable map and
        # is the dominant Matplotlib cost on a fast desktop.  The worker still
        # advects the full bounded set for exports; the display uses a stable
        # deterministic subset and can be disabled entirely.
        show_particles = bool(getattr(b, "particles_enabled", self.particles_enabled))
        if show_particles and positions.size:
            limit = int(c.particle_display_limit)
            if limit and len(positions) > limit:
                display_ix = np.linspace(0, len(positions) - 1, limit, dtype=int)
            else:
                display_ix = np.arange(len(positions))
            display_positions = positions[display_ix]
            display_trails = b.trails[display_ix] if b.trails.ndim == 3 else b.trails
            self.dots.set_offsets(display_positions)
            trail_key=(b.epoch,b.sequence,c.particle_diameter_um,len(display_positions))
            if trail_key != self.last_trail_key:
                self.dots.set_sizes([self._particle_marker_size(c.particle_diameter_um)])
                self.trail_lines.set_segments(display_trails)
                self.last_trail_key=trail_key
        else:
            self.dots.set_offsets(np.empty((0, 2)))
            self.dots.set_sizes(np.empty((0,)))
            self.trail_lines.set_segments(np.empty((0, 0, 2)))
            self.last_trail_key=None
        qu=np.ma.masked_where(g.solid[self.qy,self.qx],fields["u"][self.qy,self.qx])
        qv=np.ma.masked_where(g.solid[self.qy,self.qx],fields["v"][self.qy,self.qx])
        self.quiver.set_UVC(qu,qv)
        x=g.x*1000
        self.pressure_line.set_data(x,fields["mean_p"])
        self.wss_low.set_data(x,np.abs(fields["wss_low"]))
        self.wss_high.set_data(x,np.abs(fields["wss_high"]))
        period=60/c.heart_rate
        pulse_t=np.linspace(0,100,280)
        pulse_y=c.mean_velocity*pulse_curve(c,pulse_t/100*period)
        current=c.mean_velocity*float(pulse_curve(c,[display_time])[0])
        self.pulse_line.set_data(pulse_t,pulse_y)
        self.pulse_marker.set_data([(display_time/period%1)*100],[current])
        self.pulse_phase.set_text(f"{phase_label(c,display_time)}\n{current:.3f} m/s")
        self._expand_limits(self.pressure_ax,fields["mean_p"])
        self._expand_limits(self.wss_ax,np.r_[np.abs(fields["wss_low"]),np.abs(fields["wss_high"])],True)
        for axis,limits in ((self.pressure_ax,(0,x[-1])),(self.wss_ax,(0,x[-1])),(self.pulse_ax,(0,100))):
            if tuple(axis.get_xlim()) != limits:
                axis.set_xlim(*limits)
                if self.renderer is not None:self.renderer.invalidate()
        pulse_limits=(max(0,c.mean_velocity*(1-c.pulsatility_percent/100)*.88),c.mean_velocity*(1+c.pulsatility_percent/100)*1.12)
        if tuple(self.pulse_ax.get_ylim()) != pulse_limits:
            self.pulse_ax.set_ylim(*pulse_limits)
            if self.renderer is not None:self.renderer.invalidate()
        state,error,_=self.worker.status()
        running="RUNNING" if state=="running" and not self.paused else "PAUSED" if self.paused else state.upper()
        self.status.set_text(f"{running}   |   t = {display_time:.3f} s   |   {phase_label(c,display_time)}   |   "
                             f"Re = {c.reynolds:.0f}   |   {c.heart_rate:g} BPM (fixed/run)   |   "
                             f"pulse {current:.3f} m/s   |   seed {c.seed}   |   "
                             f"CFD Δp {fields['delta_p']:.1f} Pa   |   backflow {fields['backflow_percent']:.1f}%")
        pipe_dp=calculate_pressure_drop(c.viscosity,c.length,c.mean_velocity,c.diameter)
        pipe_wss=calculate_wall_shear_stress(c.viscosity,c.mean_velocity,c.diameter)
        womersley = fields.get("womersley", calculate_womersley_number(c.density, c.heart_rate, c.diameter, c.viscosity))
        stokes_number = fields.get("particle_stokes_number", 0.0)
        self.reference.set_text(f"Healthy circular-pipe references: Δp {pipe_dp:.2f} Pa / wall shear {pipe_wss:.3f} Pa"
                                f"     •     Womersley α≈{womersley:.2f}     •     particles {c.particle_model} {c.particle_diameter_um:g} µm"
                                f" (Stₚ≈{stokes_number:.3g}; {'on' if show_particles else 'off'})     •     Live CFD is planar; plaque-wall shear is approximate.")
        if self.selected_id is None:
            plaque_text="No plaque selected. Click a highlighted plaque, then drag its body, square width handles, or circle height handle."
        else:
            p=next((p for p in g.plaques if p["id"]==self.selected_id),None)
            report = fields["geometry_report"]
            metric = next((m for m in report.get("plaques", []) if m["id"] == self.selected_id), None)
            plaque_text=(f"Selected plaque {self.selected_id+1}: {p['wall']} wall  ·  center {p['center_mm']:.1f} mm  ·  "
                         f"width {p['length_mm']:.1f} mm  ·  requested D {p['severity_percent']:.1f}%  ·  "
                         f"realized D≈{metric['realized_raster_diameter_narrowing_percent']:.1f}%  ·  "
                         f"circular-3D area equivalent≈{metric['equivalent_circular_3d_area_reduction_percent']:.1f}%  ·  "
                         "release mouse to recompute") if p and metric else "Selected plaque is rebuilding."
        self.plaque_info.set_text(plaque_text)
        speed=display_time/max(perf_counter()-self.wall_start,1e-6)
        recording = "ON" if self.recording is not None else "OFF"
        fps_text=f"{1/self.frame_interval:.0f} FPS" if self.frame_interval else "timing display"
        self.footer.set_text(f"{self.message}\nDrag plaques to edit · Space pause · R layout · V map · P particles · S snapshot · J settings · "
                             f"Record {recording} · {speed:.3f} sim s / wall s · {b.threads} CPU threads · {fps_text}")


    def tick(self,_frame=0):
        now=perf_counter()
        if self.last_draw_time is not None:
            interval=now-self.last_draw_time
            self.frame_interval=interval if self.frame_interval is None else .9*self.frame_interval+.1*interval
        self.last_draw_time=now
        self.render()
        if self.renderer is not None:self.renderer.draw()
        if self.animation is not None:
            # Timer backends wait after the callback: subtract drawing time
            # so fast machines can reach the requested frame rate.
            self.animation.interval=max(1,round(1000/self.config.viewer_fps-1000*(perf_counter()-now)))
        self.frame+=1
        return []

    def toggle_pause(self,event=None):
        self.paused=not self.paused
        self.worker.set_running(not self.paused)
        self.button_by_name["pause"].label.set_text("Resume" if self.paused else "Pause")
        self.message="Paused; the latest fluid state is held." if self.paused else "Resumed continuous physics."
        self.render();self.fig.canvas.draw_idle()

    def toggle_particles(self,event=None):
        self.particles_enabled=not self.particles_enabled
        self.config=replace(self.config, particles_visible=self.particles_enabled)
        self.worker.set_particles_enabled(self.particles_enabled)
        self.button_by_name["particles"].label.set_text(
            "Particles: ON" if self.particles_enabled else "Particles: OFF")
        self.message=("Particle advection/rendering enabled; heatmaps continue from the same CFD state."
                      if self.particles_enabled else
                      "Particles hidden and advection paused; heatmaps continue without particle drawing cost.")
        self.render();self.fig.canvas.draw_idle()

    def toggle_recording(self,event=None):
        if self.recording is None:
            folder=Path(__file__).resolve().parent / "exports" / "runs"
            try:
                self.recording=RunRecorder(folder, every_s=.05)
                self.worker.set_recorder(self.recording)
                self.button_by_name["record"].label.set_text("Record: ON")
                self.message=f"Recording compact observables to {self.recording.csv_path.name}."
            except (OSError, ValueError) as exc:
                self.recording=None
                self.message=f"Could not start recording: {exc}"
        else:
            recording=self.recording
            self.recording=None
            self.worker.set_recorder(None)
            recording.close()
            self.button_by_name["record"].label.set_text("Record: OFF")
            self.message="Recording stopped; the CSV, metadata, and epoch geometry files are complete."
        self.render();self.fig.canvas.draw_idle()

    def _preview_heart_rate(self,value):
        self.pending_heart_rate=float(value)
        if abs(self.pending_heart_rate-self.config.heart_rate)>1e-6:
            self.message=f"BPM preview: {self.pending_heart_rate:.0f}. Release the slider to rebuild the flow at this fixed rate."

    def _commit_heart_rate(self,value):
        value=float(np.clip(value, self.heart_rate_slider.valmin, self.heart_rate_slider.valmax))
        if abs(value-self.config.heart_rate)<1e-6:
            return
        self.pending_heart_rate=value
        self._submit_config(replace(self.config,heart_rate=value),self.selected_id)
        self.message=f"Heart rate fixed at {value:.0f} BPM for this run; pulsatility remains an amplitude, not BPM variation."

    def adjust_heart_rate(self,delta):
        value=np.clip(round(self.config.heart_rate+delta),self.heart_rate_slider.valmin,self.heart_rate_slider.valmax)
        self.heart_rate_slider.set_val(value)
        self._commit_heart_rate(value)

    def _update_particle_parameters(self,model=None,diameter_um=None,density_kg_m3=None):
        c=replace(self.config,
                  particle_model=model or self.config.particle_model,
                  particle_diameter_um=(self.config.particle_diameter_um if diameter_um is None else diameter_um),
                  particle_density_kg_m3=(self.config.particle_density_kg_m3 if density_kg_m3 is None else density_kg_m3)).validate()
        self.config=deepcopy(c)
        self.worker.set_particle_parameters(c.particle_model,c.particle_diameter_um,c.particle_density_kg_m3)
        self.button_by_name["particle_model"].label.set_text(f"Model: {c.particle_model}")
        self.message=(f"Particles: {c.particle_model}, {c.particle_diameter_um:g} µm, "
                      f"{c.particle_density_kg_m3:g} kg/m³; one-way properties do not alter CFD.")

    def adjust_particle_size(self,delta):
        value=np.clip(self.config.particle_diameter_um+delta,.05,5000.)
        self._update_particle_parameters(diameter_um=value)

    def toggle_particle_model(self,event=None):
        model="inertial" if self.config.particle_model=="tracer" else "tracer"
        self._update_particle_parameters(model=model)

    def _current_plaques(self):
        return [dict(p) for p in self.geometry.plaques]

    def _submit_config(self,config,selected=None):
        try:
            config.validate();make_geometry(config)
        except (ValueError,TypeError) as exc:
            self.message=str(exc);self.preview_plaque=None;self._draw_plaque_artists();self.render();self.fig.canvas.draw_idle();return False
        self.config=deepcopy(config)
        self.selected_id=selected
        self.preview_plaque=None
        self.message="Geometry committed. Rebuilding the flow field while the viewer remains responsive."
        self.worker.reset(self.config)
        return True

    def _custom_config(self,plaques):
        return replace(self.config,geometry_mode="custom",plaques=[dict(p) for p in plaques],
                       plaque_count=len(plaques),randomize_count=False)

    def regenerate(self,event=None):
        seed=(self.config.seed+1)%(2**32)
        c=replace(self.config,geometry_mode="random",plaques=None,seed=seed)
        self._submit_config(c,None)

    def _try_add_candidate(self):
        current=self._current_plaques()
        c=self.config
        cell=c.diameter_mm/c.cells_across
        width=min(max(6*cell,.45*c.stenosis_length_mm),.16*c.length_mm)
        severity=min(max(18,.60*c.stenosis_percent),55)
        for center in np.linspace(.18*c.length_mm,.78*c.length_mm,17):
            for wall in ("bottom","top"):
                p=dict(id=max([q["id"] for q in current],default=-1)+1,center_mm=float(center),length_mm=float(width),
                       severity_percent=float(severity),wall=wall,asymmetry=.55,shape_exponent=1.15)
                try:
                    trial=current+[p]
                    make_geometry(self._custom_config(trial))
                    return trial,p["id"]
                except ValueError:
                    continue
        return None,None

    def add_plaque(self,event=None):
        plaques,selected=self._try_add_candidate()
        if plaques is None:
            self.message="There is no valid open space for another plaque at the current size. Delete, shorten, or move one first."
            self.render();self.fig.canvas.draw_idle();return
        self._submit_config(self._custom_config(plaques),selected)

    def delete_selected(self,event=None):
        if self.selected_id is None:
            self.message="Select a plaque first."
        else:
            plaques=[p for p in self._current_plaques() if p["id"]!=self.selected_id]
            self._submit_config(self._custom_config(plaques),None)
        self.render();self.fig.canvas.draw_idle()

    def _selected_plaque(self):
        return next((dict(p) for p in self.geometry.plaques if p["id"]==self.selected_id),None)

    def adjust_selected(self,kind,delta):
        p=self._selected_plaque()
        if p is None:
            self.message="Select a plaque first."
            self.render();self.fig.canvas.draw_idle();return
        if kind=="width":
            p["length_mm"]+=delta
        else:
            p["severity_percent"]+=delta
        plaques=[p if q["id"]==self.selected_id else dict(q) for q in self._current_plaques()]
        self._submit_config(self._custom_config(plaques),self.selected_id)
        self.render();self.fig.canvas.draw_idle()

    def toggle_pulse_shape(self,event=None):
        shape="sine" if self.config.pulse_shape=="systolic" else "systolic"
        self.config=replace(self.config,pulse_shape=shape)
        self.worker.reset(self.config)
        self.button_by_name["pulse"].label.set_text(f"Pulse: {shape}")
        self.message=f"Inlet waveform set to {shape}; velocity still follows simulated systole/diastole."

    def toggle_map(self,event=None):
        if self.lower_mode=="vorticity":
            self.lower_mode="backflow"
            self.lower_map.set_clim(-self.speed_limit,self.speed_limit)
            self.lower_bar.set_label("Axial velocity (m/s)",color=MUTED,fontsize=8)
            self.lower_ax.set_title("AXIAL VELOCITY  ·  blue = backward, red = forward",loc="left",color=TEXT,fontsize=9,pad=10)
            self.button_by_name["map"].label.set_text("Map: axial")
        else:
            self.lower_mode="vorticity"
            self.lower_map.set_clim(-self.vort_limit,self.vort_limit)
            self.lower_bar.set_label("Vorticity (1/s)",color=MUTED,fontsize=8)
            self.lower_ax.set_title("VORTICITY  ·  signed local rotation (shear also contributes)",loc="left",color=TEXT,fontsize=9,pad=10)
            self.button_by_name["map"].label.set_text("Map: vorticity")
        self.render();self.fig.canvas.draw_idle()

    def _hit_plaque(self,x,y):
        if x is None:return None,"move"
        candidates=[]
        c=self.snapshot.config
        for p in self.geometry.plaques:
            x0=p["center_mm"]-p["length_mm"]*p["asymmetry"]
            x1=p["center_mm"]+p["length_mm"]*(1-p["asymmetry"])
            if x<x0-.9 or x>x1+.9:continue
            h=c.diameter_mm*p["severity_percent"]/100*float(plaque_profile(np.array([x]),p)[0])
            walls=(-1,1) if p["wall"]=="both" else (-1,) if p["wall"]=="bottom" else (1,)
            dist=min(abs(y-(w*(c.diameter_mm/2-h/(2 if p["wall"]=="both" else 1)))) for w in walls)
            candidates.append((dist+abs(x-p["center_mm"])*.08,p))
        if not candidates:return None,"move"
        candidates.sort(key=lambda z:z[0]);p=candidates[0][1]
        x0=p["center_mm"]-p["length_mm"]*p["asymmetry"]
        x1=p["center_mm"]+p["length_mm"]*(1-p["asymmetry"])
        handle=max(.55,c.diameter_mm/c.cells_across*6)
        if abs(x-p["center_mm"])<handle:
            height=c.diameter_mm*p["severity_percent"]/100
            walls=(-1,1) if p["wall"]=="both" else (-1,) if p["wall"]=="bottom" else (1,)
            if any(abs(y-(w*(c.diameter_mm/2-height/(2 if p["wall"]=="both" else 1))))<handle for w in walls):
                return p,"height"
        if abs(x-x0)<max(.7,c.diameter_mm/c.cells_across*6):return p,"left"
        if abs(x-x1)<max(.7,c.diameter_mm/c.cells_across*6):return p,"right"
        return p,"move"

    def on_press(self,event):
        if event.button!=1 or event.inaxes not in (self.flow_ax,self.lower_ax) or event.xdata is None:return
        p,mode=self._hit_plaque(event.xdata,event.ydata)
        if p is not None:
            self.selected_id=p["id"]
            self.dragging=True;self.drag_moved=False;self.drag_mode=mode
            self.drag_origin=(event.xdata,event.ydata);self.drag_base=dict(p);self.preview_plaque=dict(p)
            self._draw_plaque_artists(self.preview_plaque);self.render();self.fig.canvas.draw_idle();return
        if self.particles_enabled:
            self.worker.inject(event.xdata,event.ydata)
            self.message=f"Injected tracers near ({event.xdata:.1f}, {event.ydata:.2f}) mm."
        else:
            self.message="Particles are off; turn them on before injecting tracers."

    def on_motion(self,event):
        if not self.dragging or event.inaxes not in (self.flow_ax,self.lower_ax) or event.xdata is None:return
        p=dict(self.drag_base);x=float(event.xdata);c=self.snapshot.config;cell=c.diameter_mm/c.cells_across
        xmargin=2*cell
        if self.drag_mode=="move":
            dx=x-self.drag_origin[0]
            x0=p["center_mm"]-p["length_mm"]*p["asymmetry"]+dx
            x1=p["center_mm"]+p["length_mm"]*(1-p["asymmetry"])+dx
            if x0<xmargin:dx+=xmargin-x0
            if x1>c.length_mm-xmargin:dx-=x1-(c.length_mm-xmargin)
            p["center_mm"]+=dx
            if p["wall"] in ("bottom","top") and event.ydata is not None and abs(event.ydata)>.18*c.diameter_mm:
                p["wall"]="top" if event.ydata>0 else "bottom"
        elif self.drag_mode=="height":
            dy=float(event.ydata)-self.drag_origin[1]
            if p["wall"]=="bottom":
                change=dy*100/c.diameter_mm
            elif p["wall"]=="top":
                change=-dy*100/c.diameter_mm
            else:
                change=abs(dy)*100/c.diameter_mm
            p["severity_percent"]=float(np.clip(p["severity_percent"]+change,5,95))
        else:
            left=p["center_mm"]-p["length_mm"]*p["asymmetry"]
            right=p["center_mm"]+p["length_mm"]*(1-p["asymmetry"])
            minimum=6*cell
            if self.drag_mode=="left":left=min(x,right-minimum)
            else:right=max(x,left+minimum)
            left=max(xmargin,left);right=min(c.length_mm-xmargin,right)
            p["length_mm"]=right-left;p["center_mm"]=left+p["length_mm"]*p["asymmetry"]
        self.drag_moved=True;self.preview_plaque=p
        self._draw_plaque_artists(p);self.message=f"Previewing plaque {p['id']+1}; release to commit and recompute."
        self.render();self.fig.canvas.draw_idle()

    def on_release(self,event):
        slider = getattr(self, "heart_rate_slider", None)
        if slider is not None and event.inaxes is slider.ax:
            self._commit_heart_rate(self.pending_heart_rate)
            self.render();self.fig.canvas.draw_idle()
            return
        if not self.dragging:return
        self.dragging=False
        p=self.preview_plaque
        self.preview_plaque=None
        if self.drag_moved and p is not None:
            plaques=[p if q["id"]==p["id"] else dict(q) for q in self._current_plaques()]
            self._submit_config(self._custom_config(plaques),p["id"])
        else:
            self._draw_plaque_artists();self.message="Plaque selected. Drag its body, square width handles, or circle height handle to edit it."
        self.render();self.fig.canvas.draw_idle()

    def on_key(self,event):
        actions={" ":self.toggle_pause,"r":self.regenerate,"v":self.toggle_map,"s":self.save_snapshot,"j":self.save_settings,
                 "[":lambda:self.adjust_selected("width",-.5),"]":lambda:self.adjust_selected("width",.5),
                 "{":lambda:self.adjust_selected("height",-5),"}":lambda:self.adjust_selected("height",5),
                 "a":self.add_plaque,"delete":self.delete_selected,"p":self.toggle_particles,
                 "b":lambda:self.adjust_heart_rate(1),"n":lambda:self.adjust_heart_rate(-1)}
        if event.key in actions:actions[event.key]()

    def save_snapshot(self,event=None):
        # Use the latest published fields, so exports contain a complete solver state.
        snap=self.worker.latest() or self.snapshot
        d=snap.fields;g=snap.geometry;c=snap.config
        folder=Path(__file__).resolve().parent/"exports";folder.mkdir(exist_ok=True)
        stem=f"hemoflow_seed{c.seed}_{datetime.now():%Y%m%d_%H%M%S_%f}"
        self.renderer.savefig(folder/f"{stem}.png",dpi=150,facecolor=BG)
        np.savez_compressed(folder/f"{stem}.npz",x_m=g.x,y_m=g.y,solid=g.solid,u_m_s=d["u"],v_m_s=d["v"],pressure_pa=d["p"],vorticity_s=d["vorticity"],time_s=snap.time,dt_s=snap.dt,iteration=snap.iteration,config_json=json.dumps(asdict(c)),geometry_report_json=json.dumps(d.get("geometry_report",stenosis_report(c,g))))
        np.savetxt(folder/f"{stem}.csv",np.column_stack((g.x*1000,d["mean_u"],d["mean_p"],d["wss_low"],d["wss_high"])),delimiter=",",header="x_mm,mean_axial_velocity_m_s,mean_pressure_pa,lower_wss_estimate_pa,upper_wss_estimate_pa",comments="")
        (folder/f"{stem}_stenosis_report.json").write_text(json.dumps(d.get("geometry_report",stenosis_report(c,g)),indent=2)+"\n",encoding="utf-8")
        reproducible=replace(c,geometry_mode="custom",plaques=[dict(p) for p in g.plaques],plaque_count=len(g.plaques),randomize_count=False)
        reproducible.save(folder/f"{stem}.json")
        self.message=f"Saved dashboard, fields, profiles, and reproducible plaque settings to exports/{stem}.*"
        self.render();self.fig.canvas.draw_idle()

    def save_settings(self,event=None):
        self.config.save(Path(__file__).resolve().parent/"settings.json")
        self.message="Saved current inputs and plaque layout settings."
        self.render();self.fig.canvas.draw_idle()

    def on_close(self,event=None):
        if self.recording is not None:
            recording=self.recording
            self.recording=None
            self.worker.set_recorder(None)
            recording.close()
        self.worker.stop()
        self.worker.join(timeout=.8)
        if self.animation is not None:self.animation.stop()

    def show(self):
        self.animation=self.fig.canvas.new_timer(interval=max(1,round(1000/self.config.viewer_fps)))
        self.animation.add_callback(self.tick)
        self.animation.start()
        plt.show()
