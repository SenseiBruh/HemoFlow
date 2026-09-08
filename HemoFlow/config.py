"""Physical inputs are SI except the explicitly named mm / percent fields."""
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path


@dataclass
class Config:
    density: float = 1060.0
    viscosity: float = 0.0035
    diameter_mm: float = 4.0
    length_mm: float = 60.0
    mean_velocity: float = 0.30
    stenosis_percent: float = 50.0
    stenosis_center_mm: float = 30.0
    stenosis_length_mm: float = 15.0
    heart_rate: float = 72.0
    pulsatility_percent: float = 40.0
    plaque_count: int = 3
    seed: int = 42
    geometry_mode: str = "random"
    cells_across: int = 40
    particle_count: int = 1000
    # Particle properties are one-way probes: they do not alter the blood
    # solution.  ``tracer`` follows the solved velocity exactly; ``inertial``
    # adds a Stokes relaxation time so size/density sweeps are meaningful.
    particle_model: str = "tracer"
    particle_diameter_um: float = 8.0
    particle_density_kg_m3: float = 1060.0
    particles_visible: bool = True
    particle_display_limit: int = 650
    pulse_shape: str = "systolic"
    randomize_count: bool = True
    # The default layout stays readable: one or two plaques. Three is still
    # available as an explicit advanced setting.
    max_random_plaques: int = 2
    randomize_on_launch: bool = True
    plaques: list | None = None
    compute_threads: int = 0

    @property
    def diameter(self):
        return self.diameter_mm / 1000.0

    @property
    def length(self):
        return self.length_mm / 1000.0

    @property
    def reynolds(self):
        return self.density * self.mean_velocity * self.diameter / self.viscosity

    def validate(self):
        for name in ("density", "viscosity", "diameter_mm", "length_mm", "mean_velocity"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a finite positive number.")
        if not math.isfinite(self.heart_rate) or not 20 <= self.heart_rate <= 240:
            raise ValueError("heart_rate must be between 20 and 240 BPM.")
        for name, lo, hi in (("particle_diameter_um", .05, 5000.0),
                             ("particle_density_kg_m3", 50.0, 25_000.0)):
            value = getattr(self, name)
            if not math.isfinite(value) or not lo <= value <= hi:
                raise ValueError(f"{name} must be between {lo:g} and {hi:g}.")
        for name, lo, hi in (("stenosis_percent", 0, 95), ("pulsatility_percent", 0, 95),
                             ("stenosis_center_mm", 0, self.length_mm),
                             ("stenosis_length_mm", 0, self.length_mm)):
            value = getattr(self, name)
            if not math.isfinite(value) or not lo <= value <= hi:
                raise ValueError(f"{name} must be between {lo} and {hi}.")
        for name, lo, hi in (("plaque_count", 0, 12), ("cells_across", 24, 256),
                             ("particle_count", 100, 5000), ("particle_display_limit", 0, 5000),
                             ("seed", 0, 2**32 - 1)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not lo <= value <= hi:
                raise ValueError(f"{name} must be an integer between {lo} and {hi}.")
        if self.geometry_mode not in ("random", "single", "healthy", "custom"):
            raise ValueError("geometry_mode must be random, single, healthy, or custom.")
        if self.pulse_shape not in ("sine", "systolic"):
            raise ValueError("pulse_shape must be sine or systolic.")
        if self.particle_model not in ("tracer", "inertial"):
            raise ValueError("particle_model must be tracer or inertial.")
        if not isinstance(self.particles_visible, bool):
            raise ValueError("particles_visible must be true or false.")
        if not isinstance(self.randomize_count, bool) or not isinstance(self.randomize_on_launch, bool):
            raise ValueError("Randomization switches must be true or false.")
        if type(self.max_random_plaques) is not int or not 1 <= self.max_random_plaques <= 3:
            raise ValueError("max_random_plaques must be 1, 2, or 3.")
        if type(self.compute_threads) is not int or not 0 <= self.compute_threads <= 32:
            raise ValueError("compute_threads must be an integer from 0 (auto) to 32.")
        if self.geometry_mode == "custom" and (not isinstance(self.plaques, list) or len(self.plaques) > 12):
            raise ValueError("Custom geometry needs a list of up to 12 plaques.")
        if self.length < 4 * self.diameter:
            raise ValueError("Use a vessel length of at least four diameters for this open-channel model.")
        if self.cells_across**2 * self.length / self.diameter > 600_000:
            raise ValueError("This grid is too large for the desktop preview. Reduce length or cells_across.")
        return self

    def save(self, filename):
        Path(filename).write_text(json.dumps(asdict(self), indent=2) + "\n", encoding="utf-8")

    @classmethod
    def load(cls, filename):
        values = json.loads(Path(filename).read_text(encoding="utf-8"))
        # Retain the original pulse if an earlier settings file is explicitly loaded.
        values.setdefault("pulse_shape", "sine")
        values.setdefault("particle_model", "tracer")
        values.setdefault("particle_diameter_um", 8.0)
        values.setdefault("particle_density_kg_m3", 1060.0)
        values.setdefault("particles_visible", True)
        values.setdefault("particle_display_limit", 650)
        return cls(**values).validate()
