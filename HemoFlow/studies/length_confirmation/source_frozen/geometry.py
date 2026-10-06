"""Seeded, freely placed wall plaques and an editable explicit geometry."""
from dataclasses import dataclass
import math
import numpy as np


@dataclass
class Geometry:
    x: np.ndarray
    y: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    solid: np.ndarray
    dx: float
    plaques: list

    @property
    def fluid(self):
        return ~self.solid

    @property
    def min_gap_cells(self):
        return int(np.min(np.sum(self.fluid, axis=0)))


def normalize_plaque(plaque, config, index=0):
    p = dict(plaque)
    p.setdefault("id", index)
    p.setdefault("asymmetry", .58)
    p.setdefault("shape_exponent", 1.0)
    p.setdefault("wall", "bottom")
    for name, lo, hi in (("center_mm",0,config.length_mm), ("length_mm",0,config.length_mm),
                         ("severity_percent",0,95), ("asymmetry",.25,.75), ("shape_exponent",.65,2.5)):
        if name not in p or not isinstance(p[name],(int,float)) or not math.isfinite(p[name]) or not lo <= p[name] <= hi:
            raise ValueError(f"Plaque {index+1}: {name} must be between {lo:g} and {hi:g}.")
        p[name] = float(p[name])
    if type(p["id"]) is not int or p["id"] < 0:
        raise ValueError("Plaque IDs must be non-negative integers.")
    if p["wall"] not in ("bottom","top","both"):
        raise ValueError("Plaque wall must be bottom, top, or both.")
    cell_mm = config.diameter_mm / config.cells_across
    if p["length_mm"] < 6*cell_mm:
        raise ValueError(f"Plaques need a length of at least {6*cell_mm:.2f} mm on this grid.")
    a = p["center_mm"] - p["length_mm"]*p["asymmetry"]
    b = p["center_mm"] + p["length_mm"]*(1-p["asymmetry"])
    if a < 2*cell_mm or b > config.length_mm-2*cell_mm:
        raise ValueError("Keep plaques at least two grid cells clear of the inlet and outlet.")
    return p


def plaque_profile(x_mm, plaque):
    left = plaque["length_mm"]*plaque["asymmetry"]
    right = plaque["length_mm"]-left
    z = np.where(x_mm < plaque["center_mm"], (x_mm-plaque["center_mm"])/left,
                 (x_mm-plaque["center_mm"])/right)
    return np.where(np.abs(z)<1, np.maximum(0,.5*(1+np.cos(np.pi*np.clip(z,-1,1))))**plaque["shape_exponent"], 0)


def generate_plaques(config):
    rng = np.random.default_rng(config.seed)
    if config.plaque_count == 0 or config.stenosis_percent == 0 or config.stenosis_length_mm == 0:
        return []
    if config.randomize_count:
        choices = np.arange(1, config.max_random_plaques+1)
        weights = np.array([.32,.48,.20])[:len(choices)]
        if config.stenosis_percent >= 65:
            weights = np.array([.78,.22,0.0])[:len(choices)]
        count = int(rng.choice(choices,p=weights/weights.sum()))
    else:
        count = config.plaque_count
    cell = config.diameter_mm/config.cells_across
    max_length = min(config.stenosis_length_mm, .35*config.length_mm)
    if max_length < 6*cell:
        raise ValueError(f"Set maximum plaque length to at least {6*cell:.2f} mm.")
    a = max(.07*config.length_mm, .7*config.diameter_mm)
    b = config.length_mm-max(.14*config.length_mm,1.5*config.diameter_mm)
    plaques = []
    for index in range(count):
        placed = False
        for attempt in range(180):
            width = rng.uniform(max(6*cell,.42*max_length),max_length)
            if config.stenosis_percent >= 65:
                width = max(width,.72*max_length)
            asym = rng.uniform(.38,.72)
            low, high = a+width*asym, b-width*(1-asym)
            if low >= high:
                continue
            center = rng.uniform(low,high)
            peak_max = min(config.stenosis_percent,100*(1-7.5/config.cells_across))
            severity = rng.uniform(.56,.98)*peak_max
            wall = "bottom" if rng.random()<.5 else "top"
            p = dict(id=index,center_mm=float(center),length_mm=float(width),
                     severity_percent=float(severity),wall=wall,asymmetry=float(asym),
                     shape_exponent=float(rng.uniform(.75,1.9)))
            left,right = center-width*asym, center+width*(1-asym)
            spacing = max(2*cell,.25*config.diameter_mm)
            overlaps = any(not (right+spacing < q["center_mm"]-q["length_mm"]*q["asymmetry"] or
                                    left-spacing > q["center_mm"]+q["length_mm"]*(1-q["asymmetry"]))
                           for q in plaques)
            if not overlaps:
                plaques.append(p)
                placed = True
                break
        if not placed:
            if config.randomize_count:
                break
            raise ValueError("The fixed plaque count does not fit with this length and spacing. Use fewer or shorter plaques.")
    return sorted(plaques,key=lambda p:p["center_mm"])


def resolve_plaques(config):
    if config.geometry_mode == "healthy":
        return []
    if config.geometry_mode == "custom":
        plaques = config.plaques
    elif config.geometry_mode == "single":
        if config.stenosis_percent == 0 or config.stenosis_length_mm == 0:
            return []
        plaques = [dict(id=0,center_mm=config.stenosis_center_mm,length_mm=config.stenosis_length_mm,
                        severity_percent=config.stenosis_percent,wall="both",asymmetry=.5,shape_exponent=1.0)]
    else:
        plaques = generate_plaques(config)
    result = [normalize_plaque(p,config,i) for i,p in enumerate(plaques)]
    if len({p["id"] for p in result}) != len(result):
        raise ValueError("Plaques must have unique IDs.")
    return result


def wall_profiles(config, x_mm, plaques):
    bottom = np.zeros_like(x_mm,dtype=float)
    top = np.zeros_like(x_mm,dtype=float)
    for p in plaques:
        height = config.diameter_mm*p["severity_percent"]/100*plaque_profile(x_mm,p)
        if p["wall"] in ("bottom","both"):
            bottom = np.maximum(bottom,height/(2 if p["wall"]=="both" else 1))
        if p["wall"] in ("top","both"):
            top = np.maximum(top,height/(2 if p["wall"]=="both" else 1))
    return -config.diameter_mm/2+bottom, config.diameter_mm/2-top


def make_geometry(config):
    config.validate()
    dx = config.diameter/config.cells_across
    x = np.arange(int(round(config.length/dx))+1)*dx
    y = (np.arange(config.cells_across+2)-.5)*dx-config.diameter/2
    plaques = resolve_plaques(config)
    lower, upper = wall_profiles(config,x*1000,plaques)
    lower,upper = lower/1000,upper/1000
    solid = (y[:,None]<=lower[None,:]) | (y[:,None]>=upper[None,:])
    geometry = Geometry(x,y,lower,upper,solid,dx,plaques)
    if geometry.min_gap_cells < 6:
        raise ValueError(f"The throat has only {geometry.min_gap_cells} cells. Increase grid resolution or reduce narrowing.")
    if not np.all(np.any(geometry.fluid[:,:-1] & geometry.fluid[:,1:],axis=0)):
        raise ValueError("The plaques disconnect the fluid passage. Move or resize them.")
    return geometry


def stenosis_report(config, geometry):
    """Return reproducible requested-versus-realized stenosis measurements.

    ``severity_percent`` is the requested linear lumen-height/diameter
    narrowing used by this 2D model.  The rasterized grid can round that
    value, so the report also records the realized 2D diameter equivalent and
    the circular-3D area equivalent.  The latter is a comparison label only;
    no 3D cross-sectional area is solved here.
    """
    nominal = float(config.diameter_mm)
    dx_mm = float(geometry.dx * 1000.0)
    gap_mm = (geometry.upper - geometry.lower) * 1000.0
    overall_i = int(np.argmin(gap_mm))
    overall_gap = float(gap_mm[overall_i])
    raster_counts = geometry.fluid.sum(axis=0)
    raster_gaps = raster_counts * dx_mm
    overall_raster_i = int(np.argmin(raster_gaps))
    overall_raster_gap = float(raster_gaps[overall_raster_i])

    def one(plaque):
        left = plaque["center_mm"] - plaque["length_mm"] * plaque["asymmetry"]
        right = plaque["center_mm"] + plaque["length_mm"] * (1 - plaque["asymmetry"])
        support = np.flatnonzero((geometry.x * 1000.0 >= left - dx_mm) &
                                 (geometry.x * 1000.0 <= right + dx_mm))
        if support.size == 0:
            support = np.array([int(np.argmin(np.abs(geometry.x * 1000.0 - plaque["center_mm"])))])
        local_i = int(support[np.argmin(gap_mm[support])])
        local_gap = float(gap_mm[local_i])
        # Count fluid nodes at the throat as the measurement actually used by
        # the bounce-back solver, rather than implying sub-grid precision.
        raster_i = int(support[np.argmin(raster_gaps[support])])
        raster_cells = int(raster_counts[raster_i])
        raster_gap = raster_cells * dx_mm
        d_red = 100.0 * (1.0 - raster_gap / nominal)
        continuous_red = 100.0 * (1.0 - local_gap / nominal)
        return {
            "id": int(plaque["id"]),
            "wall": plaque["wall"],
            "center_mm": float(plaque["center_mm"]),
            "length_mm": float(plaque["length_mm"]),
            "requested_diameter_narrowing_percent": float(plaque["severity_percent"]),
            "realized_continuous_diameter_narrowing_percent": float(continuous_red),
            "realized_raster_diameter_narrowing_percent": float(d_red),
            "equivalent_circular_3d_area_reduction_percent": float(100.0 * (1.0 - (1.0 - d_red / 100.0) ** 2)),
            "minimum_continuous_lumen_mm": local_gap,
            "minimum_raster_lumen_mm": float(raster_gap),
            "throat_x_mm": float(geometry.x[raster_i] * 1000.0),
            "continuous_throat_x_mm": float(geometry.x[local_i] * 1000.0),
            "throat_cells": raster_cells,
            "asymmetry": float(plaque["asymmetry"]),
            "shape_exponent": float(plaque["shape_exponent"]),
        }

    plaques = [one(p) for p in geometry.plaques]
    overall_d = 100.0 * (1.0 - overall_raster_gap / nominal)
    return {
        "model": "2D planar rasterized lumen; area value is a circular-3D equivalent",
        "nominal_diameter_mm": nominal,
        "grid_spacing_mm": dx_mm,
        "overall_minimum_lumen_mm": overall_gap,
        "overall_minimum_continuous_lumen_mm": overall_gap,
        "overall_minimum_raster_lumen_mm": overall_raster_gap,
        "overall_throat_cells": int(raster_counts[overall_raster_i]),
        "overall_raster_throat_x_mm": float(geometry.x[overall_raster_i] * 1000.0),
        "overall_realized_continuous_diameter_narrowing_percent": float(100 * (1 - overall_gap / nominal)),
        "overall_realized_raster_diameter_narrowing_percent": float(overall_d),
        "overall_equivalent_circular_3d_area_reduction_percent": float(100.0 * (1.0 - (1.0 - overall_d / 100.0) ** 2)),
        "plaques": plaques,
    }
