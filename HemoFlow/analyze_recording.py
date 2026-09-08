"""Summarize a HemoFlow ``*_timeseries.csv`` by epoch and cardiac cycle.

This is intentionally dependency-light so a recorded run can be reduced on a
different machine without reopening the viewer.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


NUMERIC = (
    "time_s", "pulse_factor", "inlet_speed_m_s", "mean_inlet_velocity_m_s",
    "mean_outlet_velocity_m_s", "delta_p_pa", "max_velocity_m_s",
    "min_axial_velocity_m_s", "backflow_percent", "max_wss_pa",
    "mean_abs_wss_pa", "max_vorticity_s_1", "flux_error", "mach",
    "density_variation", "particles_recycled",
    "particle_relaxation_time_s", "particle_stokes_number",
)


def summarize(csv_path):
    csv_path = Path(csv_path)
    if any(ch in str(csv_path) for ch in "*?[]"):
        matches = sorted(csv_path.parent.glob(csv_path.name))
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one CSV for pattern {csv_path!s}; found {len(matches)}")
        csv_path = matches[0]
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    groups = {}
    for row in rows:
        key = (int(row["epoch"]), int(row["cycle_index"]))
        groups.setdefault(key, []).append(row)

    cycles = []
    for (epoch, cycle), items in sorted(groups.items()):
        out = {"epoch": epoch, "cycle_index": cycle, "samples": len(items)}
        for name in NUMERIC:
            if name not in items[0]:
                # v2 recordings predate particle response metrics; keep the
                # summarizer backward-compatible instead of inventing values.
                continue
            values = np.asarray([float(item[name]) for item in items], dtype=float)
            out[f"{name}_mean"] = float(np.mean(values))
            out[f"{name}_std"] = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            out[f"{name}_min"] = float(np.min(values))
            out[f"{name}_max"] = float(np.max(values))
        cycles.append(out)
    return {
        "source_csv": str(csv_path),
        "rows": len(rows),
        "epochs": sorted({int(row["epoch"]) for row in rows}),
        "cycles": cycles,
    }


def main():
    parser = argparse.ArgumentParser(description="Summarize HemoFlow time-series recording by cardiac cycle")
    parser.add_argument("csv", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = summarize(args.csv)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
