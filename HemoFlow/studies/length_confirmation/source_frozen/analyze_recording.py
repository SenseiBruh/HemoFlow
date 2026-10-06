"""Time-weighted cycle summaries; marks partial cycles and exports Excel-ready CSV."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from recorder import EXTRA_METRICS

NUMERIC = (
    "time_s", "pulse_factor", "inlet_speed_m_s", "mean_inlet_velocity_m_s",
    "mean_outlet_velocity_m_s", "inlet_flow_per_depth_m2_s", "outlet_flow_per_depth_m2_s",
    "delta_p_pa", "max_velocity_m_s", "min_axial_velocity_m_s", "backflow_percent",
    "max_wss_pa", "mean_abs_wss_pa", "max_vorticity_s_1", "flux_error", "mach",
    "density_variation", "particles_recycled", "particle_relaxation_time_s",
    "particle_stokes_number",
) + tuple(name for name in EXTRA_METRICS if name != "mass_balance_valid")


def summarize(csv_path, discard_cycles=0, complete_only=False):
    csv_path = Path(csv_path)
    if any(ch in str(csv_path) for ch in "*?[]"):
        matches = sorted(csv_path.parent.glob(csv_path.name))
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one CSV for pattern {csv_path!s}; found {len(matches)}")
        csv_path = matches[0]
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("The recording has no data rows")
    groups = {}
    for row in rows:
        groups.setdefault(int(row["epoch"]), []).append(row)
    cycles = []
    for epoch, items in sorted(groups.items()):
        ts = np.asarray([float(r["time_s"]) for r in items])
        if np.any(np.diff(ts) <= 0):
            raise ValueError(f"Epoch {epoch} has duplicate or decreasing times")
        phase = np.asarray([int(r["cycle_index"])+float(r["phase_fraction"]) for r in items])
        valid = phase > 1e-10
        if not valid.any():
            continue
        period = float(np.median(ts[valid]/phase[valid]))
        if not np.isfinite(period) or period <= 0:
            raise ValueError(f"Cannot determine the cardiac period for epoch {epoch}")
        for cycle in range(int(items[0]["cycle_index"]), int(items[-1]["cycle_index"])+1):
            if cycle < discard_cycles:
                continue
            start, stop = cycle*period, (cycle+1)*period
            full = bool(ts[0] <= start+1e-10 and ts[-1] >= stop-1e-10)
            if complete_only and not full:
                continue
            left, right = max(start,ts[0]), min(stop,ts[-1])
            if right <= left:
                continue
            inside = (ts > left) & (ts < right)
            timeline = np.r_[left,ts[inside],right]
            out = dict(epoch=epoch,cycle_index=cycle,complete_cycle=full,
                       period_s=period,coverage_fraction=float((right-left)/period),
                       samples=int(((ts>=left)&(ts<=right)).sum()),
                       statistics_start_s=float(left),statistics_end_s=float(right))
            for name in NUMERIC:
                if name not in items[0]:
                    continue
                values = np.asarray([float(r[name]) for r in items])
                if not np.isfinite(values).all():
                    continue
                y = np.interp(timeline,ts,values)
                weights = np.diff(timeline)
                mean = float(np.sum(weights*(y[:-1]+y[1:])/2)/(right-left))
                # Exact second moment of the piecewise-linear interpolant.
                second = float(np.sum(weights*(y[:-1]**2+y[:-1]*y[1:]+y[1:]**2)/3)/(right-left))
                out.update({f"{name}_mean":mean,f"{name}_std":float(np.sqrt(max(0,second-mean**2))),
                            f"{name}_min":float(y.min()),f"{name}_max":float(y.max())})
            cycles.append(out)
    return dict(source_csv=str(csv_path),rows=len(rows),epochs=sorted(groups),
                statistics_method="Time-weighted piecewise-linear interpolation; boundary samples interpolated. std is temporal population std. Partial cycles are explicitly marked. A complete cycle alone does not establish periodic convergence.",
                discard_cycles=discard_cycles,complete_only=complete_only,cycles=cycles)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv",type=Path)
    parser.add_argument("--output",type=Path,help="JSON summary")
    parser.add_argument("--csv-output",type=Path,help="Flat cycle summary that opens in Excel")
    parser.add_argument("--discard-cycles",type=int,default=0,help="Discard this many initial cycles in every epoch")
    parser.add_argument("--complete-only",action="store_true")
    args = parser.parse_args()
    if args.discard_cycles<0:parser.error("--discard-cycles must be non-negative")
    result = summarize(args.csv,args.discard_cycles,args.complete_only)
    text = json.dumps(result,indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding="utf-8")
    elif not args.csv_output:
        print(text,end="")
    if args.csv_output:
        args.csv_output.parent.mkdir(parents=True,exist_ok=True)
        keys=list(dict.fromkeys(key for row in result["cycles"] for key in row))
        with args.csv_output.open("w",newline="",encoding="utf-8-sig") as handle:
            writer=csv.DictWriter(handle,fieldnames=keys)
            writer.writeheader();writer.writerows(result["cycles"])
        print(f"Saved {len(result['cycles'])} cycle summaries to {args.csv_output}")


if __name__=="__main__":main()
