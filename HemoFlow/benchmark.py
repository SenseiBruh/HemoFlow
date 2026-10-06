"""Measure this PC's CPU solver throughput without running a full study."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import statistics
import time

from config import Config
from provenance import run_manifest
from solver import FlowSolver


def benchmark(config, grids=(40, 80, 160), steps=120):
    results = []
    for cells in grids:
        solver = FlowSolver(replace(config, cells_across=cells))
        solver.tune_threads()
        solver.step(4)
        trials = []
        for _ in range(3):
            started = time.perf_counter()
            solver.step(steps)
            trials.append(time.perf_counter()-started)
            solver.check_stability()
        seconds_per_step = statistics.median(trials)/steps
        seconds_per_cycle = seconds_per_step/solver.dt*60/config.heart_rate
        row = dict(cells_across=cells, compute_threads=solver.threads_used,
                   time_scale=config.time_scale, dt_s=solver.dt,
                   million_cell_updates_per_s=len(solver.nodes)/seconds_per_step/1e6,
                   estimated_minutes_per_cycle=seconds_per_cycle/60,
                   estimated_hours_for_20_cycles=seconds_per_cycle*20/3600,
                   seconds_per_step=seconds_per_step,
                   thread_trial_seconds=getattr(solver,"thread_timings",{}),
                   manifest=run_manifest(solver))
        results.append(row)
        print(f"{cells:3d} cells | {solver.threads_used:2d} threads | "
              f"{seconds_per_cycle/60:.2f} estimated min/cycle | "
              f"{seconds_per_cycle*20/3600:.2f} estimated hours/20 cycles",flush=True)
    return dict(scope="Short CPU-only throughput estimate after compilation; no GUI, recording, or exports. Not a stability/accuracy experiment. Throat geometry and sustained load change timings.",
                results=results)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells",type=int,nargs="+",default=[40,80,160])
    parser.add_argument("--threads",type=int,default=0)
    parser.add_argument("--time-scale",type=float,default=.5)
    parser.add_argument("--config",type=Path,help="Benchmark a particular plaque geometry")
    parser.add_argument("--output",type=Path,default=Path("experiments/pc_benchmark.json"))
    args=parser.parse_args()
    config=Config.load(args.config) if args.config else Config(geometry_mode="healthy")
    config=replace(config,compute_threads=args.threads,time_scale=args.time_scale).validate()
    result=benchmark(config,args.cells)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(f"Saved {args.output}")


if __name__=="__main__":
    main()
