"""Headless wall benchmarks and inlet-symmetry sensitivity experiments.

Case-boundary resume only: interrupted cases restart in new attempt folders.
No prior result or frozen solver file is overwritten.
"""
import argparse
import csv
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import traceback
import zipfile
import numpy as np

from experiment import ROOT,VERSION,plan,make_solver,initialize_pair,advance_to,WallSampler,reflected_populations
from measurements import Probes
from startup_controls import initialize as control_initialize, sample_time, symmetry_row
from audit_math import manufactured_checks
import strict_analysis as analysis


def stamp():return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    t=Path(str(p)+'.tmp');t.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8');t.replace(p)


class CSV:
    def __init__(self,p):self.p=Path(p);self.f=None;self.writer=None
    def row(self,r):
        if self.f is None:
            self.f=self.p.open('w',newline='',encoding='utf-8')
            self.writer=csv.DictWriter(self.f,fieldnames=list(r));self.writer.writeheader()
        self.writer.writerow(r);self.f.flush()
    def close(self):
        if self.f:self.f.close()


def verify_source():
    expected=read(ROOT/'FROZEN_SOURCE.json')
    for n,h in expected.items():
        if sha(ROOT/n)!=h:raise ValueError(f'Frozen source was changed: {n}')
    package={p.relative_to(ROOT).as_posix():sha(p) for p in ROOT.glob('*.py')}
    package.update(expected)
    return package


def preflight(folder,threads):
    folder.mkdir(parents=True,exist_ok=True)
    tests=manufactured_checks()
    write(folder/'analytical_checks.json',tests)
    if not all(r['passed'] for r in tests):raise RuntimeError('Analytical measurement tests failed')
    report=dict(version=VERSION,analytical_passes=len(tests),checks=[],runtime_estimates={})
    import numba
    report['environment']=dict(python=sys.version,numpy=np.__version__,numba=numba.__version__)
    for n in [80,120]:
        spec=next(s for s in plan() if s['kind']=='severe' and s['cells']==n and s['epsilon']>0)
        p=make_solver(spec,1);m=make_solver(spec,1)
        po,pd,pm=initialize_pair(p,.001);mo,md,mm=initialize_pair(m,-.001)
        reflection_error=float(np.max(abs(pd-md[::-1])))
        if reflection_error!=0:raise RuntimeError('Inlet pair is not an exact reflection')
        cutoff=math.ceil((60/72)/p.dt)
        advance_to(p,64,po,pd,cutoff);advance_to(m,64,mo,md,cutoff)
        mirrored=reflected_populations(p.f,p.ny,p.nx)
        error=float(np.linalg.norm((m.f-mirrored)[m.nodes])/np.linalg.norm(mirrored[m.nodes]))
        if error>1e-9:raise RuntimeError(f'Short reflection check failed at N{n}: {error}')
        # Check zero-input wrapper against the unmodified stepping path.
        a=make_solver(spec,1);b=make_solver(spec,1)
        original,driven,_=initialize_pair(a,0.)
        advance_to(a,16,original,driven,0);b.step(16)
        if not np.array_equal(a.f[a.nodes],b.f[b.nodes]):raise RuntimeError('Zero-input wrapper changed the frozen solver')
        # Benchmark a fresh solver; trial state is never reused in a science case.
        s=make_solver(spec,1)
        candidates=[min(threads,numba.config.NUMBA_NUM_THREADS)] if threads else [v for v in [1,2,4,8,16] if v<=numba.config.NUMBA_NUM_THREADS]
        timings={}
        for count in candidates:
            s.set_threads(count);s.step(16)  # JIT compilation excluded
            samples=[]
            for _ in range(3):
                start=time.perf_counter();s.step(128);samples.append((time.perf_counter()-start)/128)
            timings[count]=float(np.median(samples))
        best=min(timings,key=timings.get)
        report['runtime_estimates'][str(n)]=dict(threads=best,seconds_per_step=timings[best],
            estimated_solver_hours_per_10_cycles=timings[best]*(10*60/72)/spec['dt_s']/3600,
            caveat='Short-kernel estimate excludes recording, file compression and thermal throttling.',thread_trials=timings)
        report['checks'].append(dict(cells=n,mask_mirror_mismatches=int(np.sum(s.geometry.solid!=s.geometry.solid[::-1])),
            initial_y_mirror_roundoff_m=float(np.max(abs(s.geometry.y+s.geometry.y[::-1]))),
            inlet_pair_error=reflection_error,short_reflection_relative_l2=error,
            zero_input_bitwise_match=True,positive_input=pm,negative_input=mm,numerics=s.numerics()))
        print(f"Preflight N{n}: mirror check passed; selected {best} threads; estimated {report['runtime_estimates'][str(n)]['estimated_solver_hours_per_10_cycles']:.2f} solver hours / 10 cycles",flush=True)
    write(folder/'preflight.json',report)
    return report


def scalar_row(s,d):
    row={'time_s':s.time}
    for k,v in d.items():
        if isinstance(v,(float,int,bool,np.floating,np.integer)):
            row[k]=float(v)
    return row


def safe_metrics(s,d):
    s.check_stability()
    if d['mach']>=.1:raise RuntimeError(f"Mach limit exceeded: {d['mach']}")
    if d['density_variation']>=.01:raise RuntimeError(f"Density limit exceeded: {d['density_variation']}")
    if s.iteration and abs(d['mass_balance_relative'])>=1e-8:raise RuntimeError('Discrete mass residual exceeded 1e-8')


def final_analysis(folder,spec,s):
    d=analysis.load(folder/'diagnostics/fixed_probes.csv')
    p=analysis.load(folder/'diagnostics/signed_wall_profiles.csv')
    end=spec['duration_s']
    refs=analysis.healthy_references(asdict(s.config),.03)
    out=dict(kind=spec['kind'],duration_s=s.time,simulation_completed=True,solver_validated=False)
    jets=analysis.load(folder/'jet_probes.csv')
    jet_left=max(0.,end-(.25 if spec['kind']=='healthy' else 60/s.config.heart_rate))
    out['jet_centroid_window_s']=[jet_left,end]
    out['jet_centroid_mm']={key:analysis.window_stats(jets['time_s'],jets[key],jet_left,end)
        for key in jets if key.endswith('_forward_centroid_y_mm')}
    if spec['kind']=='healthy':
        left=end-.25
        if left<0:
            out['gate']='SMOKE_ONLY';return out
        stats={key:analysis.window_stats(d['time_s'],d[key],left,end) for key in analysis.SCALARS}
        pe=100*abs(stats['dp20_50_pa']['mean']/refs['pressure_pa']-1)
        we=100*abs(stats['roi_mean_abs_wss_pa']['mean']/refs['wss_pa']-1)
        qe=100*abs(stats['x030_flow_m2_s']['mean']/refs['flow_m2_s']-1)
        drift=100*(stats['dp20_50_pa']['max']-stats['dp20_50_pa']['min'])/refs['pressure_pa']
        out.update(reference=refs,window_s=[left,end],legacy_stats=stats,pressure_error_percent=pe,
            legacy_mean_wss_error_percent=we,flow_error_percent=qe,pressure_drift_percent=drift,
            healthy_gate_pass=bool(pe<5 and we<5 and qe<1 and drift<1))
        # Quadratic estimates at the final field are compared to planar truth;
        # diagnostics only, not a promotion to production measurement.
        rows=WallSampler(s).sample(s.diagnostics())
        out['diagnostic_wss_errors']={str(radius):100*np.linalg.norm(np.array([r['shear_pa'] for r in rows if r['radius_cells']==radius])-refs['wss_pa'])/(refs['wss_pa']*math.sqrt(sum(r['radius_cells']==radius for r in rows))) for radius in (4.,5.)}
    elif end>=3*60/72:
        T=60/72
        out['repeatability']=analysis.repeatability(d,p,asdict(s.config),spec['cycles'])
        out['last_cycle_stats']={k:analysis.window_stats(d['time_s'],d[k],end-T,end) for k in analysis.SCALARS}
        out['gate']='DIAGNOSTIC_ONLY_REVIEW_REQUIRED'
    else:out['gate']='SHORT_STARTUP_DIAGNOSTIC_ONLY'
    return out


def run_case(base,spec,threads,smoke=False):
    folder=base/'cases'/spec['key']/('attempt_'+stamp())
    folder.mkdir(parents=True)
    write(folder/'case_spec.json',spec)
    write(folder/'status.json',dict(status='running',started=stamp()))
    s=make_solver(spec,threads)
    intervention=control_initialize(s,spec.get('initialization','native'))
    write(folder/'initialization_control.json',intervention)
    original,driven,perturb=initialize_pair(s,spec['epsilon'])
    T=60/s.config.heart_rate;cutoff=math.ceil(T/s.dt) if spec['epsilon'] else 0
    perturb.update(requested_end_s=T if spec['epsilon'] else 0.,actual_end_step=cutoff,actual_end_s=cutoff*s.dt)
    write(folder/'perturbation.json',perturb)
    write(folder/'input.json',asdict(s.config));write(folder/'numerics.json',s.numerics())
    write(folder/'geometry_report.json',s.geometry_report)
    write(folder/'initialization.json',dict(initial_fluid_populations_sha256=hashlib.sha256(s.f[s.nodes].tobytes()).hexdigest(),inlet_profile_sha256=hashlib.sha256(driven.tobytes()).hexdigest(),threads=threads,epsilon=spec['epsilon'],reflected_initial_state=spec['epsilon']<0))
    np.savez_compressed(folder/'inlet_profiles.npz',y_m=s.geometry.y,original_lattice=original,driven_lattice=driven)
    probes=Probes(s,folder/'diagnostics',spec['duration_s'])
    wall=WallSampler(s)
    asymmetry=CSV(folder/'symmetry_probes.csv')
    global_csv=CSV(folder/'global_timeseries.csv');jets=CSV(folder/'jet_probes.csv');wall_csv=CSV(folder/'phase_wall_estimates.csv')
    began=time.perf_counter();last_print=began;sample=0;endstep=math.ceil(spec['duration_s']/s.dt)
    field_counter=0
    try:
        while True:
            d=None
            sample_step=math.ceil(sample_time(sample)/s.dt-1e-9)
            if s.iteration>=sample_step or s.iteration==endstep:
                d=s.diagnostics();safe_metrics(s,d)
                # No previous timestep exists at t=0, even though the scratch
                # buffer was initialized. Mark storage/flux placeholders invalid.
                if s.iteration==0:
                    for k in ['control_volume_mass_per_depth_kg_m','link_inlet_mass_flux_kg_m_s','link_outlet_mass_flux_kg_m_s','mass_storage_rate_kg_m_s','mass_balance_residual_kg_m_s','mass_balance_relative']:d[k]=0.
                d['mass_balance_sample_valid']=bool(s.iteration)
                if s.iteration==endstep:probes.next_profile=s.time
                global_csv.row(scalar_row(s,d));probes.record(s,d);asymmetry.row(symmetry_row(s,d))
                jr=dict(time_s=s.time,epsilon_applied_last_step=spec['epsilon'] if 0<s.iteration<=cutoff else 0.)
                for xmm in [20,30,35,40,45,50,55]:
                    j=int(round(xmm/1000/s.dx));mask=s.geometry.fluid[:,j];uu=d['u'][mask,j]
                    weight=np.maximum(uu,0)**2
                    jr[f'x{xmm:03}_forward_centroid_y_mm']=float((s.geometry.y[mask]@weight)/max(weight.sum(),1e-30)*1000)
                jets.row(jr);sample+=1
            if probes.phase_due(s):
                if d is None:d=s.diagnostics();safe_metrics(s,d)
                probes.snapshot(s)
                for row in wall.sample(d):wall_csv.row(row)
            if s.iteration>=endstep:break
            events=[math.ceil(sample_time(sample)/s.dt-1e-9),endstep,s.iteration+2000]
            if probes.next_frame<len(probes.phase_targets):events.append(math.ceil(probes.phase_targets[probes.next_frame]/s.dt-1e-9))
            target=min(v for v in events if v>s.iteration)
            advance_to(s,target,original,driven,cutoff)
            now=time.perf_counter()
            if now-last_print>20:
                print(f"{spec['key']}: {s.time:.4f}/{spec['duration_s']:.4f} simulated s; elapsed {(now-began)/60:.1f} min",flush=True)
                write(folder/'progress.json',dict(time_s=s.time,target_s=spec['duration_s'],elapsed_s=now-began))
                last_print=now
        probes.close();global_csv.close();jets.close();wall_csv.close();asymmetry.close()
        result=final_analysis(folder,spec,s)
        result['elapsed_s']=time.perf_counter()-began
        write(folder/'analysis.json',result)
        write(folder/'status.json',dict(status='completed',finished=stamp(),elapsed_s=result['elapsed_s']))
        receipt={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file()}
        write(folder/'COMPLETED_FILES.json',receipt)
        print(f"Completed {spec['key']} in {result['elapsed_s']/60:.1f} min",flush=True)
        return folder,result
    except BaseException as exc:
        probes.close();global_csv.close();jets.close();wall_csv.close();asymmetry.close()
        write(folder/'status.json',dict(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed',error=str(exc),time_s=s.time))
        (folder/'error.log').write_text(traceback.format_exc(),encoding='utf-8')
        raise


def completed(base,key):
    cases=sorted((base/'cases'/key).glob('attempt_*'),reverse=True)
    for folder in cases:
        receipt=folder/'COMPLETED_FILES.json'
        if receipt.exists() and read(folder/'status.json')['status']=='completed':
            for name,h in read(receipt).items():
                if sha(folder/name)!=h:raise ValueError(f'Completed result changed: {folder/name}')
            return folder,read(folder/'analysis.json')
    return None

