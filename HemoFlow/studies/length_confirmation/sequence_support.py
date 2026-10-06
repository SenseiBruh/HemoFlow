"""Gated HemoFlow wall, domain, grid and timestep sensitivity sequence.

Headless CPU runs; frozen solver; no original data overwritten. Resume operates
at case boundaries. An interrupted case restarts in a fresh attempt directory.
"""
import argparse
from dataclasses import asdict
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
import numba
from experiment import ROOT,VERSION,make_solver
from audit_math import manufactured_checks
from case_runtime import run_case,completed,sha,read,write,stamp
import sequence_analysis as S


def plan(cycles=10,smoke=False):
    specs=[]
    entries=[('healthy_N080','healthy',80,60,1.5e-6),('healthy_N120','healthy',120,60,1e-6),
        ('D75_N080_L060_dt0300','severe',80,60,.3e-6),('D75_N120_L060_dt0200','severe',120,60,.2e-6),
        ('D75_N080_L090_dt0300','severe',80,90,.3e-6),('D75_N120_L090_dt0200','severe',120,90,.2e-6),
        ('D75_N160_L090_dt0150','severe',160,90,.15e-6),('D75_N120_L090_dt0300','severe',120,90,.3e-6),
        ('D75_N120_L090_dt0150','severe',120,90,.15e-6)]
    for key,kind,n,length,dt in entries:
        specs.append(dict(key=key,kind=kind,cells=n,length_mm=length,dt_s=dt,epsilon=0. if kind=='healthy' else .001,
            cycles=None if kind=='healthy' else cycles,duration_s=.002 if smoke else (2. if kind=='healthy' else cycles*60/72)))
    return specs


PAIRS={3:[(2,3,'grid')],4:[(2,4,'domain')],5:[(3,5,'domain'),(4,5,'grid')],
       6:[(5,6,'grid')],7:[(5,7,'timestep')],8:[(5,8,'timestep')]}


def fingerprints():
    expected=read(ROOT/'FROZEN_SOURCE.json')
    for name,value in expected.items():
        if sha(ROOT/name)!=value:raise ValueError(f'Frozen solver changed: {name}')
    paths=[p for p in ROOT.iterdir() if p.is_file() and p.suffix in ('.py','.md','.txt','.cmd','.json')]
    paths += list((ROOT/'source_frozen').glob('*.py'))
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}


def required_files(folder):
    names=['case_spec.json','status.json','input.json','numerics.json','geometry_report.json','jet_probes.csv',
        'global_timeseries.csv','diagnostics/definitions.json','diagnostics/sample_status.json',
        'diagnostics/fixed_probes.csv','diagnostics/signed_wall_profiles.csv',
        'diagnostics/wall_geometry.npz','diagnostics/phase_fields.json']
    names += ['diagnostics/'+e['file'] for e in read(folder/'diagnostics/phase_fields.json')]
    return names


def verify_receipt(folder,all_files=False):
    folder=Path(folder);receipt=read(folder/'COMPLETED_FILES.json')
    for name in (receipt if all_files else required_files(folder)):
        if name not in receipt or sha(folder/name)!=receipt[name]:raise ValueError(f'Missing or altered completed file: {folder/name}')
    if read(folder/'status.json')['status']!='completed':raise ValueError('Case not complete')


def compatible(study,folder,spec):
    info=read(study/'study.json');expected=read(ROOT/'FROZEN_SOURCE.json')
    parent=read(ROOT/'PARENT_FROZEN_SOURCE.json')
    for name,value in expected.items():
        if info['source_hashes'].get(name)!=value:raise ValueError(f'Prior solver hash differs: {name}')
    if info['source_hashes'].get('measurements.py')!=parent['measurements.py']:
        raise ValueError('Prior measurement definitions differ')
    old=read(folder/'case_spec.json')
    for key in ('kind','cells','epsilon','cycles'):
        if old[key]!=spec[key]:raise ValueError(f'Prior case differs: {key}')
    for key in ('duration_s','dt_s'):
        if not math.isclose(old[key],spec[key],rel_tol=1e-12):raise ValueError(f'Prior case differs: {key}')
    solver=make_solver(spec,1);config=asdict(solver.config);previous=read(folder/'input.json')
    irrelevant={'seed','compute_threads','fps','display_time_scale','particle_display_limit','particles_visible'}
    for key,value in config.items():
        if key not in irrelevant and previous.get(key)!=value:raise ValueError(f'Prior physical/numerical input differs: {key}')
    actual=read(folder/'numerics.json');wanted=solver.numerics()
    for key in ('dx_m','dt_s','tau','numerical_sound_speed_m_s'):
        if key in wanted and not math.isclose(actual[key],wanted[key],rel_tol=1e-12):raise ValueError(f'Prior numerics differ: {key}')
    with np.load(folder/'diagnostics/wall_geometry.npz') as z:
        if not np.array_equal(z['solid'],solver.geometry.solid):raise ValueError('Prior geometry mask differs')
    verify_receipt(folder)
    return True


def reuse_prior(studies,spec,base):
    if spec['length_mm']!=60:return None
    key=f"healthy_steady_N{spec['cells']:03}" if spec['kind']=='healthy' else f"D75_N{spec['cells']:03}_plus001"
    for study in studies:
        for folder in sorted((study/'cases'/key).glob('attempt_*'),reverse=True):
            try:
                compatible(study,folder,spec)
                return folder
            except (ValueError,KeyError,FileNotFoundError) as e:
                with (base/'reuse_checks.log').open('a',encoding='utf-8') as f:f.write(f'{folder}: not reused: {e}\n')
    return None


def benchmark(spec,requested):
    s=make_solver(spec,1)
    counts=[requested] if requested else [n for n in (1,2,4,8,16) if n<=numba.config.NUMBA_NUM_THREADS]
    timings={}
    for n in counts:
        s.set_threads(n);s.step(8)
        trials=[]
        for _ in range(3):
            t=time.perf_counter();s.step(64);trials.append((time.perf_counter()-t)/64)
        timings[n]=float(np.median(trials))
    best=min(timings,key=timings.get);s.set_threads(best)
    estimate=timings[best]*spec['duration_s']/spec['dt_s']/3600
    result=dict(threads=best,seconds_per_step=timings[best],estimated_solver_hours=estimate,
        thread_trials=timings,numerics=s.numerics(),caveat='Short estimate excludes diagnostics, compression and thermal throttling.')
    print(f"{spec['key']}: {best} threads, dt {s.dt*1e6:.3f} microseconds, tau {s.tau:.9f}, estimated {estimate:.2f} solver hours",flush=True)
    return result


def pack_review(base):
    target=base/'UPLOAD_REVIEW.zip';tmp=base/'UPLOAD_REVIEW.tmp'
    with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(base.rglob('*')):
            if not p.is_file() or 'uploads' in p.relative_to(base).parts or p.name=='RUNNING.lock':continue
            if p.suffix in ('.json','.py','.md','.txt','.log') or p.name in ('fixed_radius_wall_fits.csv','phase_cross_sections.csv','pair_summary.csv'):
                z.write(p,p.relative_to(base))
    tmp.replace(target)
    return target


def pack_raw(base):
    """Normal standalone ZIPs; large individual files split into lossless chunks."""
    records=read(base/'case_records.json');out=base/'uploads'/stamp();out.mkdir(parents=True)
    index=[];part=0;limit=32*1024*1024
    for key,record in records.items():
        folder=Path(record['folder'])
        verify_receipt(folder)
        for name in required_files(folder):
            p=folder/name;original_hash=sha(p);size=p.stat().st_size;chunks=max(1,math.ceil(size/limit))
            with p.open('rb') as f:
                for i in range(chunks):
                    part+=1;zipname=f'raw_{part:04}.zip';entry=f'{key}/{name}'+(f'.bytechunk_{i+1:03}' if chunks>1 else '')
                    with zipfile.ZipFile(out/zipname,'w',zipfile.ZIP_DEFLATED) as z:z.writestr(entry,f.read(limit))
                    index.append(dict(archive=zipname,entry=entry,original=f'{key}/{name}',chunk=i+1,chunks=chunks,original_sha256=original_hash,archive_sha256=sha(out/zipname)))
    write(out/'UPLOAD_INDEX.json',dict(instructions='Extract normal ZIPs; concatenate bytechunk entries in numbered order to reconstruct each original; verify SHA256.',files=index))
    print(f'Raw upload parts: {out}')


def save_comparisons(base,rows):
    write(base/'comparisons.json',rows)
    flat=[]
    for pair in rows:
        for r in pair['metrics']:
            flat.append(dict(reference=pair['reference'],comparison=pair['comparison'],kind=pair['kind'],metric=r['metric'],statistic=r['statistic'],value_percent=r['value_percent'],limit_percent=r['limit_percent'],passed=r['passed']))
    if flat:S.write_csv(base/'pair_summary.csv',flat)

