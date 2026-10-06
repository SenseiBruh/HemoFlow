"""Longer cycle-settling and domain-sensitivity experiments with a preserved solver."""
import argparse
from pathlib import Path
import shutil
import sys
import os
import platform
import traceback
import math
import numpy as np
import numba
from experiment import ROOT,VERSION,make_solver,initialize_pair,advance_to,reflected_populations
from case_runtime import read,write,sha,stamp,run_case,completed
from startup_controls import MODES,initialize
import sequence_support as U
import startup_analysis as A
import settling_analysis as C
import grid_analysis as G
from safe_export import export_review


def plan(cycles=10,smoke=False):
    specs=[]
    for kind,n in [('healthy',80),('healthy',120),('severe',80),('severe',120),('healthy',160),('severe',160)]:
        specs.append(dict(key=f'{kind}_N{n:03}',kind=kind,cells=n,length_mm=150,
            dt_s=(1.5e-6 if kind=='healthy' else .3e-6)*80/n,
            epsilon=0. if kind=='healthy' else .001,cycles=None if kind=='healthy' else cycles,
            duration_s=.002 if smoke else (2. if kind=='healthy' else cycles*60/72),initialization='native'))
    return specs


def preflight(base,threads,specs):
    rows=[];tests=U.manufactured_checks()
    if not all(r['passed'] for r in tests):raise ValueError('Manufactured measurement check failed')
    for spec in specs:
        s=make_solver(spec,threads);report=initialize(s,spec['initialization'])
        o,d,_=initialize_pair(s,spec['epsilon']);advance_to(s,64,o,d,math.ceil((60/72)/s.dt) if spec['epsilon'] else 0)
        s.check_stability()
        rows.append(dict(key=spec['key'],initialization=report,numerics=s.numerics(),steps=64))
    for spec in [v for v in specs if v['kind']=='severe']:
        p=make_solver(spec,threads);m=make_solver(spec,threads)
        po,pd,_=initialize_pair(p,.001);mo,md,_=initialize_pair(m,-.001);cut=math.ceil((60/72)/p.dt)
        advance_to(p,64,po,pd,cut);advance_to(m,64,mo,md,cut)
        mirrored=reflected_populations(p.f,p.ny,p.nx)
        error=float(np.linalg.norm((m.f-mirrored)[m.nodes])/np.linalg.norm(mirrored[m.nodes]))
        if error>1e-9:raise RuntimeError('Short reflection check failed')
        rows.append(dict(length_mm=spec['length_mm'],short_reflection_relative_l2=error))
    result=dict(checks=rows,manufactured=tests,environment=dict(python=sys.version,numpy=np.__version__,numba=numba.__version__,platform=platform.platform(),threads=threads))
    write(base/'preflight.json',result);print('All initialization, numerical and short reflection preflights passed.',flush=True)


def pairs():
    return [('severe_N080','severe_N120'),('severe_N120','severe_N160')]


def analyze(folder,key,base,smoke):
    out=base/'analysis'/key;out.mkdir(parents=True,exist_ok=True)
    result=A.summarize(folder,out,smoke)
    if result['kind']=='severe' and not smoke:
        result['cycle_settling']=C.summarize(folder,out)
        write(out/'analysis.json',result)
    for name in ('case_spec.json','input.json','numerics.json','initialization_control.json','initialization.json','perturbation.json','COMPLETED_FILES.json','geometry_report.json'):
        shutil.copyfile(folder/name,out/name)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--threads',type=int,default=None);ap.add_argument('--cycles',type=int,default=10)
    ap.add_argument('--smoke',action='store_true');ap.add_argument('--plan-only',action='store_true')
    ap.add_argument('--resume',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--export-only',type=Path)
    args=ap.parse_args()
    if args.export_only:return 0 if export_review(args.export_only.resolve()) else 1
    if not 4<=args.cycles<=30:ap.error('cycles must be 4 through 30')
    if args.resume and args.output:ap.error('Choose resume or output')
    threads=args.threads if args.threads is not None else min(16,numba.config.NUMBA_NUM_THREADS)
    if not 1<=threads<=numba.config.NUMBA_NUM_THREADS:ap.error('Unavailable thread count')
    if args.plan_only:print(__import__('json').dumps(dict(threads=threads,cases=plan(args.cycles,args.smoke)),indent=2));return 0
    hashes=U.fingerprints()
    if args.resume:
        base=args.resume.resolve();info=read(base/'study.json')
        if info['source_hashes']!=hashes:raise ValueError('Package changed. Resume with matching source.')
        if args.threads is not None and args.threads!=info['threads']:raise ValueError('Resume cannot change threads')
        threads=info['threads'];specs=info['cases'];smoke=info['smoke']
        if threads>numba.config.NUMBA_NUM_THREADS:raise ValueError('Saved thread count unavailable')
    else:
        base=args.output.resolve() if args.output else ROOT/'experiments'/('grid_refinement_'+stamp())
        base.mkdir(parents=True,exist_ok=False);smoke=args.smoke;specs=plan(args.cycles,smoke)
        info=dict(version=VERSION,created=stamp(),cases=specs,threads=threads,smoke=smoke,source_hashes=hashes,
            purpose='Grid sensitivity at L150 with fixed dx/dt; unchanged solver; separate timestep study still required.',
            sampling='50 microseconds through 20 ms; 250 microseconds thereafter; actual solver times recorded.',
            onset_thresholds_m_s=A.ONSET_LEVELS,comparison_windows_s=A.WINDOWS)
        write(base/'study.json',info)
        for name in hashes:
            dst=base/'source_snapshot'/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,dst)
    lock=base/'RUNNING.lock'
    try:
        with lock.open('x') as f:f.write(f'PID {os.getpid()}')
    except FileExistsError:raise RuntimeError('Study lock exists; confirm old process has stopped before removing a stale lock.')
    status=dict(status='running',cases=[],comparisons=[],solver_validated=False,smoke_only=smoke,threads=threads)
    records=read(base/'case_records.json') if (base/'case_records.json').exists() else {}
    folders={};code=0
    print(f'Results: {base}\nFixed threads: {threads}. Ten cycles per severe case by default. Keep the PC awake.',flush=True)
    try:
        preflight(base,threads,specs)
        trials={};remaining=0.
        for spec in specs:
            if completed(base,spec['key']):continue
            key=(spec['kind'],spec['cells'],spec['dt_s'])
            if key not in trials:trials[key]=U.benchmark(spec,threads)
            remaining+=trials[key]['estimated_solver_hours']
        write(base/'runtime_estimates.json',dict(estimated_remaining_solver_hours=remaining,trials={str(k):v for k,v in trials.items()},notes='Excludes recording, analysis, compression and thermal throttling.'))
        print(f'Estimated {remaining:.2f} solver hours plus recording and analysis.',flush=True)
        for spec in specs:
            if spec['cells']==160 and not smoke and not G.allow_fine(status['comparisons']):
                status['status']='grid_review_required_before_N160'
                print('80/120 comparison did not pass all gates; 160-cell stage deferred.',flush=True)
                break
            key=spec['key'];done=completed(base,key)
            if done:
                folder=done[0]
                if read(folder/'case_spec.json')!=spec or read(folder/'numerics.json')['compute_threads']!=threads:raise ValueError('Completed case inputs differ')
                print('Reusing completed case: '+key,flush=True)
            else:folder,_=run_case(base,spec,threads,smoke)
            records[key]=dict(folder=str(folder.resolve()));write(base/'case_records.json',records)
            result=analyze(folder,key,base,smoke);folders[key]=folder
            status['cases'].append(dict(key=key,folder=str(folder),analysis=result));write(base/'batch_status.json',status)
            if not result['basic_pass']:raise RuntimeError('Numerical health check failed: '+key)
            if spec['kind']=='healthy' and not smoke and not result['healthy_gate_pass']:
                status['status']='healthy_review_required';break
            for ka,kb in pairs():
                if ka not in folders or kb not in folders:continue
                label=ka+'__'+kb
                if any(r['pair']==label for r in status['comparisons']):continue
                comparison=G.compare(folders[ka],folders[kb],base/'analysis',label,smoke)
                status['comparisons'].append(comparison)
            write(base/'batch_status.json',status)
        else:status['status']='smoke_completed_not_scientific' if smoke else 'grid_completed_review_required'
    except BaseException as exc:
        code=130 if isinstance(exc,KeyboardInterrupt) else 1
        status.update(status='interrupted' if code==130 else 'failed',error=str(exc))
        (base/'error.log').write_text(traceback.format_exc(),encoding='utf-8');print('Stopped: '+str(exc),flush=True)
    finally:
        try:
            write(base/'batch_status.json',status)
            (base/'RESULTS.md').write_text('# Grid-refinement batch\n\nStatus: '+status['status']+'\n\nReview cycle_settling.json and pair-specific grid_comparison.json. Passing is limited to recorded observables; grid/timestep independence and physical validation remain unestablished. Upload the newest successfully printed ZIP in exports. Raw data and full fields remain in cases.\n',encoding='utf-8')
            archive=export_review(base)
            print(f'Status: {status["status"]}\nResume:\n& "{sys.executable}" "{ROOT / "run_grid.py"}" --resume "{base}"',flush=True)
        finally:lock.unlink(missing_ok=True)
    return code

if __name__=='__main__':sys.exit(main())
