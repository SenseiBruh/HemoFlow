"""Phase-aligned cycle diagnostics, with complete-window checks and no filtering.

All percentage L2 scores use the earlier cycle (or shorter domain) as reference.
A 1% repeatability screen must hold for the last three consecutive transitions.
This is a finite-duration diagnostic, never a physical validation verdict.
"""
from pathlib import Path
import numpy as np
import strict_analysis as A
from case_runtime import read,write
from sequence_analysis import write_csv

REGIONS={'full':(20,55),'throat':(28,32),'downstream':(37.5,55),'far':(50,55)}

def streams(folder):
    folder=Path(folder)
    return (A.load(folder/'diagnostics/fixed_probes.csv'),A.load(folder/'diagnostics/signed_wall_profiles.csv'))

def cycle_rows(d,p,T,cycles):
    metrics=[];repeat=[]
    for c in range(1,cycles+1):
        l,r=(c-1)*T,c*T
        for name in A.SCALARS:
            metrics.append(dict(cycle=c,metric=name,**A.window_stats(d['time_s'],d[name],l,r)))
        if c<2:continue
        row=dict(cycle=c,reference_cycle=c-1)
        for name in A.SCALARS:
            row[name+'_l2_percent']=A.relative_l2(d['time_s'],d[name],d['time_s']+T,d[name],l,r)
        shifted=dict(p,time_s=p['time_s']+T)
        for name,region in REGIONS.items():
            row['signed_wss_'+name+'_l2_percent']=A.profiles_l2(p,shifted,l,r,region)
        repeat.append(row)
    return metrics,repeat

def summarize(folder,out):
    spec=read(folder/'case_spec.json');T=A.period(read(folder/'input.json'))
    d,p=streams(folder);metrics,rows=cycle_rows(d,p,T,spec['cycles'])
    write_csv(out/'cycle_metrics.csv',metrics);write_csv(out/'cycle_repeatability.csv',rows)
    tail=rows[-3:]
    worst=max(v for row in tail for k,v in row.items() if k.endswith('_percent'))
    result=dict(cycles=spec['cycles'],comparison_count=len(rows),last_three_transitions=tail,
        worst_last_three_l2_percent=worst,threshold_percent=1.,
        recorded_observables_repeatable=bool(len(tail)==3 and worst<1.),
        solver_validated=False,scope='Pressure, flow and signed WSS in fixed regions only; does not prove the entire field is periodic or grid independent.')
    write(out/'cycle_settling.json',result)
    return result

def compare(a,b,out,kind="domain"):
    ca,cb=read(a/'input.json'),read(b/'input.json')
    A.matches_input(ca,cb,kind,read(a/'numerics.json'),read(b/'numerics.json'))
    sa,sb=read(a/'case_spec.json'),read(b/'case_spec.json')
    if sa['cycles']!=sb['cycles'] or sa['epsilon']!=sb['epsilon']:raise ValueError('Forcing/duration mismatch')
    T=A.period(ca);da,pa=streams(a);db,pb=streams(b);rows=[]
    for c in range(1,sa['cycles']+1):
        l,r=(c-1)*T,c*T;row=dict(cycle=c)
        for name in A.SCALARS:
            va=A.window_stats(da['time_s'],da[name],l,r);vb=A.window_stats(db['time_s'],db[name],l,r)
            row[name+'_reference_mean']=va['mean'];row[name+'_comparison_mean']=vb['mean']
            row[name+'_mean_difference_percent']=100*(vb['mean']-va['mean'])/abs(va['mean'])
            row[name+'_l2_percent']=A.relative_l2(db['time_s'],db[name],da['time_s'],da[name],l,r)
        for name,region in REGIONS.items():row['signed_wss_'+name+'_l2_percent']=A.profiles_l2(pb,pa,l,r,region)
        rows.append(row)
    limits={'dp20_50_pa_mean_difference_percent':2.,'dp20_50_pa_l2_percent':5.,
            'roi_mean_abs_wss_pa_mean_difference_percent':5.,'roi_mean_abs_wss_pa_l2_percent':5.,
            'roi_peak_abs_wss_pa_l2_percent':10.,'x030_flow_m2_s_mean_difference_percent':1.,
            'x030_flow_m2_s_l2_percent':2.,**{'signed_wss_'+n+'_l2_percent':10. for n in REGIONS}}
    passed=len(rows)>=3 and all(abs(row[k])<=lim for row in rows[-3:] for k,lim in limits.items())
    result=dict(reference=f"N{sa['cells']} L{sa['length_mm']} mm",comparison=f"N{sb['cells']} L{sb['length_mm']} mm",cycles=rows,limits_percent=limits,
        last_three_cycle_domain_screen_passed=bool(passed),solver_validated=False,
        scope='Direct comparison only; interpret with within-case repeatability and sampled velocity fields. No automatic promotion to broad severity studies.')
    write_csv(out/'domain_cycles.csv',rows);write(out/'domain_cycles.json',result)
    return result
