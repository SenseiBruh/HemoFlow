"""Startup sensitivity diagnostics. No convergence verdict from short runs."""
from pathlib import Path
import numpy as np
import sequence_analysis as S
import jet_analysis as J
from case_runtime import read,write

WINDOWS=((0.,.02),(.02,.1),(.1,.3),(.3,.8),(.8,1.6))
ONSET_LEVELS=(1e-9,1e-6,1e-3)

def summarize(folder,out,smoke=False):
    spec,d,p,g=S.inventory(folder)
    maxima=dict(mach=float(g['mach'].max()),density_departure_percent=float(g['density_variation'].max()*100),mass_residual_relative=float(abs(g['mass_balance_relative'][1:]).max()))
    basic=all(maxima[k]<S.A.LIMITS[k] for k in maxima)
    fit=S.fixed_wall_fits(folder,out)
    jets=J.jet_history(folder,out)
    symmetry=S.A.load(folder/'symmetry_probes.csv');onsets=[]
    if not np.array_equal(symmetry['time_s'],g['time_s']):raise ValueError('Symmetry timestamps differ')
    for name,values in symmetry.items():
        if name=='time_s':continue
        for limit in ONSET_LEVELS:
            idx=np.flatnonzero(values>=limit)
            onsets.append(dict(station=name,threshold_m_s=limit,first_sample_s=float(symmetry['time_s'][idx[0]]) if len(idx) else None,
                previous_sample_s=float(symmetry['time_s'][idx[0]-1]) if len(idx) and idx[0]>0 else None,max_asymmetry_m_s=float(values.max())))
    S.write_csv(out/'startup_onsets.csv',onsets)
    result=dict(kind=spec['kind'],initialization=spec['initialization'],complete=True,basic_pass=basic,maxima=maxima,
        smoke_only=smoke,solver_validated=False,convergence_assessed=False,jet_history=jets,fixed_radius_fits=fit,
        notes='Thresholds are diagnostic resolutions, not instability criteria. No filtering. Positive inlet skew intentionally breaks symmetry.')
    if spec['kind']=='healthy' and not smoke:
        gate=S.summarize(folder,False);errors=fit.get('final_healthy_rms_error_percent')
        result['healthy_gate_pass']=bool(gate['gate_pass'] and errors and max(errors.values())<5 and fit['invalid_rows']==0)
        result['healthy_benchmark']=gate
    write(out/'analysis.json',result)
    return result

def compare(a,b,out,label):
    sa,sb=read(a/'case_spec.json'),read(b/'case_spec.json')
    if sa['epsilon']!=sb['epsilon'] or sa['dt_s']!=sb['dt_s'] or sa['duration_s']!=sb['duration_s']:raise ValueError('Uncontrolled forcing, timestep or duration')
    ca,cb=read(a/'input.json'),read(b/'input.json')
    S.A.matches_input(ca,cb,'domain' if sa['length_mm']!=sb['length_mm'] else 'repeat',read(a/'numerics.json'),read(b/'numerics.json'))
    if ca['compute_threads']!=cb['compute_threads']:raise ValueError('Thread mismatch')
    da,db=[S.A.load(p/'diagnostics/fixed_probes.csv') for p in (a,b)]
    rows=[]
    for left,right in WINDOWS:
        if right>sa['duration_s']+1e-10:continue
        for metric in ['dp20_50_pa','roi_mean_abs_wss_pa','x030_flow_m2_s','x020_pressure_pa','x040_pressure_pa','x055_pressure_pa']+[f'x{x:03}_pressure_pa' for x in range(60,min(sa['length_mm'],sb['length_mm']),5)]:
            ta,tb=da['time_s'],db['time_s'];x=da[metric];y=db[metric]
            va=S.A.window_stats(ta,x,left,right);vb=S.A.window_stats(tb,y,left,right)
            rows.append(dict(pair=label,metric=metric,left_s=left,right_s=right,reference_mean=va['mean'],comparison_mean=vb['mean'],mean_difference=vb['mean']-va['mean'],
                waveform_l2_percent=S.A.relative_l2(tb,y,ta,x,left,right),diagnostic_only=True))
    if rows:S.write_csv(out/(label+'_startup_windows.csv'),rows)
    old_regions=J.REGIONS
    short=min(sa['length_mm'],sb['length_mm'])
    J.REGIONS={**old_regions,'common_tail':(55,short-1),'short_outlet_neighborhood':(short-10,short-1)}
    try:
        reflection=J.reflection_pair(a,b,out,label)
        reflection['region_bounds_mm']=J.REGIONS.copy()
    finally:J.REGIONS=old_regions
    return dict(pair=label,reference_mode=sa['initialization'],comparison_mode=sb['initialization'],reference_length_mm=sa['length_mm'],comparison_length_mm=sb['length_mm'],
        windows=rows,field_diagnostics=reflection,diagnostic_only=True,notes='Direct and reflected fields retained. Reflection and initialization changes do not correct or validate the solver.')
