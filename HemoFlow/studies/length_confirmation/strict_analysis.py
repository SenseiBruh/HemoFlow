"""Complete-window, metadata-driven checks. Raw data are never filtered or edited."""
import csv
import json
import math
import re
from pathlib import Path
import numpy as np

LIMITS = dict(mach=0.1, density_departure_percent=1.0, mass_residual_relative=1e-8,
              mean_flow_error_percent=1.0, repeat_l2_percent=1.0,
              healthy_pressure_error_percent=5.0, healthy_wss_error_percent=5.0,
              outlet_shear_l2_percent=10.0, outlet_u_profile_l2_percent=10.0,
              outlet_transverse_over_mean_percent=5.0,
              pressure_mean_percent=2.0, pressure_waveform_l2_percent=5.0,
              mean_wss_percent=5.0, peak_wss_percent=10.0, profile_l2_percent=10.0,
              flow_mean_percent=1.0, flow_waveform_l2_percent=2.0)
PROFILE = re.compile(r'^(lower|upper)_x(\d+(?:\.\d+)?)_pa$')
SCALARS = ('dp20_50_pa','roi_mean_abs_wss_pa','roi_peak_abs_wss_pa','x030_flow_m2_s')


class InvalidComparison(ValueError):
    pass


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save_json(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temp.replace(path)


def load(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f)
        try: header=next(reader)
        except StopIteration: raise InvalidComparison(f'Empty CSV: {path}')
        if len(set(header))!=len(header) or 'time_s' not in header:
            raise InvalidComparison(f'Duplicate headers or missing time_s: {path}')
        values=[]
        for i,row in enumerate(reader,2):
            if len(row)!=len(header):raise InvalidComparison(f'Row {i} has wrong width: {path}')
            try:values.append([float(v) for v in row])
            except ValueError:raise InvalidComparison(f'Non-numeric row {i}: {path}')
    a=np.asarray(values,dtype=float)
    if a.ndim!=2 or len(a)<2 or not np.isfinite(a).all():
        raise InvalidComparison(f'Need at least two finite samples: {path}')
    d=dict(zip(header,a.T)); validate_time(d['time_s'])
    if 'epoch' in d and len(np.unique(d['epoch']))!=1:
        raise InvalidComparison('Geometry/settings changed within the recording')
    return d


def validate_time(t):
    t=np.asarray(t,dtype=float)
    if t.ndim!=1 or len(t)<2 or not np.isfinite(t).all() or not np.all(np.diff(t)>0):
        raise InvalidComparison('Timestamps must be finite and strictly increasing')
    return t


def require_window(t,left,right):
    t=validate_time(t)
    if not np.isfinite([left,right]).all() or right<=left:
        raise InvalidComparison('Invalid requested interval')
    # Floating-point tolerance only, not a missing-sample/cadence allowance.
    tol=1e-10
    if left<t[0]-tol or right>t[-1]+tol:
        raise InvalidComparison(f'Incomplete window [{left:.12g}, {right:.12g}]; '
                                f'available [{t[0]:.12g}, {t[-1]:.12g}]')


def knots(ta,tb,left,right):
    require_window(ta,left,right); require_window(tb,left,right)
    t=np.union1d(ta,tb)
    return np.r_[left,t[(t>left)&(t<right)],right]


def square(t,y):
    return float(np.sum(np.diff(t)*(y[:-1]**2+y[:-1]*y[1:]+y[1:]**2)/3))


def relative_l2(ta,a,tb,b,left,right):
    tt=knots(ta,tb,left,right)
    aa=np.interp(tt,ta,a);bb=np.interp(tt,tb,b)
    if not np.isfinite(aa).all() or not np.isfinite(bb).all():
        raise InvalidComparison('Nonfinite comparison values')
    den=square(tt,bb)
    if den<=1e-25:raise InvalidComparison('Reference norm is zero or unresolved')
    return 100*math.sqrt(square(tt,aa-bb)/den)


def window_stats(t,y,left,right):
    tt=knots(t,t,left,right); yy=np.interp(tt,t,y)
    return dict(mean=float(np.trapezoid(yy,tt)/(right-left)),
                min=float(yy.min()),max=float(yy.max()))


def profile_columns(d):
    result={}
    for key in d:
        if key=='time_s':continue
        m=PROFILE.fullmatch(key)
        if not m:raise InvalidComparison(f'Unexpected wall-profile column: {key}')
        location=(m[1],float(m[2]))
        if location in result:raise InvalidComparison('Duplicate wall/location')
        result[location]=key
    lower={x for w,x in result if w=='lower'};upper={x for w,x in result if w=='upper'}
    if not lower or lower!=upper:raise InvalidComparison('Both walls must have identical physical stations')
    return result


def profiles_l2(a,b,left,right,region=None,swap_reference_walls=False):
    ca,cb=profile_columns(a),profile_columns(b)
    if ca.keys()!=cb.keys():raise InvalidComparison('Wall-profile physical locations do not match')
    tt=knots(a['time_s'],b['time_s'],left,right)
    num=den=0.; count=0
    for (wall,x),col in ca.items():
        if region is not None and not region[0]-1e-8<=x<=region[1]+1e-8:continue
        other=('upper' if wall=='lower' else 'lower') if swap_reference_walls else wall
        aa=np.interp(tt,a['time_s'],a[col])
        bb=np.interp(tt,b['time_s'],b[cb[(other,x)]])
        num+=square(tt,aa-bb);den+=square(tt,bb);count+=1
    if not count or den<=1e-25:raise InvalidComparison('No resolved reference profile in region')
    return 100*math.sqrt(num/den)


def period(config):
    bpm=float(config['heart_rate'])
    if not np.isfinite(bpm) or bpm<=0:raise InvalidComparison('Invalid heart rate')
    return 60/bpm


def healthy_references(config,separation_m):
    h=float(config['diameter_mm'])/1000
    u=float(config['mean_velocity']);mu=float(config['viscosity'])
    return dict(pressure_pa=12*mu*u*separation_m/h**2,wss_pa=6*mu*u/h,
                flow_m2_s=u*h)


def matches_input(a,b,kind,na=None,nb=None):
    # Numerical settings are intentionally allowed only for their named study.
    keys=['density','viscosity','diameter_mm','mean_velocity','heart_rate',
          'pulsatility_percent','pulse_shape','geometry_mode','plaques',
          'stenosis_center_mm','stenosis_length_mm','stenosis_percent']
    if kind!='domain':keys.append('length_mm')
    if kind!='grid':keys.append('cells_across')
    if kind!='outlet':keys.append('boundary_model')
    for key in keys:
        if key not in a or key not in b or a[key]!=b[key]:
            raise InvalidComparison(f'Unexpected input difference or missing input: {key}')
    if na is not None and nb is not None:
        if kind=='grid':
            va=na['dx_m']/na['dt_s'];vb=nb['dx_m']/nb['dt_s']
            if not math.isclose(va,vb,rel_tol=1e-10):
                raise InvalidComparison('Grid comparison changes numerical wave speed')
        elif kind in ('outlet','domain') and not math.isclose(na['dt_s'],nb['dt_s'],rel_tol=1e-10):
            raise InvalidComparison('Boundary/domain comparison changes timestep')


def repeatability(d,p,config,cycles):
    T=period(config);rows=[]
    # Do not use the startup cycle as the reference. At 3 cycles: compare 3 vs 2.
    for i in range(max(2,cycles-3),cycles):
        left,right=i*T,(i+1)*T
        row={'cycle':i+1,'reference_cycle':i}
        for k in ('dp20_50_pa','roi_mean_abs_wss_pa','x030_flow_m2_s'):
            row[k]=relative_l2(d['time_s'],d[k],d['time_s']+T,d[k],left,right)
        shifted=dict(p,time_s=p['time_s']+T)
        row['signed_profiles']=profiles_l2(p,shifted,left,right)
        rows.append(row)
    if not rows:raise InvalidComparison('Need >=3 cycles for the chosen repeatability check')
    worst=max(v for r in rows for k,v in r.items() if k not in ('cycle','reference_cycle'))
    return dict(comparisons=rows,worst_l2_percent=worst,passed=worst<LIMITS['repeat_l2_percent'])


def assess_case(attempt,spec,smoke=False):
    attempt=Path(attempt);c=read_json(attempt/'input.json');T=period(c)
    d=load(attempt/'diagnostics/fixed_probes.csv')
    p=load(attempt/'diagnostics/signed_wall_profiles.csv')
    o=load(attempt/'diagnostics/outlet_probes.csv')
    sample=read_json(attempt/'diagnostics/sample_status.json')
    final=read_json(attempt/'final/metrics.json')
    if not sample['closed'] or sample['rows']!=len(d['time_s']) or sample['profile_rows']!=len(p['time_s']):
        raise InvalidComparison('Incomplete recorder closure or row count mismatch')
    if not np.array_equal(d['time_s'],o['time_s']):raise InvalidComparison('Outlet/fixed timestamps differ')
    raw=list((attempt/'raw').glob('*timeseries.csv'))
    if len(raw)!=1:raise InvalidComparison('Need one original scalar recording per case')
    with raw[0].open(encoding='utf-8-sig') as f:
        global_times=np.array([float(r['time_s']) for r in csv.DictReader(f)])
    if not np.array_equal(d['time_s'],global_times):raise InvalidComparison('Global/fixed timestamps differ')
    target=float(spec['duration_s'])
    if final.get('stopped_by_user') or final['time_s']<target-1e-10:
        raise InvalidComparison('Incomplete simulation')
    definitions=read_json(attempt/'diagnostics/definitions.json')
    actual=definitions['numerics']
    if actual['boundary_model']!=c['boundary_model'] or not math.isclose(actual['dt_s'],spec['requested_dt_s'],rel_tol=1e-10):
        raise InvalidComparison('Actual solver boundary/timestep differs from requested case')
    fields=read_json(attempt/'diagnostics/phase_fields.json')
    if sample['phase_frames']!=len(fields) or len(fields)!=len(definitions['phase_targets_s']):
        raise InvalidComparison('Phase-field inventory incomplete')
    for f in fields:
        if not (attempt/'diagnostics'/f['file']).is_file() or not -1e-10<=f['error_s']<=final['dt_s']+1e-10:
            raise InvalidComparison('Phase snapshot missing or mistimed')
    result={'case':spec['key'],'complete':True,'scope':'Regression screen, not physical validation',
            'smoke_only':smoke,'gate_pass':False,'numerics':definitions['numerics'],
            'rows':len(d['time_s']),'profile_rows':len(p['time_s'])}
    if smoke:
        result['note']='Short setup test only; no numerical acceptance verdict'
        return result
    cycles=int(spec['cycles']);left=(cycles-2)*T;right=cycles*T
    for stream in (d,p,o):require_window(stream['time_s'],left,right)
    result['window_s']=[left,right]
    result['fixed']={k:window_stats(d['time_s'],d[k],left,right) for k in SCALARS}
    numerical=dict(mach=float(o['mach'].max()),
                   density_departure_percent=float(o['density_departure_percent'].max()),
                   mass_residual_relative=float(np.abs(o['mass_residual_relative'][1:]).max()))
    result['numerical_maxima']=numerical
    result['basic_pass']=all(numerical[k]<LIMITS[k] for k in numerical)
    ref=healthy_references(c,.030)
    flow_error=100*abs(result['fixed']['x030_flow_m2_s']['mean']/ref['flow_m2_s']-1)
    result['mean_flow_error_percent']=flow_error
    result['repeatability']=repeatability(d,p,c,cycles)
    if spec['stage']=='healthy':
        shear={}
        for wall in ('lower','upper'):
            for where in ('outlet','inside1cell'):
                key=f'{where}_{wall}_wss_pa'
                shear[key]=relative_l2(o['time_s'],o[key],o['time_s'],o[f'reference5mm_{wall}_wss_pa'],left,right)
        result['outlet_shear_l2_percent']=shear
        result['outlet_u_profile_l2_max_percent']=window_stats(o['time_s'],o['outlet_u_profile_l2_percent'],left,right)['max']
        result['outlet_transverse_over_mean_max_percent']=window_stats(o['time_s'],o['outlet_transverse_over_mean_percent'],left,right)['max']
        outlet_pass=(max(shear.values())<=LIMITS['outlet_shear_l2_percent'] and
            result['outlet_u_profile_l2_max_percent']<=LIMITS['outlet_u_profile_l2_percent'] and
            result['outlet_transverse_over_mean_max_percent']<=LIMITS['outlet_transverse_over_mean_percent'])
        result['outlet_regression_pass']=outlet_pass
        analytical=True
        if c['pulsatility_percent']==0:
            pe=100*abs(result['fixed']['dp20_50_pa']['mean']/ref['pressure_pa']-1)
            we=100*abs(result['fixed']['roi_mean_abs_wss_pa']['mean']/ref['wss_pa']-1)
            analytical=pe<LIMITS['healthy_pressure_error_percent'] and we<LIMITS['healthy_wss_error_percent']
            result['analytical']={'reference':ref,'pressure_error_percent':pe,'wss_error_percent':we,'passed':analytical}
        result['gate_pass']=bool(result['basic_pass'] and flow_error<LIMITS['mean_flow_error_percent'] and
                result['repeatability']['passed'] and outlet_pass and analytical)
    else:
        # A short stenosis pilot is evidence for review, never full verification.
        result['gate_pass']=bool(result['basic_pass'] and result['repeatability']['passed'])
    failed=[]
    for key,val in numerical.items():
        if val>=LIMITS[key]:failed.append(f'{key}: {val:.6g} >= {LIMITS[key]:g}')
    if not result['repeatability']['passed']:
        failed.append(f"cycle repeatability: {result['repeatability']['worst_l2_percent']:.6g}% >= {LIMITS['repeat_l2_percent']:g}%")
    if spec['stage']=='healthy':
        if flow_error>=LIMITS['mean_flow_error_percent']:failed.append(f'mean flow error: {flow_error:.6g}%')
        if not analytical:failed.append('steady planar analytical pressure/WSS check')
        if not outlet_pass:
            failed.append(f"outlet check: shear L2 {max(shear.values()):.6g}%, velocity L2 {result['outlet_u_profile_l2_max_percent']:.6g}%, transverse {result['outlet_transverse_over_mean_max_percent']:.6g}%")
    result['failed_checks']=failed
    return result


def compare_cases(a,b,kind,cycles):
    a,b=Path(a),Path(b);ca=read_json(a/'input.json');cb=read_json(b/'input.json')
    na=read_json(a/'diagnostics/definitions.json')['numerics'];nb=read_json(b/'diagnostics/definitions.json')['numerics']
    matches_input(ca,cb,kind,na,nb)
    T=period(ca);left=(cycles-2)*T;right=cycles*T
    da=load(a/'diagnostics/fixed_probes.csv');db=load(b/'diagnostics/fixed_probes.csv')
    pa=load(a/'diagnostics/signed_wall_profiles.csv');pb=load(b/'diagnostics/signed_wall_profiles.csv')
    rows=[]
    for key,stat,lim in [('dp20_50_pa','mean','pressure_mean_percent'),
            ('roi_mean_abs_wss_pa','mean','mean_wss_percent'),('roi_peak_abs_wss_pa','max','peak_wss_percent'),
            ('x030_flow_m2_s','mean','flow_mean_percent')]:
        x=window_stats(da['time_s'],da[key],left,right)[stat]
        y=window_stats(db['time_s'],db[key],left,right)[stat]
        if abs(y)<1e-25:raise InvalidComparison('Unresolved relative difference denominator')
        diff=100*(x-y)/abs(y)
        rows.append(dict(metric=key,statistic=stat,value=x,reference_value=y,difference_percent=diff,
                         limit_percent=LIMITS[lim],passed=abs(diff)<=LIMITS[lim]))
    for key,lim in [('dp20_50_pa','pressure_waveform_l2_percent'),('x030_flow_m2_s','flow_waveform_l2_percent')]:
        val=relative_l2(da['time_s'],da[key],db['time_s'],db[key],left,right)
        rows.append(dict(metric=key,statistic='waveform_l2_percent',value=val,limit_percent=LIMITS[lim],passed=val<=LIMITS[lim]))
    direct=profiles_l2(pa,pb,left,right)
    rows.append(dict(metric='signed_wall_profiles',statistic='space_time_l2_percent',value=direct,
                     limit_percent=LIMITS['profile_l2_percent'],passed=direct<=LIMITS['profile_l2_percent']))
    diagnostics={}
    for region in [(20,28),(28,32),(32,40),(40,55)]:
        diagnostics[f'{region[0]}_{region[1]}mm']=profiles_l2(pa,pb,left,right,region=region)
    swapped=profiles_l2(pa,pb,left,right,swap_reference_walls=True)
    return dict(case=a.name,reference=b.name,kind=kind,window_s=[left,right],metrics=rows,
                passed=all(r['passed'] for r in rows),profile_regions_l2_percent=diagnostics,
                wall_swap_diagnostic_l2_percent=swapped,
                wall_swap_note='Orientation diagnostic only; never substitutes for the direct score or changes acceptance.')


def readiness(case_results,required,comparisons):
    missing=[k for k in required if k not in case_results or not case_results[k].get('complete')]
    failed=[k for k in required if not case_results.get(k,{}).get('gate_pass',False)]
    return dict(required=required,missing_or_incomplete=missing,failed_gates=failed,
                all_required_cases_pass=not missing and not failed,
                all_pair_comparisons_pass=bool(comparisons) and all(c.get('passed',False) for c in comparisons),
                ready_for_broad_severity_study=False,
                scope='Passing these short checks does not establish grid independence or physical validation.')
