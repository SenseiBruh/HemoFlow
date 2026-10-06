"""Raw jet histories and explicitly diagnostic reflection comparisons."""
from pathlib import Path
import numpy as np
import sequence_analysis as S
from case_runtime import read,write

THRESHOLD_MM=.05
REGIONS={'all':(20,55),'stenosis':(22.5,37.5),'downstream':(37.5,55),'far_downstream':(50,55)}


def branch_events(t,y,threshold=THRESHOLD_MM):
    """Hysteresis: retain last side while inside the central band."""
    last=0;events=[]
    for time,value in zip(t,y):
        side=1 if value>threshold else -1 if value<-threshold else 0
        if side and side!=last:
            events.append(dict(time_s=float(time),from_side=last,to_side=side,centroid_mm=float(value)))
            last=side
    return events


def jet_history(folder,out):
    d=S.A.load(Path(folder)/'jet_probes.csv');t=d['time_s']
    if not np.all(np.diff(t)>0):raise ValueError('Jet timestamps are not increasing')
    spec=read(Path(folder)/'case_spec.json');T=60/72;rows=[];events={}
    for name,y in d.items():
        if not name.endswith('_forward_centroid_y_mm'):continue
        events[name]=branch_events(t,y)
        for cycle in range(int(round(spec['duration_s']/T))):
            left,right=cycle*T,(cycle+1)*T
            if right>t[-1]+1e-10:continue
            st=S.A.window_stats(t,y,left,right)
            rows.append(dict(station=name,cycle=cycle+1,left_s=left,right_s=right,**st))
    if rows:S.write_csv(Path(out)/'jet_cycles.csv',rows)
    result=dict(threshold_mm=THRESHOLD_MM,definition='Forward axial-velocity-squared weighted centroid; threshold crossings are diagnostics, not a physical instability criterion.',
        samples=len(t),last_time_s=float(t[-1]),events=events,
        sampling_caveat='Events detected at recorded samples, not exact crossing times. Previously completed runs can have fewer stations and startup fields.')
    write(Path(out)/'jet_history.json',result)
    return result


def load_field(folder,entry):
    with np.load(Path(folder)/'diagnostics'/entry['file']) as z:
        return {k:z[k].copy() for k in ('x_m','y_m','solid','u_m_s','v_m_s','pressure_pa','time_s','lower_wss_pa','upper_wss_pa')}


def vector_error(a,b,mask,reflect=False):
    bu=b['u_m_s'][::-1] if reflect else b['u_m_s']
    bv=-b['v_m_s'][::-1] if reflect else b['v_m_s']
    den=np.sum(a['u_m_s'][mask]**2+a['v_m_s'][mask]**2)
    if den<=1e-30:return None
    return float(100*np.sqrt(np.sum((a['u_m_s'][mask]-bu[mask])**2+(a['v_m_s'][mask]-bv[mask])**2)/den))


def reflection_pair(a,b,out,label):
    """Compare exact common nodes. Pressure is gauge centered for this diagnostic."""
    a,b=Path(a),Path(b)
    ia,ib=[read(p/'diagnostics/phase_fields.json') for p in (a,b)]
    rows=[];spec=read(a/'case_spec.json');T=60/72
    for ea in ia:
        matches=[eb for eb in ib if abs(eb['requested_time_s']-ea['requested_time_s'])<1e-9]
        if not matches:continue
        eb=matches[0];fa,fb=load_field(a,ea),load_field(b,eb)
        nx=min(len(fa['x_m']),len(fb['x_m']))
        if not np.allclose(fa['x_m'][:nx],fb['x_m'][:nx],rtol=0,atol=1e-12) or not np.allclose(fa['y_m'],-fb['y_m'][::-1],rtol=0,atol=1e-12):raise ValueError('Reflection diagnostic requires matching physical nodes')
        for f in (fa,fb):
            for key in ('solid','u_m_s','v_m_s','pressure_pa'):f[key]=f[key][:,:nx]
        if not np.array_equal(fa['solid'],fb['solid']) or not np.array_equal(fa['solid'],fb['solid'][::-1]):raise ValueError('Common geometry is not identical and symmetric')
        X=fa['x_m'][:nx]*1000
        for region,(lo,hi) in REGIONS.items():
            mask=(~fa['solid'])&((X>=lo-1e-9)&(X<=hi+1e-9))[None,:]
            if not mask.any():raise ValueError('Empty diagnostic region')
            for reflected in (False,True):
                bp=fb['pressure_pa'][::-1] if reflected else fb['pressure_pa']
                p=fa['pressure_pa'][mask];q=bp[mask];p=p-p.mean();q=q-q.mean()
                w=(X>=lo-1e-9)&(X<=hi+1e-9)
                wa=np.r_[fa['lower_wss_pa'][:nx][w],fa['upper_wss_pa'][:nx][w]]
                wb=np.r_[fb['upper_wss_pa' if reflected else 'lower_wss_pa'][:nx][w],fb['lower_wss_pa' if reflected else 'upper_wss_pa'][:nx][w]]
                rows.append(dict(pair=label,requested_time_s=ea['requested_time_s'],actual_time_a_s=float(fa['time_s']),actual_time_b_s=float(fb['time_s']),
                    last_cycle=ea['requested_time_s']>=spec['duration_s']-T-1e-9,region=region,orientation='reflected' if reflected else 'direct',
                    velocity_l2_percent=vector_error(fa,fb,mask,reflected),
                    gauge_centered_pressure_l2_percent=float(100*np.linalg.norm(p-q)/np.linalg.norm(p)) if np.linalg.norm(p)>1e-12 else None,
                    signed_wss_l2_percent=float(100*np.linalg.norm(wa-wb)/np.linalg.norm(wa)) if np.linalg.norm(wa)>1e-12 else None))
    if not rows:raise ValueError('No common phase fields')
    S.write_csv(Path(out)/(label+'_reflection.csv'),rows)
    summary=[]
    for region in REGIONS:
        for orientation in ('direct','reflected'):
            r=[v for v in rows if v['last_cycle'] and v['region']==region and v['orientation']==orientation]
            vals=[v['velocity_l2_percent'] for v in r if v['velocity_l2_percent'] is not None]
            summary.append(dict(region=region,orientation=orientation,last_cycle_phases=len(r),max_velocity_l2_percent=max(vals) if vals else None))
    return dict(pair=label,summary=summary,diagnostic_only=True,notes='Reflection maps u(y) to u(-y), v(y) to -v(-y), and swaps upper/lower +x wall shear. It does not replace direct acceptance scores or establish domain independence.')
