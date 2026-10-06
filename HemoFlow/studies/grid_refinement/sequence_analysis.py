"""Predeclared comparisons of raw recordings; no filtering or field mutation."""
import csv
import math
from pathlib import Path
import numpy as np
import strict_analysis as A
from audit_math import build_stencil

RADII_MM=(.20,.25)
FIELD_LIMIT_PERCENT=10.0
REGIONS={'all':(20,55),'upstream':(20,22),'stenosis':(22.5,37.5),'downstream':(37.5,55)}


def write_csv(path,rows):
    if not rows:raise ValueError('No rows to write')
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def inventory(folder):
    folder=Path(folder);spec=A.read_json(folder/'case_spec.json')
    d=A.load(folder/'diagnostics/fixed_probes.csv');p=A.load(folder/'diagnostics/signed_wall_profiles.csv')
    g=A.load(folder/'global_timeseries.csv')
    sample=A.read_json(folder/'diagnostics/sample_status.json')
    defs=A.read_json(folder/'diagnostics/definitions.json')
    fields=A.read_json(folder/'diagnostics/phase_fields.json')
    if not sample['closed'] or sample['rows']!=len(d['time_s']) or sample['profile_rows']!=len(p['time_s']):
        raise A.InvalidComparison('Recorder incomplete or row count mismatch')
    if not np.array_equal(g['time_s'],d['time_s']):raise A.InvalidComparison('Global and probe timestamps differ')
    for stream in (d,p):A.require_window(stream['time_s'],0,spec['duration_s'])
    if len(fields)!=len(defs['phase_targets_s']) or sample['phase_frames']!=len(fields):
        raise A.InvalidComparison('Incomplete phase inventory')
    for row,target in zip(fields,defs['phase_targets_s']):
        if abs(row['requested_time_s']-target)>1e-10 or not -1e-10<=row['error_s']<=spec['dt_s']+1e-10:
            raise A.InvalidComparison('Mistimed phase field')
        if not (folder/'diagnostics'/row['file']).exists():raise A.InvalidComparison('Missing phase field')
    return spec,d,p,g


def summarize(folder,smoke=False):
    spec,d,p,g=inventory(folder);c=A.read_json(Path(folder)/'input.json')
    result=dict(complete=True,gate_pass=False,smoke_only=smoke,solver_validated=False)
    result['maxima']=dict(mach=float(g['mach'].max()),density_departure_percent=float(g['density_variation'].max()*100),
        mass_residual_relative=float(np.abs(g['mass_balance_relative'][1:]).max()))
    result['basic_pass']=all(result['maxima'][k]<A.LIMITS[k] for k in result['maxima'])
    if smoke:return result
    end=spec['duration_s'];T=A.period(c)
    left=end-.25 if spec['kind']=='healthy' else end-3*T
    result['window_s']=[left,end]
    result['stats']={k:A.window_stats(d['time_s'],d[k],left,end) for k in A.SCALARS}
    ref=A.healthy_references(c,.03)
    if spec['kind']=='healthy':
        values=result['stats']
        errors=dict(pressure=100*abs(values['dp20_50_pa']['mean']/ref['pressure_pa']-1),
            wss=100*abs(values['roi_mean_abs_wss_pa']['mean']/ref['wss_pa']-1),
            flow=100*abs(values['x030_flow_m2_s']['mean']/ref['flow_m2_s']-1),
            pressure_drift=100*(values['dp20_50_pa']['max']-values['dp20_50_pa']['min'])/ref['pressure_pa'])
        result.update(reference=ref,analytical_errors_percent=errors)
        result['gate_pass']=bool(result['basic_pass'] and errors['pressure']<5 and errors['wss']<5 and errors['flow']<1 and errors['pressure_drift']<1)
    else:
        result['repeatability']=A.repeatability(d,p,c,spec['cycles'])
        result['gate_pass']=bool(result['basic_pass'] and result['repeatability']['passed'])
    return result


def fixed_wall_fits(folder,out):
    """Free-intercept fit on smooth target wall, at fixed physical radii.

    This is a diagnostic measurement of a staircase-wall solution, not a new
    wall treatment. Invalid estimates are retained as invalid, never zeroed.
    """
    folder=Path(folder);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    c=A.read_json(folder/'input.json');index=A.read_json(folder/'diagnostics/phase_fields.json')
    with np.load(folder/'diagnostics/wall_geometry.npz') as z:
        x=z['x_mm']/1000;low=z['lower_mm']/1000;high=z['upper_mm']/1000;solid=z['solid'].copy()
    with np.load(folder/'diagnostics'/index[0]['file']) as z:y=z['y_m'].copy()
    dx=x[1]-x[0];stencils=[]
    for wall,surface in [('lower',low),('upper',high)]:
        slope=np.gradient(surface,dx)
        for xmm in np.arange(20.,55.00001,.5):
            j=int(np.argmin(abs(x-xmm/1000)))
            for radius in RADII_MM:
                st=build_stencil(x,y,solid,x[j],surface[j],slope[j],wall,radius/1000/dx)
                stencils.append((wall,xmm,j,radius,st))
    rows=[];sections=[]
    for entry in index:
        with np.load(folder/'diagnostics'/entry['file']) as archive:
            # NpzFile does not cache decoded arrays. Materialize once per frame,
            # not once per stencil or cross-section sample.
            z={name:archive[name] for name in ('time_s','u_m_s','v_m_s','lower_wss_pa','upper_wss_pa','pressure_pa','density_kg_m3')}
            for wall,xmm,j,radius,st in stencils:
                fit=st.evaluate(z['u_m_s'],z['v_m_s'],c['viscosity'])
                rows.append(dict(requested_time_s=entry['requested_time_s'],time_s=float(z['time_s']),wall=wall,x_mm=xmm,
                    actual_x_mm=x[j]*1000,radius_mm=radius,radius_cells=radius/1000/dx,
                    legacy_pa=float(z['lower_wss_pa' if wall=='lower' else 'upper_wss_pa'][j]),**fit))
            for xmm in sorted(set(range(5,int(round(c['length_mm'])),5))|{c['length_mm']-1}):
                j=int(np.argmin(abs(x-xmm/1000)))
                for i in np.flatnonzero(~solid[:,j]):
                    sections.append(dict(requested_time_s=entry['requested_time_s'],time_s=float(z['time_s']),
                        x_mm=xmm,actual_x_mm=x[j]*1000,y_mm=y[i]*1000,u_m_s=float(z['u_m_s'][i,j]),
                        v_m_s=float(z['v_m_s'][i,j]),pressure_pa=float(z['pressure_pa'][i,j]),density_kg_m3=float(z['density_kg_m3'][i,j])))
    write_csv(out/'fixed_radius_wall_fits.csv',rows);write_csv(out/'phase_cross_sections.csv',sections)
    valid=[r for r in rows if r['valid']]
    result=dict(radii_mm=RADII_MM,rows=len(rows),invalid_rows=len(rows)-len(valid),
        interpretation='Measurement diagnostic only; does not change solver or original acceptance scores.')
    if valid:
        result['max_wall_velocity_residual_m_s']=max(math.hypot(r['wall_tangent_velocity_m_s'],r['wall_normal_velocity_m_s']) for r in valid)
    if c['geometry_mode']=='healthy':
        expected=A.healthy_references(c,.03)['wss_pa'];last=index[-1]['requested_time_s']
        final=[r for r in rows if r['requested_time_s']==last]
        result['final_healthy_rms_error_percent']={str(radius):float(100*np.sqrt(np.mean([(r['shear_pa']/expected-1)**2 for r in final if r['radius_mm']==radius]))) for radius in RADII_MM} if all(r['valid'] for r in final) else None
    A.save_json(out/'fixed_radius_wall_fits.json',result)
    return result


def scalar_comparison(da,db,pa,pb,left,right):
    """B relative to A, retaining the original acceptance definitions."""
    rows=[]
    for key,stat,limit in [('dp20_50_pa','mean',2.),('roi_mean_abs_wss_pa','mean',5.),
            ('roi_peak_abs_wss_pa','max',10.),('x030_flow_m2_s','mean',1.)]:
        a=A.window_stats(da['time_s'],da[key],left,right)[stat];b=A.window_stats(db['time_s'],db[key],left,right)[stat]
        if abs(a)<1e-25:raise A.InvalidComparison('Unresolved reference mean/peak')
        value=100*abs(b-a)/abs(a)
        rows.append(dict(metric=key,statistic=stat,reference=a,comparison=b,value_percent=value,limit_percent=limit,passed=value<limit))
    for key,limit in [('dp20_50_pa',5.),('x030_flow_m2_s',2.)]:
        value=A.relative_l2(db['time_s'],db[key],da['time_s'],da[key],left,right)
        rows.append(dict(metric=key,statistic='waveform_l2',value_percent=value,limit_percent=limit,passed=value<limit))
    value=A.profiles_l2(pb,pa,left,right)
    rows.append(dict(metric='signed_wall_profiles',statistic='space_time_l2',value_percent=value,limit_percent=10.,passed=value<10))
    return rows


def interpolation_map(a,b):
    """B interpolated onto A; only common fluid nodes with four fluid corners."""
    xa,ya=a['x_m'],a['y_m'];xb,yb=b['x_m'],b['y_m']
    ix=np.clip(np.searchsorted(xb,xa,side='right')-1,0,len(xb)-2)
    iy=np.clip(np.searchsorted(yb,ya,side='right')-1,0,len(yb)-2)
    X,Y=np.meshgrid(xa,ya);I,J=np.meshgrid(ix,iy)
    wx=(X-xb[I])/(xb[I+1]-xb[I]);wy=(Y-yb[J])/(yb[J+1]-yb[J])
    mask=(~a['solid']) & (X>=xb[0]) & (X<=xb[-1]) & (Y>=yb[0]) & (Y<=yb[-1])
    for dj,di in ((0,0),(0,1),(1,0),(1,1)):mask &= ~b['solid'][J+dj,I+di]
    def interp(values):return ((1-wy)*((1-wx)*values[J,I]+wx*values[J,I+1])+wy*((1-wx)*values[J+1,I]+wx*values[J+1,I+1]))
    return X,mask,interp


def phase_pairs(a,b):
    sa=A.read_json(a/'case_spec.json');T=A.period(A.read_json(a/'input.json'));end=sa['duration_s']
    ai=A.read_json(a/'diagnostics/phase_fields.json');bi=A.read_json(b/'diagnostics/phase_fields.json')
    lookup={round(e['requested_time_s'],9):e for e in bi}
    pairs=[]
    for entry in ai:
        t=entry['requested_time_s']
        if not end-T-1e-9<=t<end-1e-9:continue
        other=lookup.get(round(t,9))
        if other is not None:pairs.append((entry,other))
    if len(pairs)<4:raise A.InvalidComparison('Need at least four common final-cycle phase fields')
    return pairs


def fields_pair(a,b):
    rows=[]
    for ea,eb in phase_pairs(a,b):
        with np.load(a/'diagnostics'/ea['file']) as za,np.load(b/'diagnostics'/eb['file']) as zb:
            keys=('x_m','y_m','solid','u_m_s','v_m_s','time_s','dt_s')
            fa={k:za[k] for k in keys};fb={k:zb[k] for k in keys}
            if abs(float(fa['time_s'])-float(fb['time_s']))>max(float(fa['dt_s']),float(fb['dt_s']))+1e-10:
                raise A.InvalidComparison('Phase fields not time matched')
            X,valid,interp=interpolation_map(fa,fb)
            bu,bv=interp(fb['u_m_s']),interp(fb['v_m_s'])
            for region,(lo,hi) in REGIONS.items():
                mask=valid&(X>=lo/1000-1e-12)&(X<=hi/1000+1e-12)
                if not mask.any():raise A.InvalidComparison('No common-fluid samples in region')
                au,av=fa['u_m_s'][mask],fa['v_m_s'][mask]
                den=np.sum(au**2+av**2)
                if den<=1e-25:raise A.InvalidComparison('Unresolved field reference norm')
                value=100*np.sqrt(np.sum((bu[mask]-au)**2+(bv[mask]-av)**2)/den)
                rows.append(dict(requested_time_s=ea['requested_time_s'],region=region,common_nodes=int(mask.sum()),
                    velocity_l2_percent=float(value),passed=bool(value<FIELD_LIMIT_PERCENT)))
    return rows


def compare_fits(a,b,analysis_a,analysis_b):
    pairs=phase_pairs(a,b);times={round(ea['requested_time_s'],9) for ea,eb in pairs}
    def load(path):
        with path.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
    def key(r):return (round(float(r['requested_time_s']),9),r['wall'],float(r['x_mm']),float(r['radius_mm']))
    ra={key(r):r for r in load(analysis_a/'fixed_radius_wall_fits.csv') if key(r)[0] in times}
    rb={key(r):r for r in load(analysis_b/'fixed_radius_wall_fits.csv') if key(r)[0] in times}
    if ra.keys()!=rb.keys():raise A.InvalidComparison('Fixed-radius wall-fit locations differ')
    rows=[]
    for radius in RADII_MM:
        for region,(lo,hi) in REGIONS.items():
            keys=[k for k in ra if k[3]==radius and lo<=k[2]<=hi]
            valid=all(ra[k]['valid']=='True' and rb[k]['valid']=='True' for k in keys)
            val=None
            if valid and keys:
                aa=np.array([float(ra[k]['shear_pa']) for k in keys]);bb=np.array([float(rb[k]['shear_pa']) for k in keys])
                den=np.linalg.norm(aa)
                val=float(100*np.linalg.norm(bb-aa)/den) if den>1e-25 else None
            rows.append(dict(radius_mm=radius,region=region,samples=len(keys),valid=valid,l2_percent=val))
    return rows


def compare(a,b,kind,summary_a,summary_b,analysis_a,analysis_b):
    a,b=Path(a),Path(b)
    ca,cb=A.read_json(a/'input.json'),A.read_json(b/'input.json')
    sa,sb=A.read_json(a/'case_spec.json'),A.read_json(b/'case_spec.json')
    if sa['epsilon']!=sb['epsilon'] or sa['duration_s']!=sb['duration_s']:
        raise A.InvalidComparison('Forcing or duration differs')
    A.matches_input(ca,cb,kind,A.read_json(a/'numerics.json'),A.read_json(b/'numerics.json'))
    end=sa['duration_s'];left=end-3*A.period(ca)
    da,db=[A.load(p/'diagnostics/fixed_probes.csv') for p in (a,b)]
    pa,pb=[A.load(p/'diagnostics/signed_wall_profiles.csv') for p in (a,b)]
    rows=scalar_comparison(da,db,pa,pb,left,end)
    fields=fields_pair(a,b)
    centroids=[]
    for folder in (a,b):
        jet=A.load(folder/'jet_probes.csv');centroids.append(A.window_stats(jet['time_s'],jet['x040_forward_centroid_y_mm'],end-A.period(ca),end)['mean'])
    same_branch=not (centroids[0]*centroids[1]<0 and min(abs(x) for x in centroids)>.05)
    passed=all(r['passed'] for r in rows) and all(r['passed'] for r in fields if r['region'] in ('all','downstream')) and same_branch and summary_a['gate_pass'] and summary_b['gate_pass']
    return dict(reference=a.parent.name,comparison=b.parent.name,kind=kind,window_s=[left,end],metrics=rows,
        field_comparisons=fields,field_limit_percent=FIELD_LIMIT_PERCENT,
        branch_centroids_x40_mm=centroids,same_branch_screen=same_branch,
        both_case_gates_pass=summary_a['gate_pass'] and summary_b['gate_pass'],passed=bool(passed),
        wall_fit_diagnostics=compare_fits(a,b,Path(analysis_a),Path(analysis_b)),
        notes='Direct same-orientation comparisons. Wall fits and wall swaps do not replace original scores.')


def domain_pass(comparisons):
    rows=[r for r in comparisons if r['kind']=='domain']
    return len(rows)==2 and all(r['passed'] for r in rows)
