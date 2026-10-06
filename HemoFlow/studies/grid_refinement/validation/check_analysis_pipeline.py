from pathlib import Path
import sys,tempfile,json,csv
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import asdict
from experiment import Config
from case_runtime import write
import grid_analysis as G
from strict_analysis import SCALARS,InvalidComparison

def csvsave(path,d):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(d);w.writerows(zip(*d.values()))
with tempfile.TemporaryDirectory() as temp:
    root=Path(temp);out=root/'analysis';out.mkdir(); folders=[]
    for n in (80,120):
        key=f'severe_N{n:03}';p=root/key;p.mkdir();(p/'diagnostics').mkdir();a=out/key;a.mkdir();folders.append(p)
        write(p/'case_spec.json',dict(key=key,epsilon=.001,duration_s=4.,cycles=4,cells=n,length_mm=150))
        write(p/'input.json',asdict(Config(length_mm=150.,cells_across=n,heart_rate=60)))
        write(p/'numerics.json',dict(dx_m=.004/n,dt_s=.3e-6*80/n))
        write(p/'geometry_report.json',dict(note='manufactured analysis fixture, not solver data'))
        write(a/'analysis.json',dict(basic_pass=True,cycle_settling=dict(recorded_observables_repeatable=True)))
        t=np.linspace(0,4,401);y=2+np.sin(2*np.pi*t)
        csvsave(p/'diagnostics/fixed_probes.csv',{'time_s':t,**{k:y for k in SCALARS}})
        csvsave(p/'diagnostics/signed_wall_profiles.csv',{'time_s':t,**{f'{w}_x{x}_pa':y for w in ('lower','upper') for x in (20,21,28,30,32,38,40,50,55)}})
        csvsave(p/'jet_probes.csv',{'time_s':t,'x040_forward_centroid_y_mm':t*0+.1})
        fields=[]
        for j,tt in enumerate((3.,3.25,3.5,3.75)):
            name=f'phase{j}.npz';x=np.linspace(.02,.055,71 if n==80 else 141);yy=np.linspace(-.001,.001,5 if n==80 else 9);xx,yy=np.meshgrid(x,yy)
            np.savez(p/'diagnostics'/name,x_m=x,y_m=yy[:,0],solid=np.zeros(xx.shape,bool),u_m_s=1+xx+yy,v_m_s=xx*0,time_s=tt,dt_s=.3e-6*80/n)
            fields.append(dict(file=name,requested_time_s=tt))
        write(p/'diagnostics/phase_fields.json',fields)
        (a/'fixed_radius_wall_fits.csv').write_text('requested_time_s,wall,x_mm,radius_mm,valid,shear_pa\n')
    result=G.compare(*folders,out,'identical_manufactured')
    assert result['passed'],result
    p=folders[1]/'diagnostics/signed_wall_profiles.csv'
    rows=list(csv.reader(p.open()));indices=[i for i,k in enumerate(rows[0]) if '_x50_' in k or '_x55_' in k]
    for row in rows[1:]:
        for i in indices:row[i]=str(float(row[i])*1.3)
    with p.open('w',newline='') as f:csv.writer(f).writerows(rows)
    result=G.compare(*folders,out,'far_wall_fault')
    assert not result['passed']
    write(folders[1]/'numerics.json',dict(dx_m=.004/120,dt_s=.3e-6))
    try:G.compare(*folders,out,'wave_speed_fault')
    except InvalidComparison:pass
    else:raise AssertionError('wave speed mismatch accepted')
print('PASS: full non-smoke analysis pipeline; affine matched grids pass; far-wall fault fails; wave-speed mismatch rejected. Manufactured data only.')
