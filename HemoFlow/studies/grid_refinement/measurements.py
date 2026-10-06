"""Read-only physical probes. Unfiltered samples; no solver-state mutation."""
import csv
import json
import math
import os
from pathlib import Path
import numpy as np

PROBES_MM=(5,15,20,25,30,35,40,45,50,55)
PROFILE_MM=np.arange(20.,55.0001,.5)
METRICS=('dp20_50_pa','dp20_55_pa','dp30_50_pa','roi_mean_abs_wss_pa',
         'roi_peak_abs_wss_pa','x030_flow_m2_s','backflow40_55_nodes_percent')


def json_write(path,value):
    temp=Path(str(path)+'.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    os.replace(temp,path)


class Probes:
    def __init__(self,solver,folder,duration):
        self.folder=Path(folder) if folder else None
        self.rows=0;self.profile_rows=0;self.last=-1.;self.next_profile=0.
        self.handles=[];self.writers={};self.frames=[]
        self.x=solver.geometry.x*1000
        stations=sorted(set(PROBES_MM)|set(range(60,int(round(self.x[-1])),5)))
        self.targets=list(stations)+[self.x[-1]-1.,self.x[-1]-solver.dx*1000]
        self.labels=[f'x{int(x):03}' for x in stations]+['outlet_minus1mm','outlet_minus1cell']
        self.roi=(self.x>=20-1e-8)&(self.x<=55+1e-8)
        self.backroi=(self.x>=40-1e-8)&(self.x<=55+1e-8)
        self.weights=np.ones(self.roi.sum())*solver.dx
        self.weights[[0,-1]]*=.5
        self.wlow=self.weights*np.sqrt(1+solver.lower_slope[self.roi]**2)
        self.whigh=self.weights*np.sqrt(1+solver.upper_slope[self.roi]**2)
        self.period=60/solver.config.heart_rate
        count=int(round(duration/self.period))
        candidates=[0,.0003,.0006,.001,.002,.005,.01,.02,.025,.05,.075,.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.4,.5,.625,.75,self.period,1.,1.25,1.5]
        candidates += [(cycle+ph/8)*self.period for cycle in range(count) for ph in range(8)] if count>=2 else []
        # Deduplicate by solver step to avoid two requested frames on one iteration.
        unique={}
        for value in sorted(t for t in candidates+[duration] if 0<=t<=duration):
            unique.setdefault(math.ceil(value/solver.dt-1e-9),value)
        self.phase_targets=sorted(unique.values())
        self.next_frame=0
        if self.folder:
            self.folder.mkdir(parents=True,exist_ok=True)
            outlet_parameters={'treatment':'fixed-density regularized outlet'}
            if solver.config.boundary_model=='characteristic':
                from outlet_candidate import SIGMA,RELAX_LENGTH_IN_GAPS
                outlet_parameters=dict(treatment='experimental characteristic outlet',sigma=SIGMA,
                                       relaxation_length_in_nominal_gaps=RELAX_LENGTH_IN_GAPS,
                                       relaxation_length_m=RELAX_LENGTH_IN_GAPS*solver.config.diameter,
                                       incoming_transverse_reservoir_velocity_m_s=0.)
            elif solver.config.boundary_model=='section_impedance':
                outlet_parameters=dict(treatment='experimental section-mean acoustic impedance',
                    revision='2026-09-13.2',reference_mean_velocity_m_s=solver.config.mean_velocity,
                    acoustic_impedance_pa_s_m=solver.config.density*solver.velocity_scale/math.sqrt(3),
                    pressure_relation='p_out=rho_ref*c_s*(mean_u_adjacent_old-U_ref)',
                    velocity_and_stress='Same-time streamed interior non-equilibrium extrapolation; not clipped.',
                    limits='Plane-mode condition at a straight full-height outlet. Not a physiological load or proof of vortex/backflow accuracy.')
            location=[]
            for name,x in zip(self.labels,self.targets):
                right=int(np.searchsorted(self.x,x,side='right'));right=min(right,len(self.x)-1);left=max(0,right-1)
                w=(x-self.x[left])/(self.x[right]-self.x[left])
                location.append(dict(name=name,target_x_mm=x,left_x_mm=float(self.x[left]),right_x_mm=float(self.x[right]),right_weight=float(w)))
            json_write(self.folder/'definitions.json',dict(
                model='2D planar rigid-wall Newtonian; area labels are circular equivalents only',
                numerics=solver.numerics(),outlet_parameters=outlet_parameters,probe_locations=location,profile_x_mm=PROFILE_MM.tolist(),
                phase_targets_s=self.phase_targets,definitions={
                    'pressure_pa':'Arithmetic mean over fluid nodes of section, then linear interpolation in x; gauge pressure.',
                    'density_kg_m3':'Arithmetic cross-section mean lattice density times reference density.',
                    'flow_m2_s':'dx * sum(u) over section; volume flow per unit out-of-plane depth, NOT mL/min.',
                    'reverse_flow_m2_s':'Positive magnitude dx*sum(max(-u,0)); separate from node backflow percentage.',
                    'wall_shear_pa':'Existing approximate first-fluid-node wall-tangential estimate. Positive along increasing x at BOTH walls.',
                    'roi_mean_abs_wss_pa':'20-55 mm arclength-weighted mean absolute shear on both native walls.',
                    'roi_peak_abs_wss_pa':'Largest native-node absolute shear in 20-55 mm; sampling and wall-grid dependent.',
                    'backflow40_55_nodes_percent':'Fluid node fraction u < -0.01*mean inlet speed in 40-55 mm.',
                    'profiles':'Signed shear at 0.5 mm fixed locations, nominal 1 kHz; native wall geometry and phase fields also saved.',
                    'phase_fields':'Full unsmoothed float64 fields at first lattice step at/after requested phase; actual timing recorded.',
                    'sampling':'50 microseconds through 20 ms, then 250 microseconds. Actual times recorded. No filtering.'}))
            np.savez_compressed(self.folder/'wall_geometry.npz',x_mm=self.x,lower_mm=solver.geometry.lower*1000,
                                upper_mm=solver.geometry.upper*1000,solid=solver.geometry.solid)

    def values(self,s,d):
        g=s.geometry;mask=g.fluid;u=d['u'];rho=s.rho.reshape(mask.shape)*s.config.density
        columns=dict(pressure_pa=d['mean_p'],density_kg_m3=(rho*mask).sum(0)/mask.sum(0),
                     mean_u_m_s=d['mean_u'],flow_m2_s=(u*mask).sum(0)*s.dx,
                     mass_flow_kg_m_s=(rho*u*mask).sum(0)*s.dx,
                     forward_flow_m2_s=(np.maximum(u,0)*mask).sum(0)*s.dx,
                     reverse_flow_m2_s=(np.maximum(-u,0)*mask).sum(0)*s.dx,
                     lower_wss_pa=d['wss_low'],upper_wss_pa=d['wss_high'])
        row=dict(time_s=float(s.time),epoch=0)
        for key,values in columns.items():
            row.update({name+'_'+key:float(v) for name,v in zip(self.labels,np.interp(self.targets,self.x,values))})
        for a,b in ((20,50),(20,55),(30,50)):
            row[f'dp{a}_{b}_pa']=row[f'x{a:03}_pressure_pa']-row[f'x{b:03}_pressure_pa']
        low,high=d['wss_low'][self.roi],d['wss_high'][self.roi]
        row['roi_mean_abs_wss_pa']=float((abs(low)@self.wlow+abs(high)@self.whigh)/(self.wlow.sum()+self.whigh.sum()))
        row['roi_peak_abs_wss_pa']=float(max(abs(low).max(),abs(high).max()))
        region=mask[:,self.backroi]
        row['backflow40_55_nodes_percent']=float(100*(u[:,self.backroi][region]<-.01*s.config.mean_velocity).mean())
        ur=u[:,self.roi];valid=mask[:,self.roi]&mask[::-1,self.roi]
        row['axial_mirror_asymmetry_l2_percent']=float(100*np.linalg.norm((ur-ur[::-1])[valid])/max(np.linalg.norm(ur[valid]),1e-30))
        if not all(math.isfinite(v) for v in row.values()):
            raise ValueError('Nonfinite fixed probe; original samples retained, study stopped')
        return row

    def append(self,name,row):
        if not self.folder:return
        if name not in self.writers:
            f=(self.folder/name).open('w',newline='',encoding='utf-8')
            self.handles.append(f);w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();self.writers[name]=(w,f)
        w,f=self.writers[name];w.writerow(row);f.flush()

    def record(self,s,d):
        if s.time<=self.last:return
        self.append('fixed_probes.csv',self.values(s,d));self.rows+=1;self.last=s.time
        if s.time>=self.next_profile-1e-12:
            row=dict(time_s=float(s.time))
            for label,key in (('lower','wss_low'),('upper','wss_high')):
                for x,value in zip(PROFILE_MM,np.interp(PROFILE_MM,self.x,d[key])):
                    row[f'{label}_x{x:04.1f}_pa']=float(value)
            self.append('signed_wall_profiles.csv',row);self.profile_rows+=1
            tail=dict(time_s=float(s.time))
            for label,key in (('lower','wss_low'),('upper','wss_high')):
                for x in np.arange(55.,self.x[-1]-1.+1e-9,.5):
                    tail[f'{label}_x{x:04.1f}_pa']=float(np.interp(x,self.x,d[key]))
            self.append('downstream_wall_profiles.csv',tail)
            self.next_profile=(math.floor((s.time+1e-12)/.001)+1)*.001

    def phase_due(self,s):
        return self.next_frame<len(self.phase_targets) and s.time>=self.phase_targets[self.next_frame]-1e-12

    def snapshot(self,s):
        if not self.phase_due(s):return
        requested=self.phase_targets[self.next_frame]
        if self.folder:
            d=s.diagnostics();g=s.geometry
            name=f'phase_{self.next_frame:02d}.npz';path=self.folder/name
            with Path(str(path)+'.tmp').open('wb') as f:
                np.savez_compressed(f,x_m=g.x,y_m=g.y,solid=g.solid,u_m_s=d['u'],v_m_s=d['v'],
                                    pressure_pa=d['p'],density_kg_m3=s.rho.reshape(g.fluid.shape)*s.config.density,
                                    lower_wss_pa=d['wss_low'],upper_wss_pa=d['wss_high'],
                                    requested_time_s=requested,time_s=s.time,dt_s=s.dt,iteration=s.iteration,
                                    boundary_model=s.config.boundary_model)
            os.replace(str(path)+'.tmp',path)
            self.frames.append(dict(file=name,requested_time_s=requested,time_s=s.time,error_s=s.time-requested))
            json_write(self.folder/'phase_fields.json',self.frames)
        self.next_frame+=1

    def close(self):
        for f in self.handles:f.close()
        if self.folder:json_write(self.folder/'sample_status.json',dict(rows=self.rows,profile_rows=self.profile_rows,last_time_s=self.last,closed=True,phase_frames=len(self.frames)))
