"""Initialization-only interventions. No state projection during evolution."""
import hashlib
import numpy as np
from experiment import reflected_populations

MODES=('native','symmetric','symmetric_uniform_density')

def initialize(s,mode):
    if mode not in MODES:raise ValueError('Unknown initialization mode')
    before=s.f.copy()
    if mode!='native':
        s.f=.5*(s.f+reflected_populations(s.f,s.ny,s.nx))
        if mode=='symmetric_uniform_density':
            rho=s.f.sum(axis=1).reshape(s.ny,s.nx)
            rho=.5*(rho+rho[::-1])
            s.f=s.f/rho.ravel()[:,None]
        # Make reported macroscopic moments consistent with the altered initial populations.
        from solver import CX,CY
        s.rho=s.f.sum(axis=1)
        s.ux=(s.f@CX)/s.rho;s.uy=(s.f@CY)/s.rho
        s.other=s.f.copy()
    fluid=s.nodes;mirror=reflected_populations(s.f,s.ny,s.nx)
    return dict(mode=mode,applied_at_time_s=0.,during_evolution=False,
        description='Native; or arithmetic population/reflection average; optionally normalize each node to unit lattice density while preserving its velocity.',
        population_relative_change=float(np.linalg.norm((s.f-before)[fluid])/np.linalg.norm(before[fluid])),
        max_population_mirror_error=float(np.max(abs((s.f-mirror)[fluid]))),
        density_min=float(s.rho[fluid].min()),density_max=float(s.rho[fluid].max()),
        initial_populations_sha256=hashlib.sha256(s.f[fluid].tobytes()).hexdigest())

def sample_time(index):
    return index*.00005 if index<=400 else .020+(index-400)*.00025

def symmetry_row(s,d):
    result={'time_s':s.time}
    for x in sorted(set(range(5,int(round(s.config.length_mm)),5))|{s.config.length_mm-10,s.config.length_mm-5,s.config.length_mm-1}):
        j=int(round(x/1000/s.dx));mask=s.geometry.fluid[:,j]
        u=d['u'][:,j];v=d['v'][:,j]
        result[f'x{int(x):03}_asymmetry_rms_m_s']=float(np.sqrt(np.mean(((u-u[::-1])**2+(v+v[::-1])**2)[mask])))
    return result
