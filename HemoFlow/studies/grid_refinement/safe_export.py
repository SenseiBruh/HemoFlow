"""Unique exports never overwrite an open ZIP. Export failure is nonfatal."""
from pathlib import Path
import os
import zipfile
import uuid
from datetime import datetime,timezone

CSV_NAMES={'cycle_metrics.csv','cycle_repeatability.csv','domain_cycles.csv','signed_wall_profiles.csv','global_timeseries.csv','downstream_wall_profiles.csv','fixed_probes.csv','jet_probes.csv','symmetry_probes.csv','jet_cycles.csv','phase_cross_sections.csv','fixed_radius_wall_fits.csv','startup_windows.csv','startup_onsets.csv'}

def build_archive(base,target):
    # The timestamp/UUID target has never existed. Avoid rename-overwrite entirely.
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(Path(base).rglob('*')):
            if not p.is_file() or 'exports' in p.relative_to(base).parts:continue
            if p.suffix in ('.json','.py','.md','.txt','.cmd','.log') or p.name in CSV_NAMES or p.name.endswith(('_reflection.csv','_startup_windows.csv')):
                z.write(p,p.relative_to(base))

def export_review(base,builder=build_archive):
    base=Path(base)
    try:
        out=base/'exports';out.mkdir(exist_ok=True)
        token=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')+'_'+uuid.uuid4().hex[:8]
        target=out/f'UPLOAD_GRID_REVIEW_{token}.zip'
        builder(base,target)
        with zipfile.ZipFile(target) as check:
            damaged=check.testzip()
            if damaged is not None:raise OSError('Export integrity check failed: '+damaged)
        summary=out/f'UPLOAD_GRID_SUMMARY_{token}.zip'
        with zipfile.ZipFile(summary,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for f in sorted(base.rglob('*')):
                rel=f.relative_to(base)
                if not f.is_file() or rel.parts[0] in ('exports','cases'):continue
                if f.suffix in ('.json','.py','.md','.txt','.cmd','.log') or f.name in ('cycle_metrics.csv','cycle_repeatability.csv','domain_cycles.csv'):
                    z.write(f,rel)
        with zipfile.ZipFile(summary) as check:
            if check.testzip() is not None:raise OSError('Summary integrity check failed')
        print('UPLOAD THIS: '+str(target),flush=True)
        print('If too large, upload summary first: '+str(summary),flush=True)
        return str(target)
    except Exception as exc:
        message='Export failed; simulation results remain on disk. Retry with --export-only. '+repr(exc)
        print(message,flush=True)
        try:
            with (base/'export_warnings.log').open('a',encoding='utf-8') as f:f.write(message+'\n')
        except OSError:pass
        return None
