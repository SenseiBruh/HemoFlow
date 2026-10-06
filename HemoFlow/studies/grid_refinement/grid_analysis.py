"""Controlled grid comparison. No reflection or filtering changes acceptance."""
import math
from pathlib import Path
from case_runtime import read,write
import sequence_analysis as S
import settling_analysis as C


def allow_fine(comparisons):
    return any(c.get('pair')=='severe_N080__severe_N120' and c.get('passed') is True and not c.get('smoke_only',False) for c in comparisons)


def compare(a,b,out,label,smoke=False):
    a,b,out=Path(a),Path(b),Path(out)
    dest=out/label;dest.mkdir(exist_ok=True)
    ca,cb=read(a/'input.json'),read(b/'input.json');sa,sb=read(a/'case_spec.json'),read(b/'case_spec.json')
    na,nb=read(a/'numerics.json'),read(b/'numerics.json')
    S.A.matches_input(ca,cb,'grid',na,nb)
    if ca['compute_threads']!=cb['compute_threads']:raise ValueError('Thread count changed')
    if sa['epsilon']!=sb['epsilon'] or sa['duration_s']!=sb['duration_s']:raise ValueError('Input or duration mismatch')
    result=dict(pair=label,smoke_only=smoke,passed=False,solver_validated=False,
        numerical_speed_m_s=[na['dx_m']/na['dt_s'],nb['dx_m']/nb['dt_s']],
        geometry_reports=[read(a/'geometry_report.json'),read(b/'geometry_report.json')])
    if smoke:
        result['note']='Short plumbing test; no scientific gate verdict.'
    else:
        aa,ab=out/sa['key'],out/sb['key']
        ra,rb=read(aa/'analysis.json'),read(ab/'analysis.json')
        ga=dict(gate_pass=ra['basic_pass'] and ra['cycle_settling']['recorded_observables_repeatable'])
        gb=dict(gate_pass=rb['basic_pass'] and rb['cycle_settling']['recorded_observables_repeatable'])
        # Existing full-field interpolation excludes nodes whose interpolation
        # stencil crosses a solid wall. WSS remains a separate mandatory test.
        old=S.REGIONS
        S.REGIONS={**old,'far_downstream':(50,55)}
        try:
            detail=S.compare(a,b,'grid',ga,gb,aa,ab)
        finally:S.REGIONS=old
        cycles=C.compare(a,b,dest,kind='grid')
        field_pass=all(r['passed'] for r in detail['field_comparisons'])
        result.update(aggregate=detail,cycle_comparison=cycles,
            passed=bool(detail['passed'] and field_pass and cycles['last_three_cycle_domain_screen_passed']),
            interpolation_note='Fine velocity interpolated to coarse physical nodes with four fluid corners. Excluded near-wall nodes are not verified by this field score; signed WSS is checked separately.',
            scope='Finite-resolution sensitivity screen only. No observed order/GCI or physical validation claim; timestep sensitivity is still outstanding.')
    write(dest/'grid_comparison.json',result)
    return result
