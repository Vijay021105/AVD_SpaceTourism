"""Redraw professional figures from published CSV/NPZ products; no solver rerun."""
from pathlib import Path
import argparse,json
import numpy as np
from types import SimpleNamespace
from tara.io import table
from tara import figures


def redraw(out,media=False):
    out=Path(out);data=out/'data';c=json.loads((out/'effective_configuration.json').read_text());s=json.loads((out/'mission_summary.json').read_text())
    read=lambda name:table(data/f'{name}.csv')
    ascent=read('ascent')
    for row in ascent:row['stage']=int(row['stage'])
    cases=read('cfd_cases') if (data/'cfd_cases.csv').exists() else [];solutions=[]
    for case in cases:
        if case['status']=='EXCLUDED':continue
        case['converged']=str(case['converged']).lower()=='true'
        row=min(ascent,key=lambda r:abs(r['Mach']-case['actual_Mach']));label=case['case']
        with np.load(out/'fields'/f'{label}.npz') as field:
            mesh=SimpleNamespace(X=field['X'],R=field['R'],xc=field['xc'],rc=field['rc'],x=field['X'][:,0])
            result=dict(mesh=mesh,U=field['U'],p=field['p'],Mach=field['Mach'],
                history=table(data/'cfd'/f'{label}_residual.csv'),surface=table(data/'cfd'/f'{label}_surface.csv'),summary=case)
        solutions.append((label,row,result))
    study=read('cfd_grid_sensitivity') if (data/'cfd_grid_sensitivity.csv').exists() else []
    cp=read('windward_leeward_Cp_proxy') if cases else []
    flowdata=(cases,solutions,study,cp)
    modal=read('modal_frequencies');wet=[r for r in modal if r['state']=='wet']
    dyn=dict(wet=dict(frequency_Hz=np.array([r['frequency_Hz'] for r in wet])),shapes=read('mode_shapes'),
        spectral=read('acoustic_buffet_response'),pogo=read('pogo_frequency_response'),pogostat=json.loads((data/'pogo_assessment.json').read_text()),
        ground=read('ground_wind_loads'),slosh=read('slosh_modes'),separation=read('separation_clearance'),
        redundancy=read('redundancy_assessment'),pogo_transient=read('pogo_transient'))
    rv=read('rendezvous');rvdata=dict(history=rv,approach=[r for r in rv if r['phase']=='prescribed_docking_approach'],aerodb=read('aerodynamic_database'),polars=read('aero_polars'))
    th=(read('thermal_history'),read('TPS_history'),read('capsule_heat_distribution'),read('cryo_heat_leak_and_vents'),read('interstage_vents'),s['thermal'])
    trade=read('structural_fraction_trade')
    for row in trade:row['feasible']=str(row['feasible']).lower()=='true'
    names=figures.make(c,out,ascent,read('events'),s['stages'],trade,read('engine_layout'),read('loads_envelope'),read('structural_margins'),read('dry_mass_closure'),dyn,th,rvdata,flowdata)
    if media:figures.animations(c,out,ascent,rvdata,flowdata)
    print(f'Redrawn {len(names)} figures from stored numeric products.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path,nargs='?',default=Path(__file__).parent/'results');p.add_argument('--animations',action='store_true');a=p.parse_args();redraw(a.output,a.animations)
