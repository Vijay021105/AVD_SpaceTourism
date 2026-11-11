"""One dependency-ordered mission pipeline and visible acceptance gates."""
from pathlib import Path
import json
import platform
import time
import shutil
import hashlib
import numpy as np
import scipy
import matplotlib
from .io import load_inputs,write_csv,json_write,G0
from . import sizing,geometry,ascent,rendezvous,structures,dynamics,thermal,flow,figures


def physics_hashes():
    root=Path(__file__).parent
    names=['ascent.py','orbits.py','environment.py','aerodynamics.py','sizing.py','geometry.py','structures.py']
    return {name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in names}


def coupled_sizing(c,out):
    rows,status,hist,trade=sizing.run(c,out);floors=None;iterations=[]
    for k in range(35):
        if floors is not None:
            rows,status,hist=sizing.optimize(c,floors=floors)
            for s,r in zip(c['stages'],rows):s.update(r)
        geometry.build(c)
        items=structures.hardware_mass(c)
        new=np.array([sum(r['mass_kg'] for r in items if r['stage']==j) for j in [1,2,3]])
        budget=np.array([r['dry_kg'] for r in rows]);error=max(new-budget)
        iterations.append(dict(iteration=k,liftoff_mass_kg=rows[0]['ignition_mass_kg'],
            sizing_feasible=status['success'],max_hardware_over_budget_kg=max(0,error),
            stage1_hardware_kg=new[0],stage2_hardware_kg=new[1],stage3_hardware_kg=new[2],
            stage1_dry_budget_kg=budget[0],stage2_dry_budget_kg=budget[1],stage3_dry_budget_kg=budget[2]))
        print(f'Sizing coupling {k}: {rows[0]["ignition_mass_kg"]/1000:.2f} t, mass deficit {max(error,0):.2f} kg',flush=True)
        if error<.05:break
        floors=new+1.0 # one kg numerical closure buffer; physical uncertainty handled separately
    status['hardware_floor_converged']=bool(error<.05);status['coupling_iterations']=k+1
    geometry.build(c,out);write_csv(out/'data/stages.csv',rows)
    write_csv(out/'data/sizing_iterations.csv',hist);write_csv(out/'data/mass_geometry_coupling.csv',iterations)
    return rows,status,trade,iterations


def run(inputs,out,profile='standard',make_animations=True,run_cfd=True,verify_cfd=True,reuse_ascent=False):
    started=time.time();out=Path(out);out.mkdir(parents=True,exist_ok=True)
    marker=out/'.tara-generated-output'
    if any(out.iterdir()) and not marker.exists():
        raise ValueError('Output directory is not a recognized TARA_LV results folder; choose an empty directory')
    marker.write_text('TARA_LV owns data, figures, fields, gifs and videos in this directory.\n')
    # These folders are generated products. Clear media/flow products so skipped
    # or excluded cases from an older configuration cannot appear in a new run.
    for sub in ['figures','fields','gifs','videos']:
        if (out/sub).exists():shutil.rmtree(out/sub)
    if (out/'data/cfd').exists():shutil.rmtree(out/'data/cfd')
    for name in ['cfd_cases.csv','cfd_grid_sensitivity.csv','windward_leeward_Cp_proxy.csv']:
        (out/'data'/name).unlink(missing_ok=True)
    for sub in ['data','figures','fields','gifs','videos']:(out/sub).mkdir(exist_ok=True)
    c=load_inputs(inputs);stages,sizing_status,trade,coupling=coupled_sizing(c,out)
    layout,installation=geometry.engine_layout(c,out)
    aerodb,polars=__import__('tara.aerodynamics',fromlist=['database']).database(c,out)
    cache=out/'data/ascent_checkpoint.json'
    if reuse_ascent and cache.exists():
        saved=json.loads(cache.read_text())
        if saved['input_hashes']!=c['input_hashes']:raise ValueError('Ascent cache input hashes differ; rerun targeting')
        if saved.get('physics_hashes')!=physics_hashes():raise ValueError('Ascent physics code changed; rerun targeting')
        print('Reusing verified nominal guidance',flush=True)
        rows,events,inj,_=ascent.simulate(c,saved['guidance'],resolution=300)
        flight=saved['status'];guide=np.array(saved['guidance']);write_csv(out/'data/ascent.csv',rows);write_csv(out/'data/events.csv',events)
    else:
        print('Targeting 3-D ascent',flush=True)
        rows,events,inj,flight,guide=ascent.target(c,out,quick=(profile=='quick'))
        json_write(cache,dict(input_hashes=c['input_hashes'],physics_hashes=physics_hashes(),guidance=guide,status=flight,injection=inj))
    print('Injection:',flight['final_altitude_m']/1000,'km; target pass:',flight['orbit_target_pass'],flush=True)
    print('Solving actual-state J2 phasing',flush=True)
    rvhist,rvstatus,burns,approach=rendezvous.run(c,inj,out)
    print('Rendezvous delta-v:',rvstatus['total_spacecraft_dv_mps'],'m/s; miss:',rvstatus['arrival_position_error_m'],'m',flush=True)
    lastmass=rows[-1]['mass_kg'];unused=max(0,flight['remaining_stage3_propellant_kg']-c['stages'][2]['reserve_kg'])
    launcher_reserve=G0*c['stages'][2]['engine']['isp_vac_s']*np.log(lastmass/(lastmass-unused))
    write_csv(out/'data/velocity_reserves.csv',[
        dict(owner='launcher_S3',reserve_mps=launcher_reserve,propellant_kg=unused,note='unused launcher ideal capability at injection; stays with separated S3'),
        dict(owner='crew_spacecraft',reserve_mps=rvstatus['dv_reserve_mps'],propellant_kg=rvstatus['available_propellant_kg']-rvstatus['required_propellant_kg'],note='remaining allocated maneuver budget; payload propulsion only')])
    print('Computing flight loads, dynamics and thermal screens',flush=True)
    env,margins,masslist,closure,loadstatus=structures.run(c,rows,out)
    write_csv(out/'data/engine_out_thrust_screen.csv',[dict(t_s=r['t_s'],stage=r['stage'],
        nominal_thrust_N=r['thrust_N'],one_engine_unavailable_thrust_N=r['thrust_N']*(c['stages'][int(r['stage'])-1]['engine_count']-1)/c['stages'][int(r['stage'])-1]['engine_count'],
        mass_kg=r['mass_kg'],note='instantaneous thrust screen at nominal mass and throttle; no failure trajectory or control validation') for r in rows if r['throttle']>0])
    dyn=dynamics.run(c,rows,events,out);th=thermal.run(c,rows,out)
    flows=flow.run(c,rows,out,profile,verify_cfd) if run_cfd else ([],[],[],[])
    checks=[]
    def gate(name,passed,value,criterion,scope='conceptual numerical/model check'):
        checks.append(dict(check=name,status='PASS' if passed else 'OPEN',value=value,criterion=criterion,scope=scope))
    m=c['mission']
    gate('constrained_sizing',sizing_status['success'],sizing_status['achieved_dv_mps'],f'launcher ideal delta-v >= {sizing_status["target_launcher_dv_mps"]:.1f} m/s and selected ignition T/W')
    gate('dry_hardware_floor',all(r['closure_pass'] for r in closure),min(r['unallocated_dry_budget_kg'] for r in closure),'estimated hardware <= dry budgets; missing hardware still requires design')
    gate('tank_capacity',all(abs(t['capacity_mass_kg']-t['loaded_mass_kg'])<.01 for t in c['tanks']),max(abs(t['capacity_mass_kg']-t['loaded_mass_kg']) for t in c['tanks']),'mass-volume equality <0.01 kg')
    gate('diameter',max(t['installed_outer_radius_m'] for t in c['tanks'])*2<=m['diameter_m'],max(t['installed_outer_radius_m'] for t in c['tanks'])*2,'installed tanks inside 4 m envelope')
    gate('engine_installation',installation['gimbal_installation_pass'],installation['gimbal_envelope_diameter_m'],'4 m radial and pairwise gimbal clearance; assumed supplier dimensions')
    gate('orbit_injection',flight['orbit_target_pass'],flight['final_altitude_m'],'altitude ±2 km, speed ±2 m/s, radial speed <2 m/s, inclination ±0.1 deg')
    gate('ascent_integration_refinement',flight['integration_refinement_position_m']<1000 and flight['integration_refinement_velocity_mps']<2,
        flight['integration_refinement_position_m'],'doubled fraction resolution: <1000 m and <2 m/s')
    gate('dynamic_pressure',flight['max_q_Pa']<=m['q_limit_Pa']*1.01,flight['max_q_Pa'],f'<= {m["q_limit_Pa"]} Pa')
    gate('q_alpha',flight['max_qalpha_Pa_rad']<=m['qalpha_limit_Pa_rad']*1.01,flight['max_qalpha_Pa_rad'],f'<= {m["qalpha_limit_Pa_rad"]} Pa rad')
    gate('crew_acceleration',flight['max_proper_g']<=m['proper_g_limit']*1.01,flight['max_proper_g'],f'<= {m["proper_g_limit"]} g; physiological qualification not established')
    gate('phasing_shooting',rvstatus['shooting_success'],rvstatus['arrival_position_error_m'],'<10 m endpoint miss before velocity-match impulse')
    gate('phasing_orbit_clearance',rvstatus['phasing_min_altitude_m']>180000,rvstatus['phasing_min_altitude_m'],'screening perigee >180 km')
    gate('spacecraft_delta_v',rvstatus['total_spacecraft_dv_mps']<=m['spacecraft_dv_allocation_mps'],rvstatus['total_spacecraft_dv_mps'],'within allocated 300 m/s; no double counting launcher ideal reserve')
    gate('spacecraft_propellant',rvstatus['required_propellant_kg']<=rvstatus['available_propellant_kg'],rvstatus['required_propellant_kg'],'required propellant <= included payload propellant')
    gate('approach_thrust',rvstatus['max_approach_thrust_N']<=m['approach_thrust_available_N'],rvstatus['max_approach_thrust_N'],'prescribed-path demand <= assumed available thrust')
    finite_margins=[float(r[k]) for r in margins for k in ['yield_margin','ultimate_margin','buckling_margin'] if r.get(k,'')!='']
    gate('structural_screens',min(finite_margins)>=0,min(finite_margins),'all calculated metallic and buckling margins >=0; qualification work remains')
    gate('TPS_ascent',th[-1]['TPS_screen_pass'],th[-1]['max_bondline_K'],'bondline <= assumed 423 K for ascent only')
    gate('cryo_boiloff_reserve',th[-1]['cryo_reserve_pass'],max(t['boiloff_fraction'] for t in th[3]),'pad+ascent boiloff <= selected 1% liquid reserve; operational loading must account for losses')
    gate('POGO_assumed_loop',dyn['pogostat']['status']=='NO_RHP_ROOTS_DETECTED',dyn['pogostat']['status'],'no unstable roots in supplied scalar delay loop; not flight validation')
    if flows[0]:
        admissible=[r for r in flows[0] if r['status']!='EXCLUDED']
        gate('Euler_iterations',all(r['status']=='CONVERGED' for r in admissible),sum(r['status']=='CONVERGED' for r in admissible),f'{len(admissible)} admissible cases must meet residual + force criteria')
    gate('human_rating',False,'NOT_ESTABLISHED','hazard analysis, abort/escape, FDIR, ECLSS, qualification, fracture control, contact dynamics and test evidence required','program verification gate')
    gate('viscous_separated_CFD',False,'NOT_IMPLEMENTED','3-D transition/RANS/URANS, base/plume and incidence validation required','requested higher fidelity analysis')
    gate('aeroelastic_stability',False,'NOT_ESTABLISHED','two-way aerodynamic/structural/control coupling and flutter/divergence analysis required','requested higher fidelity analysis')
    rvstatus['kinematic_mission_budget_pass']=bool(flight['orbit_target_pass'] and rvstatus['shooting_success'] and rvstatus['total_spacecraft_dv_mps']<=m['spacecraft_dv_allocation_mps'] and rvstatus['phasing_min_altitude_m']>180000)
    write_csv(out/'data/design_checks.csv',checks)
    print('Rendering mission figures',flush=True)
    rvdata=dict(history=rvhist,approach=approach,aerodb=aerodb,polars=polars)
    names=figures.make(c,out,rows,events,stages,trade,layout,env,margins,closure,dyn,th,rvdata,flows)
    animation_status=figures.animations(c,out,rows,rvdata,flows) if make_animations else []
    summary=dict(project='TARA_LV',version='1.0.0',mission='Space tourism crew launch; two nominal periods; J2 phasing; docking approach',
        model_readiness='CONCEPTUAL_RESEARCH_NOT_HUMAN_RATED',input_hashes=c['input_hashes'],
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
        sizing=sizing_status,installation=installation,flight=flight,rendezvous=rvstatus,thermal=th[-1],dynamics=dyn['status'],loads=loadstatus,
        liftoff_mass_kg=stages[0]['ignition_mass_kg'],height_m=c['tip_z_m'],diameter_m=c['diameter_m'],
        stages=stages,checks=checks,figures=names,animations=animation_status,runtime_s=time.time()-started)
    json_write(out/'mission_summary.json',summary)
    json_write(out/'effective_configuration.json',c)
    from .report import write_report
    write_report(c,out,summary,closure,flows[0])
    print(f'Finished: {len(names)} figures; {sum(r["status"]=="PASS" for r in checks)}/{len(checks)} checks passed. See {out}/MISSION_REPORT.md',flush=True)
    return summary
