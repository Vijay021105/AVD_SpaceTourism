"""Computed Euler snapshots with applicability/convergence gates and grid sensitivity."""
import numpy as np
from .cfd import solve,primitive
from .io import write_csv


def run(c,ascent,out,profile='standard',verify=True):
    c['cfd']=dict(gamma=1.4,upstream_m=10,outer_radius_m=4*c['diameter_m'],
        cfl=.5,local_cfl=.85,flux_scheme='rusanov',reconstruction='first_order',consecutive_checkpoints=3,
        residual_tolerance=c['aerodynamics']['cfd_residual_tolerance'],
        drag_change_tolerance=c['aerodynamics']['cfd_drag_change_tolerance'])
    nx,nr,it=(60,16,1500) if profile=='quick' else (100,28,12000)
    cases=[];solutions=[]
    for mach in [.3,.8,.98,1.2,2.,3.,5.,8.,15.,25.]:
        row=min(ascent,key=lambda r:abs(r['Mach']-mach));m=row['Mach'];totalT=row['T_K']*(1+.2*m*m)
        kn=row['mean_free_path_m']/c['diameter_m'];reason=[]
        if abs(m-mach)/mach>.12:reason.append('Mach_not_reached')
        if totalT>c['aerodynamics']['cfd_max_total_temp_K']:reason.append('perfect_gas_temperature_limit')
        if kn>.01:reason.append('rarefied_flow')
        if not row['atmosphere_valid']:reason.append('upper_atmosphere_extrapolated')
        label=f'M{mach:05.2f}_S{row["stage"]}';case=dict(case=label,requested_Mach=mach,actual_Mach=m,stage=row['stage'],altitude_m=row['h_m'],
            total_temperature_K=totalT,Kn_D=kn,model='axisymmetric perfect-gas Euler; zero incidence; pressure only')
        if reason:
            case.update(status='EXCLUDED',reason=';'.join(reason));cases.append(case);continue
        print('CFD:',label,flush=True)
        result=solve(c,row,nx,nr,it);summary=result['summary'];case.update(summary)
        case.update(status='CONVERGED' if summary['converged'] else 'NOT_CONVERGED',reason='')
        folder=out/'fields';folder.mkdir(exist_ok=True)
        rho,u,v,p=primitive(result['U']);mesh=result['mesh']
        np.savez_compressed(folder/f'{label}.npz',U=result['U'],X=mesh.X,R=mesh.R,xc=mesh.xc,rc=mesh.rc,
                            rho=rho,u=u,v=v,p=p,Mach=result['Mach'],freestream_Mach=m)
        write_csv(out/f'data/cfd/{label}_residual.csv',result['history']);write_csv(out/f'data/cfd/{label}_surface.csv',result['surface'])
        solutions.append((label,row,result));cases.append(case)
    write_csv(out/'data/cfd_cases.csv',cases)
    study=[]
    if verify:
        admissible=[r for _,r,_ in solutions]
        if admissible:
            row=min(admissible,key=lambda r:abs(r['Mach']-2))
            for nxx,nrr in ([(40,12),(60,16),(90,24)] if profile=='quick' else [(60,16),(100,28),(160,44)]):
                matching=next((res for _,rr,res in solutions if rr is row and res['summary']['Nx']==nxx and res['summary']['Nr']==nrr),None)
                res=matching or solve(c,row,nxx,nrr,it)
                record=dict(Mach=row['Mach'],**res['summary']);study.append(record)
                write_csv(out/f'data/cfd/grid_{nxx}_{nrr}_residual.csv',res['history'])
            write_csv(out/'data/cfd_grid_sensitivity.csv',study)
    # Pressure-vs-incidence diagnostic is explicitly analytical, never an axisymmetric CFD result.
    cp=[]
    for m in [.3,.98,2,5]:
        for side in ['windward','leeward']:
            sign=1 if side=='windward' else -1
            for x in np.linspace(0,1,100):
                local=np.deg2rad(15*(1-x)+sign*c['aerodynamics']['alpha_reference_deg'])
                value=2*max(np.sin(local),0)**2 if m>=2 else 2*local/np.sqrt(max(.15,abs(1-m*m)))
                cp.append(dict(Mach=m,side=side,x_over_capsule_length=x,Cp=value,
                    model='analytic exposure proxy; modified Newtonian for M>=2; compressibility estimate below; not CFD incidence pressure'))
    write_csv(out/'data/windward_leeward_Cp_proxy.csv',cp)
    return cases,solutions,study,cp
