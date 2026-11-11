"""Actual-state coast, multi-revolution J2 shooting and prescribed final approach.
All impulses reference the achieved ascent state; BAS is never reset to the crew.
The final approach is an inverse-dynamics feasibility path, not autonomous GNC.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from .io import MU,RE,G0,write_csv
from .orbits import propagate,rtn,gravity,elements,groundtrack


def run(c,injection,out):
    m=c['mission'];ti=float(injection[6]);crew0=injection[:6].copy();R=RE+m['orbit_altitude_m'];period=2*np.pi*np.sqrt(R**3/MU)
    wait=m['coast_periods']*period;tdepart=ti+wait
    tswait=np.linspace(ti,tdepart,420);crew_wait=propagate(crew0,tswait)
    bas_depart=propagate(c['bas_initial_state'],[0,tdepart])[-1]
    depart=crew_wait[-1];basis=rtn(depart)
    rb=basis.T@bas_depart[:3];phi=np.arctan2(rb[1],rb[0])
    N=int(m['phasing_revolutions']);turns=N+m['phasing_fraction']
    transfer_time=(turns-phi/(2*np.pi))*period
    aphase=(MU*(transfer_time/turns/(2*np.pi))**2)**(1/3)
    rad=np.linalg.norm(depart[:3]);vt=np.sqrt(MU*(2/rad-1/aphase))
    guess=basis[:,1]*vt-depart[3:6]
    arrival_t=tdepart+transfer_time
    bas_arr=propagate(c['bas_initial_state'],[0,arrival_t])[-1]
    bb=rtn(bas_arr);offset=bb@np.array([0,-m['standoff_m'],0.])
    omega=np.cross(bas_arr[:3],bas_arr[3:])/np.linalg.norm(bas_arr[:3])**2
    target_r=bas_arr[:3]+offset;target_v=bas_arr[3:]+np.cross(omega,offset)
    iterations=[]
    def propagate_end(dv):
        return propagate(np.r_[depart[:3],depart[3:]+dv],[tdepart,arrival_t])[-1]
    def residual(dv):
        end=propagate_end(dv);error=end[:3]-target_r
        iterations.append(dict(evaluation=len(iterations),miss_m=np.linalg.norm(error),departure_dv_mps=np.linalg.norm(dv)))
        return error/1000
    fit=least_squares(residual,guess,x_scale=100,diff_step=1e-5,max_nfev=70,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    end=propagate_end(fit.x);arrival_dv=target_v-end[3:]
    tsh=np.linspace(tdepart,arrival_t,1600)
    phased=propagate(np.r_[depart[:3],depart[3:]+fit.x],tsh)
    # BAS is propagated once on the union time grid; no interpolation across a reset.
    tapp=np.linspace(arrival_t,arrival_t+m['approach_duration_s'],400)
    alltimes=np.unique(np.r_[0,tswait,tsh,tapp]);basall=propagate(c['bas_initial_state'],alltimes)
    blookup={float(t):y for t,y in zip(alltimes,basall)}
    history=[];approach=[];coastdv=[]
    def record(t,crew,phase,acc=0.):
        bas=blookup[float(t)];br=rtn(bas);dr=crew[:3]-bas[:3];dv=crew[3:]-bas[3:]
        om=np.cross(bas[:3],bas[3:])/np.linalg.norm(bas[:3])**2
        rr=br.T@dr;vv=br.T@(dv-np.cross(om,dr))
        lo,la=groundtrack(c,t,crew)
        return dict(t_s=t,phase=phase,relative_distance_m=np.linalg.norm(dr),relative_speed_mps=np.linalg.norm(vv),
            radial_m=rr[0],alongtrack_m=rr[1],crosstrack_m=rr[2],
            rx_m=crew[0],ry_m=crew[1],rz_m=crew[2],vx_mps=crew[3],vy_mps=crew[4],vz_mps=crew[5],
            BAS_rx_m=bas[0],BAS_ry_m=bas[1],BAS_rz_m=bas[2],BAS_vx_mps=bas[3],BAS_vy_mps=bas[4],BAS_vz_mps=bas[5],
            altitude_m=np.linalg.norm(crew[:3])-RE,lon_deg=lo,lat_deg=la,commanded_accel_mps2=acc)
    for t,y in zip(tswait,crew_wait):history.append(record(t,y,'two_period_coast'))
    for t,y in zip(tsh,phased):history.append(record(t,y,'J2_phasing'))
    # Smooth along-track path with zero RTN velocity at both ends; basis derivatives
    # and relative acceleration evaluated numerically to include J2 frame rotation.
    s=(tapp-arrival_t)/m['approach_duration_s'];rel=-m['standoff_m']*(1-3*s*s+2*s**3)
    basapp=np.array([blookup[float(t)] for t in tapp]);off=np.array([rtn(b)[:,1]*q for b,q in zip(basapp,rel)])
    off_v=np.gradient(off,tapp,axis=0,edge_order=2);off_a=np.gradient(off_v,tapp,axis=0,edge_order=2)
    rel_rate=-m['standoff_m']*(-6*s+6*s*s)/m['approach_duration_s']
    for k,(t,b) in enumerate(zip(tapp,basapp)):
        om=np.cross(b[:3],b[3:])/np.linalg.norm(b[:3])**2
        # Analytical velocity makes endpoint velocity match exactly.
        vo=np.cross(om,off[k])+rtn(b)[:,1]*rel_rate[k]
        yy=np.r_[b[:3]+off[k],b[3:]+vo]
        command=gravity(b[:3])+off_a[k]-gravity(yy[:3])
        rec=record(t,yy,'prescribed_docking_approach',np.linalg.norm(command));history.append(rec);approach.append(rec);coastdv.append(np.linalg.norm(command))
    approach_dv=float(np.trapezoid(coastdv,tapp))
    impulses=[dict(event='phasing_departure',t_s=tdepart,dvx_mps=fit.x[0],dvy_mps=fit.x[1],dvz_mps=fit.x[2],dv_mps=np.linalg.norm(fit.x)),
        dict(event='arrival_velocity_match',t_s=arrival_t,dvx_mps=arrival_dv[0],dvy_mps=arrival_dv[1],dvz_mps=arrival_dv[2],dv_mps=np.linalg.norm(arrival_dv))]
    total=float(np.linalg.norm(fit.x)+np.linalg.norm(arrival_dv)+approach_dv)
    fuel=m['payload_kg']*(1-np.exp(-total/(G0*m['spacecraft_isp_s'])))
    available=sum(p['mass_kg'] for p in c['payload'] if p['component']=='spacecraft_propellant')
    status=dict(nominal_period_s=period,two_period_wait_s=wait,actual_phase_at_departure_deg=np.rad2deg(phi),
        phasing_revolutions=N,phasing_fraction=m['phasing_fraction'],phasing_duration_s=transfer_time,departure_dv_mps=np.linalg.norm(fit.x),arrival_dv_mps=np.linalg.norm(arrival_dv),
        approach_dv_mps=approach_dv,total_spacecraft_dv_mps=total,dv_allocation_mps=m['spacecraft_dv_allocation_mps'],
        dv_reserve_mps=m['spacecraft_dv_allocation_mps']-total,required_propellant_kg=fuel,available_propellant_kg=available,
        arrival_position_error_m=np.linalg.norm(end[:3]-target_r),shooting_success=bool(fit.success and np.linalg.norm(end[:3]-target_r)<10),
        phasing_min_altitude_m=float(np.min(np.linalg.norm(phased[:,:3],axis=1))-RE),
        max_approach_thrust_N=max(coastdv)*m['payload_kg'],docking_time_s=tapp[-1],
        final_distance_m=history[-1]['relative_distance_m'],final_relative_speed_mps=history[-1]['relative_speed_mps'],
        docking_model='prescribed smooth relative path; inverse-dynamics acceleration screen; no sensor/contact/GNC simulation',
        kinematic_mission_budget_pass=False)
    write_csv(out/'data/rendezvous.csv',history);write_csv(out/'data/rendezvous_burns.csv',impulses)
    write_csv(out/'data/rendezvous_targeting_iterations.csv',iterations)
    write_csv(out/'data/BAS_ephemeris.csv',[dict(t_s=t,**dict(zip(['rx_m','ry_m','rz_m','vx_mps','vy_mps','vz_mps'],y)),**elements(y)) for t,y in zip(alltimes,basall)])
    write_csv(out/'data/injection_state.csv',[dict(epoch_utc=m['epoch_utc'],t_since_epoch_s=ti,frame='ECI; GMST spherical Earth convention',
        **dict(zip(['rx_m','ry_m','rz_m','vx_mps','vy_mps','vz_mps'],crew0)),**elements(crew0))])
    return history,status,impulses,approach
