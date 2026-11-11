"""Rotating-atmosphere, 3-D ECI point-mass ascent and numerical targeting.
RK4 in consumed-propellant fraction gives exact mass cutoffs. Final sampling
is finer than optimization. No claimed 6-DOF control or flight software.
"""
import numpy as np
from scipy.optimize import least_squares
from .io import G0,MU,RE,OMEGA,write_csv,json_write
from .orbits import gravity,frame,rotation_z,setup,elements,groundtrack
from .environment import atmosphere,wind
from .aerodynamics import coefficients


def seed(c):
    g=c['guidance']
    return np.array([g[k] for k in ['pitch1_deg','pitch2_deg','pitch3_deg','stage2_end_deg',
        'stage3_start_deg','stage3_end_deg','coast23_s','stage3_use_fraction','plane_trim_deg']])


def simulate(c,x=None,resolution=100,save=True):
    if x is None:x=seed(c)
    g=c['guidance'];mission=c['mission'];r,v,_=setup(c)
    # r,v,time,integrated thrust dv,drag loss,gravity projection loss.
    y=np.r_[r,v,0.,0.,0.,0.];rows=[];events=[];initial_mass=c['stages'][0]['ignition_mass_kg'];coast12=g['coast12_s']
    mass=initial_mass
    def evaluate(y,j,m,f,powered=True):
        r=y[:3];v=y[3:6];t=y[6];s=c['stages'][j];h=np.linalg.norm(r)-RE
        env=atmosphere(h);er,ee,en=frame(r);we,wn=wind(c,h)
        va=v-np.cross([0,0,OMEGA],r)-we*ee-wn*en;speed=np.linalg.norm(va);unit=va/max(speed,1e-9)
        M=speed/env['a_mps'];q=.5*env['rho_kgm3']*speed*speed
        trim=np.deg2rad(x[8]) if len(x)>8 else 0.
        axis=c['site_r_eci']/np.linalg.norm(c['site_r_eci'])
        command_normal=c['plane_normal']*np.cos(trim)+np.cross(axis,c['plane_normal'])*np.sin(trim)
        normal=rotation_z(c['node_precession_radps']*t)@command_normal
        et=np.cross(normal,er);et/=np.linalg.norm(et)
        if j==0:
            hold=min(.15,max(.01,g['vertical_time_s']*s['mdot_full_kgps']/s['usable_propellant_kg']))
            pitch=np.interp(f,[0,hold,g['pitch_fraction_1'],g['pitch_fraction_2'],1],[90,90,x[0],x[1],x[2]])
        elif j==1:pitch=np.interp(f,[0,1],[x[2],x[3]])
        else:pitch=np.interp(f,[0,x[7]],[x[4],x[5]])
        pr=np.deg2rad(pitch)
        ndot=c['node_precession_radps']*np.cross([0,0,1],normal)
        yaw=np.clip(-4e-6*np.dot(r,normal)-.004*(np.dot(v,normal)+np.dot(r,ndot)),-.25,.25)*np.clip((h-20000)/40000,0,1)
        direction=np.sin(pr)*er+np.cos(pr)*et+yaw*normal;direction/=np.linalg.norm(direction)
        raw=np.arccos(np.clip(np.dot(unit,direction),-1,1)) if speed>10 else 0
        # Vector load-relief steering bounds incidence at appreciable q.
        alpha_limit=min(np.deg2rad(7),mission['qalpha_limit_Pa_rad']/max(q,1))
        if q>1000 and raw>alpha_limit:
            transverse=direction-np.dot(direction,unit)*unit;transverse/=max(np.linalg.norm(transverse),1e-12)
            direction=np.cos(alpha_limit)*unit+np.sin(alpha_limit)*transverse
        alpha=np.arccos(np.clip(np.dot(unit,direction),-1,1)) if speed>10 else 0.
        ac=coefficients(c,M,alpha,j+1,env['rho_kgm3'],speed,env['mu_Pas'])
        area=np.pi*c['diameter_m']**2/4;drag=q*area*ac['Cd']
        eng=s['engine'];ae=(eng['thrust_vac_N']-eng['thrust_sl_N'])/101325 if eng['thrust_sl_N'] else np.pi*eng['nozzle_exit_diameter_m']**2/4
        full=max(0,s['engine_count']*(eng['thrust_vac_N']-env['p_Pa']*ae))
        th=min(1,(mission['proper_g_limit']*G0*m+drag)/max(full,1))
        if j==0 and q>mission['q_limit_Pa']*g['q_throttle_on_fraction']:
            th=min(th,1-g['q_controller_gain']*(q/(mission['q_limit_Pa']*g['q_throttle_on_fraction'])-1))
        th=float(np.clip(th,eng['min_throttle'],1)) if powered else 0.
        thrust=th*full;grav=gravity(r)
        # Normal force opposes lateral flow in the plane of incidence.
        transverse=direction-np.dot(direction,unit)*unit
        nu=transverse/max(np.linalg.norm(transverse),1e-12)
        normal_force=q*area*ac['Cn'];faero=-drag*unit+normal_force*nu
        acc=grav+(thrust*direction+faero)/m
        deriv=np.r_[v,acc,1,thrust/m,drag/m,-np.dot(grav,v/max(np.linalg.norm(v),1))]
        rec=dict(t_s=t,stage=j+1,h_m=h,speed_inertial_mps=np.linalg.norm(v),v_mps=speed,vr_mps=np.dot(v,er),
            mass_kg=m,Mach=M,q_Pa=q,alpha_deg=np.rad2deg(alpha),qalpha_Pa_rad=q*alpha,
            commanded_pitch_deg=pitch,pitch_deg=np.rad2deg(np.arcsin(np.clip(np.dot(direction,er),-1,1))),
            thrust_N=thrust,drag_N=drag,normal_N=normal_force,proper_g=np.linalg.norm(thrust*direction+faero)/(m*G0),
            proper_axial_mps2=np.dot(thrust*direction+faero,direction)/m,throttle=th,
            active_prop_fraction=f,integrated_thrust_dv_mps=y[7],drag_loss_mps=y[8],gravity_loss_projection_mps=y[9],
            Ca=ac['Ca'],Cd=ac['Cd'],Cn=ac['Cn'],Cnalpha_per_rad=ac['Cnalpha_per_rad'],
            Xcp_from_nose_m=ac['Xcp_from_nose_m'],east_wind_mps=we,north_wind_mps=wn,**env)
        rec.update(dict(zip(['rx_m','ry_m','rz_m','vx_mps','vy_mps','vz_mps'],y[:6])))
        rec['lon_deg'],rec['lat_deg']=groundtrack(c,t,y)
        return deriv,rec,th
    for j,s in enumerate(c['stages']):
        start_mass=mass;usable=s['usable_propellant_kg'];fend=x[7] if j==2 else 1.
        N=max(40,int(resolution*fend));fs=np.linspace(0,fend,N+1)
        events.append(dict(event=f'S{j+1}_ignition',t_s=y[6],mass_before_kg=mass,mass_after_kg=mass,
            rx_m=y[0],ry_m=y[1],rz_m=y[2],vx_mps=y[3],vy_mps=y[4],vz_mps=y[5]))
        def rhs(f,state):
            mm=start_mass-usable*f;dy,_,th=evaluate(state,j,mm,f)
            return dy*usable/(s['mdot_full_kgps']*max(th,1e-8))
        for f0,f1 in zip(fs[:-1],fs[1:]):
            if save:rows.append(evaluate(y,j,start_mass-usable*f0,f0)[1])
            df=f1-f0;k1=rhs(f0,y);k2=rhs(f0+df/2,y+k1*df/2);k3=rhs(f0+df/2,y+k2*df/2);k4=rhs(f1,y+k3*df)
            y=y+df*(k1+2*k2+2*k3+k4)/6
            if np.linalg.norm(y[:3])<RE-10000 or not np.isfinite(y).all():
                raise RuntimeError(f'Ground intersection S{j+1} f={f1:.3f} h={np.linalg.norm(y[:3])-RE:.0f} t={y[6]:.0f}')
        mass=start_mass-usable*fend
        if save:rows.append(evaluate(y,j,mass,fend)[1])
        event=dict(event=f'S{j+1}_cutoff',t_s=y[6],mass_before_kg=mass,mass_after_kg=mass,
            rx_m=y[0],ry_m=y[1],rz_m=y[2],vx_mps=y[3],vy_mps=y[4],vz_mps=y[5]);events.append(event)
        if j<2:
            mass-=s['dry_kg']+s['reserve_kg']
            events.append(dict(event=f'S{j+1}_separation',t_s=y[6],mass_before_kg=event['mass_before_kg'],mass_after_kg=mass,
                rx_m=y[0],ry_m=y[1],rz_m=y[2],vx_mps=y[3],vy_mps=y[4],vz_mps=y[5]))
            coast=coast12 if j==0 else x[6]
            nt=max(1,int(coast/2)+1);dt=coast/nt
            for _ in range(nt):
                def cr(state):return evaluate(state,j+1,mass,0,False)[0]
                k1=cr(y);k2=cr(y+k1*dt/2);k3=cr(y+k2*dt/2);k4=cr(y+k3*dt)
                y+=dt*(k1+2*k2+2*k3+k4)/6
                if save:rows.append(evaluate(y,j+1,mass,0,False)[1])
    if save:
        # Remove same-stage duplicate event-time samples; preserve discontinuities across stages.
        for stage in (1,2,3):
            ii=[k for k,row in enumerate(rows) if row['stage']==stage]
            ts=np.array([rows[k]['t_s'] for k in ii]);pp=np.array([rows[k]['pitch_deg'] for k in ii])
            unique=np.r_[True,np.diff(ts)>1e-8]
            rate=np.gradient(pp[unique],ts[unique]) if sum(unique)>1 else np.zeros(sum(unique))
            lookup=dict(zip(ts[unique],rate))
            for k in ii:rows[k]['pitch_rate_degps']=float(lookup[rows[k]['t_s']])
    return rows,events,y,mass


def target(c,out,quick=False):
    x0=seed(c);hist=[];R=RE+c['mission']['orbit_altitude_m'];V=np.sqrt(MU/R)
    def residual(x):
        try:
            rows,_,y,_=simulate(c,x,resolution=70,save=True);r=y[:3];v=y[3:6];er=r/np.linalg.norm(r)
            el=elements(y);q=max(z['q_Pa'] for z in rows);minh=min(z['h_m'] for z in rows[20:])
            result=np.r_[(np.linalg.norm(r)-R)/2000,np.dot(v,er)/20,
                (np.linalg.norm(v-np.dot(v,er)*er)-V)/20,
                (el['inclination_deg']-c['mission']['inclination_deg'])/.1,
                max(0,q/c['mission']['q_limit_Pa']-1)*5,max(0,-minh)/2000,
                (x[:6]-x0[:6])*.0005]
            hist.append(dict(evaluation=len(hist),scaled_residual_norm=np.linalg.norm(result[:6]),
                altitude_m=np.linalg.norm(r)-RE,radial_speed_mps=np.dot(v,er),speed_mps=np.linalg.norm(v),max_q_Pa=q))
            if len(hist)%30==0:
                print(f'Ascent targeting: evaluation {len(hist)}, residual={np.linalg.norm(result[:6]):.4g}',flush=True)
                json_write(out/'data/targeting_checkpoint.json',dict(guidance=x,input_hashes=c['input_hashes'],last=hist[-1]))
            return result
        except (RuntimeError,ValueError,FloatingPointError):return np.ones(12)*1e4
    bounds=([45,20,0,-10,-20,-15,0,.25,-2],[88,80,55,35,35,15,450,1,2])
    sol=least_squares(residual,x0,bounds=bounds,x_scale=[20,20,20,15,15,15,80,.2,.2],
        diff_step=2e-4,max_nfev=40 if quick else 70,ftol=2e-6,xtol=2e-6,gtol=2e-6)
    json_write(out/'data/coarse_guidance.json',dict(guidance=sol.x))
    # Correct discretization at finer resolution using orbital constraints only.
    refined=300
    active=np.array([4,5,7,8])
    def refine(vary):
        x=sol.x.copy();x[active]=vary
        _,_,y,_=simulate(c,x,resolution=refined,save=False);r=y[:3];v=y[3:6];er=r/np.linalg.norm(r)
        return np.r_[(np.linalg.norm(r)-R)/2000,np.dot(v,er)/20,(np.linalg.norm(v-np.dot(v,er)*er)-V)/20,(elements(y)['inclination_deg']-c['mission']['inclination_deg'])/.1]
    fine=least_squares(refine,sol.x[active],bounds=(np.array(bounds[0])[active],np.array(bounds[1])[active]),diff_step=1e-4,max_nfev=18,ftol=1e-8,xtol=1e-7,gtol=1e-7)
    finalx=sol.x.copy();finalx[active]=fine.x
    rows,events,y,mass=simulate(c,finalx,resolution=refined)
    if hist:write_csv(out/'data/ascent_targeting_iterations.csv',hist)
    write_csv(out/'data/ascent.csv',rows);write_csv(out/'data/events.csv',events)
    names=['S1_pitch_30_deg','S1_pitch_65_deg','S1_end_deg','S2_end_deg','S3_start_deg','S3_end_deg','coast23_s','S3_use_fraction','plane_trim_deg']
    write_csv(out/'data/optimized_guidance.csv',[dict(parameter=k,value=v) for k,v in zip(names,finalx)])
    el=elements(y);last=rows[-1]
    status=dict(**el,final_altitude_m=last['h_m'],final_speed_mps=last['speed_inertial_mps'],final_radial_speed_mps=last['vr_mps'],
        injection_time_s=y[6],max_q_Pa=max(r['q_Pa'] for r in rows),max_qalpha_Pa_rad=max(r['qalpha_Pa_rad'] for r in rows),
        max_proper_g=max(r['proper_g'] for r in rows),optimizer_success=bool(fine.success),
        integrated_thrust_dv_mps=y[7],drag_loss_mps=y[8],gravity_projection_loss_mps=y[9],
        remaining_stage3_propellant_kg=c['stages'][2]['propellant_kg']-c['stages'][2]['usable_propellant_kg']*finalx[7],
        crew_separation_mass_kg=c['mission']['payload_kg'])
    status['orbit_target_pass']=bool(abs(status['final_altitude_m']-c['mission']['orbit_altitude_m'])<c['mission']['orbit_alt_tolerance_m']
        and abs(status['final_speed_mps']-V)<c['mission']['orbit_speed_tolerance_mps'] and abs(status['final_radial_speed_mps'])<2
        and abs(status['inclination_deg']-c['mission']['inclination_deg'])<c['mission']['orbit_incl_tolerance_deg'])
    # Independent doubled-resolution integration check of the same guidance.
    _,_,yy,_=simulate(c,finalx,resolution=600,save=False)
    status['integration_refinement_position_m']=np.linalg.norm(yy[:3]-y[:3]);status['integration_refinement_velocity_mps']=np.linalg.norm(yy[3:6]-y[3:6])
    return rows,events,y,status,finalx
