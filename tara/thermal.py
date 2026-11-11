"""Ascent heating correlations, implicit 1-D TPS, cryogenic heat leak and vents."""
import numpy as np
from scipy.linalg import solve_banded
from .io import write_csv

SIGMA=5.670374419e-8


def conduction(c,times,flux,thickness,keep=False):
    th=c['thermal'];mat=c['materials'][1];N=31;dx=thickness/(N-1)
    cap=np.full(N,mat['density_kgm3']*mat['cp_JkgK']*dx);cap[[0,-1]]*=.5
    cap[-1]+=.003*c['materials'][0]['density_kgm3']*c['materials'][0]['cp_JkgK']
    conduct=mat['k_WmK']/dx;T=np.full(N,th['initial_K']);history=[];maxback=T[-1]
    for k in range(1,len(times)):
        span=times[k]-times[k-1]
        if span<=0:continue
        steps=max(1,int(np.ceil(span/2)));dt=span/steps
        for j in range(steps):
            f=(j+1)/steps;heat=(1-f)*flux[k-1]+f*flux[k]
            diag=cap/dt+2*conduct;diag[[0,-1]]-=conduct
            rhs=cap/dt*T
            radlin=4*th['emissivity']*SIGMA*T[0]**3
            diag[0]+=radlin
            rhs[0]+=heat+th['emissivity']*SIGMA*(3*T[0]**4+th['ambient_K']**4)
            ab=np.zeros((3,N));ab[0,1:]=-conduct;ab[2,:-1]=-conduct;ab[1]=diag
            T=solve_banded((1,1),ab,rhs)
            maxback=max(maxback,T[-1])
        if keep:history.append(dict(t_s=times[k],surface_K=T[0],bondline_K=T[-1],midwall_K=T[N//2],heat_flux_Wm2=flux[k]))
    return maxback,history,T


def orifice_flux(p0,p1,T,R=287.05287,gamma=1.4,Cd=.65):
    """One-way ideal-gas mass flux per area; choked or isentropic subcritical."""
    if p0<=p1 or p0<=0:return 0.
    ratio=max(0,p1/p0);critical=(2/(gamma+1))**(gamma/(gamma-1))
    if ratio<=critical:return Cd*p0/np.sqrt(R*T)*np.sqrt(gamma)*(2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
    return Cd*p0/np.sqrt(R*T)*np.sqrt(2*gamma/(gamma-1)*(ratio**(2/gamma)-ratio**((gamma+1)/gamma)))


def run(c,ascent,out):
    # Collapse staging discontinuities for the thermal time coordinate.
    rows={float(r['t_s']):r for r in ascent};rows=[rows[t] for t in sorted(rows)]
    times=np.array([r['t_s'] for r in rows]);th=c['thermal'];nose=c['geometry']['nose_radius_m']
    heat=np.array([th['sutton_graves_k']*np.sqrt(max(r['rho_kgm3'],0)/nose)*r['v_mps']**3 for r in rows])
    heat[np.array([r['h_m']>150000 for r in rows])]=0 # continuum heating correlation not extended into rarefied orbit
    history=[]
    for r,h in zip(rows,heat):
        stage=c['stages'][int(r['stage'])-1];power=stage['mdot_full_kgps']*r['throttle']*th['combustion_energy_Jkg']
        base=th['base_heat_fraction']*power/(np.pi*c['diameter_m']**2/4)
        history.append(dict(t_s=r['t_s'],stagnation_flux_Wm2=h,base_heat_proxy_Wm2=base,
            radiative_equilibrium_K=(h/(th['emissivity']*SIGMA)+th['ambient_K']**4)**.25,
            cumulative_stagnation_heat_Jm2=float(np.trapezoid(heat[:len(history)+1],times[:len(history)+1])),
            model='Sutton-Graves ascent correlation; plume heating is an assumed power fraction'))
    lo,hi=th['tps_min_m'],th['tps_max_m'];iterations=[]
    minimum,_hist,_=conduction(c,times,heat,lo)
    maximum,_,_=conduction(c,times,heat,hi)
    if minimum<=th['bondline_limit_K']:chosen=lo
    elif maximum>th['bondline_limit_K']:chosen=hi
    else:
        for k in range(18):
            mid=(lo+hi)/2;back,_,_=conduction(c,times,heat,mid)
            iterations.append(dict(iteration=k,thickness_m=mid,max_bondline_K=back))
            if back>th['bondline_limit_K']:lo=mid
            else:hi=mid
        chosen=hi
    back,tpshistory,final=conduction(c,times,heat,chosen,True)
    # Spatial ascent capsule heat envelope: modified-Newtonian-type exposure proxy.
    theta=np.linspace(0,np.pi/2,121)
    spatial=[dict(surface_angle_deg=np.rad2deg(x),arc_m=nose*x,heat_flux_Wm2=max(heat)*np.cos(x)**1.5,
        local_TPS_thickness_m=chosen,model='prescribed exposure distribution; no viscous CFD heat-transfer solution') for x in theta]
    cryo=[]
    for t in c['tanks']:
        if t['fluid']=='RP1':continue
        fluid=next(f for f in c['fluids'] if f['fluid']==t['fluid'])
        leak=th['cryo_foam_k_WmK']*t['area_m2']*(th['ambient_K']-t['temperature_K'])/t['insulation_m']
        hold=th['pad_hold_s']+times[-1];boilrate=leak/fluid['latent_Jkg'];boilmass=boilrate*hold
        R=th[t['fluid']+'_vent_R_JkgK']
        # Orifice sized conservatively against near-vacuum downstream; use warm upstream gas assumption.
        area=1.5*boilrate/max(orifice_flux(t['pressure_abs_Pa'],101325,max(t['temperature_K'],100),R,Cd=th['vent_Cd']),1e-20)
        cryo.append(dict(tank=t['name'],heat_leak_W=leak,pad_plus_flight_boiloff_kg=boilmass,
            boiloff_fraction=boilmass/t['loaded_mass_kg'],available_reserve_kg=t['loaded_mass_kg']*.01,
            reserve_covers_boiloff=bool(boilmass<t['loaded_mass_kg']*.01),vent_area_m2=area,
            equivalent_vent_diameter_m=np.sqrt(4*area/np.pi),note='Isothermal latent-heat upper screen; no two-phase flow or pressurization control'))
    # Bay depressurization requirement from ambient dP/dt, sized at allowable overpressure.
    p=np.array([r['p_Pa'] for r in rows]);dpdt=np.gradient(p,times)
    bays=[];venthistory=[]
    for s in c['stages']:
        volume=np.pi*c['diameter_m']**2/4*s['interstage_m'];required=volume/(287.05287*th['bay_temperature_K'])*np.maximum(0,-dpdt)
        flux=np.array([orifice_flux(pa+th['vent_diff_pressure_Pa'],pa,th['bay_temperature_K'],Cd=th['vent_Cd']) for pa in p])
        area=float(np.max(required/np.maximum(flux,1e-20))*1.5)
        bays.append(dict(stage=s['stage'],bay_volume_m3=volume,vent_area_m2=area,equivalent_diameter_m=np.sqrt(4*area/np.pi),
            allowable_diff_pressure_Pa=th['vent_diff_pressure_Pa'],note='Quasi-steady ideal-gas vent sizing; detailed transient valves not modeled'))
        venthistory.extend(dict(stage=s['stage'],t_s=t,ambient_pressure_Pa=pa,required_mass_flow_kgps=mf,capacity_mass_flow_kgps=ff*area) for t,pa,mf,ff in zip(times,p,required,flux))
    for name,data in [('thermal_history',history),('TPS_history',tpshistory),('capsule_heat_distribution',spatial),
                      ('cryo_heat_leak_and_vents',cryo),('interstage_vents',bays),('vent_sizing_history',venthistory)]:write_csv(out/f'data/{name}.csv',data)
    if not iterations:iterations=[dict(iteration=0,thickness_m=chosen,max_bondline_K=back)]
    write_csv(out/'data/TPS_sizing_iterations.csv',iterations)
    status=dict(peak_stagnation_flux_Wm2=max(heat),integrated_heat_Jm2=np.trapezoid(heat,times),TPS_thickness_m=chosen,
        max_bondline_K=back,bondline_limit_K=th['bondline_limit_K'],TPS_screen_pass=bool(back<=th['bondline_limit_K']),
        TPS_nose_area_mass_kg=2*np.pi*c['geometry']['capsule_base_radius_m']**2*chosen*c['materials'][1]['density_kgm3'],
        cryo_reserve_pass=all(r['reserve_covers_boiloff'] for r in cryo),
        scope='Ascent only; insulation and plume screens; no re-entry thermal qualification; TPS mass audited against payload allowance')
    return history,tpshistory,spatial,cryo,bays,status
