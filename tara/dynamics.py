"""Beam modes, assumed acoustic/buffet response, POGO, slosh and separation."""
from math import comb
import numpy as np
from scipy.linalg import eigh
from .beam import BeamSection,modes
from .pogo import assess
from .io import G0,write_csv,json_write
from .structures import mass_distribution


def run(c,ascent,events,out):
    mat=c['materials'][0];st=c['structures'];d=c['dynamics'];r=c['diameter_m']/2
    area=2*np.pi*r*st['effective_skin_m'];I=np.pi*r**3*st['effective_skin_m'];E=mat['E_Pa'];G=E/(2*(1+mat['nu']))
    sections=[]
    for s in c['stages']:
        L=s['length_m'];dry=s['dry_kg']/L;wet=(s['dry_kg']+s['propellant_kg'])/L
        sections.append(BeamSection(L,E*area,E*I,.5*G*area,dry,wet,dry*r*r/2,wet*r*r/2,8))
    L=c['tip_z_m']-c['stages'][2]['forward_z_m'];mass=c['mission']['payload_kg']/L
    sections.append(BeamSection(L,E*area,E*I,.5*G*area,mass,mass,mass*r*r/2,mass*r*r/2,4))
    wet=modes(sections,True,count=8);dry=modes(sections,False,count=8)
    # Independently clamp all three base DOFs for the pad case.
    vals,_=eigh(wet['K'][3:,3:],wet['M'][3:,3:],subset_by_index=(0,5))
    groundfreq=np.sqrt(np.maximum(vals,0))/(2*np.pi)
    modal=[]
    for state,result in [('wet',wet),('dry',dry)]:
        for k,f in enumerate(result['frequency_Hz']):
            modal.append(dict(state=state,mode=k+1,frequency_Hz=f,axial_kinetic_fraction=result['axial_kinetic_fraction'][k],
                eigen_residual=result['eigen_residual'][k],boundary='free-free',note='homogenized beam; liquids rigidly attached'))
    modal.extend(dict(state='ground_wet',mode=k+1,frequency_Hz=f,boundary='cantilever') for k,f in enumerate(groundfreq))
    shapes=[]
    for k in range(6):
        pp=wet['modes'][:,1,k];pp=pp/max(abs(pp))
        shapes.extend(dict(mode=k+1,z_m=z,lateral_normalized=val,axial_normalized=axx/max(1e-30,max(abs(wet['modes'][:,0,k])))) for z,val,axx in zip(wet['z_m'],pp,wet['modes'][:,0,k]))
    f=np.geomspace(.1,2500,1500);w=2*np.pi*f
    prms=20e-6*10**(d['OASPL_dB']/20);shape=np.exp(-.5*(np.log(f/100)/.8)**2)/f;psd=prms**2*shape/np.trapezoid(shape,f)
    bending=[x for x,frac in zip(wet['frequency_Hz'],wet['axial_kinetic_fraction']) if frac<.5]
    fn=bending[0] if bending else wet['frequency_Hz'][0];wn=2*np.pi*fn;zeta=d['modal_damping'];modalmass=.25*ascent[0]['mass_kg'];forcearea=np.pi*r*r
    H=(-w*w)/(wn*wn-w*w+2j*zeta*wn*w)*forcearea/modalmass
    q=max(row['q_Pa'] for row in ascent);buffet_rms=d['buffet_rms_q_fraction']*q
    bfshape=np.exp(-.5*(np.log(f/15)/.65)**2)/f;bfpsd=buffet_rms**2*bfshape/np.trapezoid(bfshape,f)
    spectral=[dict(frequency_Hz=ff,pressure_PSD_Pa2_Hz=pp,buffet_PSD_Pa2_Hz=bb,acceleration_PSD_m2s4_Hz=hh*pp,
        buffet_acceleration_PSD_m2s4_Hz=hh*bb) for ff,pp,bb,hh in zip(f,psd,bfpsd,abs(H)**2)]
    ws=2*np.pi*d['pogo_struct_Hz'];wf=2*np.pi*d['pogo_feed_Hz'];gain=d['pogo_loop_gain']
    loop=gain*ws*ws*wf*wf/((ws*ws-w*w+2j*d['pogo_struct_damping']*ws*w)*(wf*wf-w*w+2j*d['pogo_feed_damping']*wf*w))*np.exp(-1j*w*d['pogo_delay_s'])
    pogostat=assess([ws*ws],[1,2*d['pogo_struct_damping']*ws,ws*ws],
        [gain*wf*wf],[1,2*d['pogo_feed_damping']*wf,wf*wf],d['pogo_delay_s'])
    pg=[dict(frequency_Hz=ff,loop_real=ll.real,loop_imag=ll.imag,loop_gain_dB=20*np.log10(max(abs(ll),1e-30)),
        loop_phase_deg=phase,closed_loop_magnitude=1/abs(1+ll)) for ff,ll,phase in zip(f,loop,np.rad2deg(np.unwrap(np.angle(loop))))]
    # Method-of-steps transient of the SAME assumed scalar delayed loop. The
    # 0.1 mm initial displacement is a test disturbance, not measured flight data.
    dt=.001;ptime=np.arange(0,15+dt/2,dt);Y=np.zeros((len(ptime),4));Y[0,0]=1e-4
    stop=len(ptime)-1
    for k in range(len(ptime)-1):
        def dd(tt,yy):
            delayed=float(np.interp(tt-d['pogo_delay_s'],ptime[:k+1],Y[:k+1,2],left=0))
            q,qd,feed,fd=yy
            return np.array([qd,-2*d['pogo_struct_damping']*ws*qd-ws*ws*q-gain*ws*ws*delayed,
                fd,-2*d['pogo_feed_damping']*wf*fd-wf*wf*feed+wf*wf*q])
        a=dd(ptime[k],Y[k]);b=dd(ptime[k]+dt/2,Y[k]+dt*a/2);cc=dd(ptime[k]+dt/2,Y[k]+dt*b/2);ee=dd(ptime[k]+dt,Y[k]+dt*cc)
        Y[k+1]=Y[k]+dt*(a+2*b+2*cc+ee)/6
        if max(abs(Y[k+1,0]),abs(Y[k+1,2]))>=.01:
            stop=k+1;break # no extrapolation of this linear model to huge physical displacements
    ptime=ptime[:stop+1];Y=Y[:stop+1]
    transient=[dict(t_s=t,structure_displacement_m=q,feed_equivalent_displacement_m=feed,
        scope='assumed delayed loop; 0.1 mm disturbance; stops at 10 mm linear-range screen; not measured flight') for t,q,feed in zip(ptime[::10],Y[::10,0],Y[::10,2])]
    slosh=[]
    for tank in c['tanks']:
        rr=tank['inner_radius_m'];h=tank['z1_m']-tank['z0_m'];wave=1.841/rr
        for fill in [.1,.25,.5,.75,1.]:
            freq=np.sqrt(G0*wave*np.tanh(wave*h*fill))/(2*np.pi)
            slosh.append(dict(tank=tank['name'],fill_fraction=fill,frequency_Hz=freq,
                unbaffled_damping_assumption=.005,baffled_damping_assumption=.03,baffle_count=st['baffle_count'],
                note='cylindrical first lateral slosh mode at 1g; dome and acceleration history not resolved'))
    z=np.linspace(0,c['tip_z_m'],301);U=d['ground_wind_mps']*(np.maximum(z,2)/10)**.14
    line=.5*1.225*U*U*1.2*c['diameter_m']*d['gust_factor'];dz=z[1]-z[0]
    shear=np.cumsum(line[::-1])[::-1]*dz;moment=np.array([np.trapezoid(line[k:]*(z[k:]-zz),z[k:]) for k,zz in enumerate(z)])
    slope=np.cumsum(moment/(E*I))*dz;deflect=np.cumsum(slope)*dz
    gw=[dict(z_m=zz,wind_mps=u,line_load_Nm=p,shear_N=v,bending_Nm=b,deflection_m=y) for zz,u,p,v,b,y in zip(z,U,line,shear,moment,deflect)]
    separation=[]
    for event in events:
        if 'separation' not in event['event']:continue
        ma=event['mass_before_kg']-event['mass_after_kg'];mb=event['mass_after_kg'];force=d['separation_push_N'];burn=d['separation_push_duration_s']
        for t in np.linspace(0,15,200):
            duration=min(t,burn);relacc=force*(1/ma+1/mb)
            distance=.5*relacc*duration**2+relacc*burn*max(0,t-burn)
            tip=np.deg2rad(d['separation_tipoff_degps'])*t
            separation.append(dict(event=event['event'],time_after_sep_s=t,axial_clearance_m=distance,
                tipoff_deg=np.rad2deg(tip),lateral_tip_sweep_m=c['diameter_m']/2*np.sin(tip),
                required_axial_clearance_m=d['separation_clearance_m'],note='two-body collinear push; no plume/contact/6DOF'))
    redundancy=[]
    for r0 in c['human_rating']:
        n=int(r0['channels']);k=int(r0['required_healthy']);p=r0['unit_failure_probability'];beta=r0['common_cause_beta']
        independent=sum(comb(n,j)*(1-p)**j*p**(n-j) for j in range(k))
        fail=beta*p+(1-beta)*independent
        redundancy.append(dict(**r0,independent_function_failure_probability=independent,
            illustrative_common_cause_function_failure_probability=fail,qualification_status='NOT_A_CREW_LOSS_PROBABILITY'))
    for name,data in [('modal_frequencies',modal),('mode_shapes',shapes),('acoustic_buffet_response',spectral),
                      ('pogo_frequency_response',pg),('slosh_modes',slosh),('ground_wind_loads',gw),('separation_clearance',separation),
                      ('redundancy_assessment',redundancy),('pogo_transient',transient)]:write_csv(out/f'data/{name}.csv',data)
    json_write(out/'data/pogo_assessment.json',pogostat)
    return dict(wet=wet,dry=dry,modal=modal,shapes=shapes,spectral=spectral,pogo=pg,pogostat=pogostat,
        slosh=slosh,ground=gw,separation=separation,redundancy=redundancy,pogo_transient=transient,
        status=dict(first_wet_bending_Hz=fn,first_ground_mode_Hz=groundfreq[0],
            shedding_Hz=.2*d['ground_wind_mps']/c['diameter_m'],ground_tip_deflection_m=deflect[-1],
            ground_base_moment_Nm=moment[0],acoustic_RMS_g=np.sqrt(np.trapezoid(abs(H)**2*psd,f))/G0,
            pogo_status=pogostat['status'],assumed_acoustic_OASPL_dB=d['OASPL_dB'],
            pogo_transient_end_s=float(ptime[-1]),pogo_transient_linear_limit_reached=bool(stop<len(np.arange(0,15+dt/2,dt))-1),
            note='Linear response to assumed spectra and loop; no validated acoustic source or coupled aeroservoelastic solution'))
