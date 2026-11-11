"""Publication-style PNG figures; numeric sources remain in results/data."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Ellipse,Rectangle,Polygon
from matplotlib import animation
from .io import RE,table
from .geometry import radius

C=['#126782','#ee8c31','#639f72','#8d5c9e','#c65353','#63778b']


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,
        'axes.titleweight':'bold','axes.labelcolor':'#243746','text.color':'#243746','axes.edgecolor':'#a8b8c4',
        'axes.prop_cycle':plt.cycler(color=C),'grid.alpha':.22,'axes.grid':True,'figure.facecolor':'white',
        'savefig.facecolor':'white','legend.frameon':False,'lines.linewidth':1.8})


def values(rows,key):return np.array([r[key] for r in rows],float)


def save(fig,out,name,title,note='Conceptual design | engineering models and assumed inputs | see model status'):
    fig.suptitle('TARA_LV  /  '+title,fontsize=17,fontweight='bold',x=.03,ha='left')
    fig.text(.03,.015,note,fontsize=8,color='#63778b')
    fig.tight_layout(rect=(.015,.04,.99,.94))
    path=out/'figures'/f'{name}.png';path.parent.mkdir(exist_ok=True)
    fig.savefig(path,dpi=190);plt.close(fig)


def line(ax,x,y,title,xlabel,ylabel,label=None):
    ax.plot(x,y,label=label);ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
    if label:ax.legend(fontsize=8)


def make(c,out,ascent,events,stages,trade,layout,envelope,margins,closure,dyn,thermal,rv,flows):
    style();t=values(ascent,'t_s');hm=values(ascent,'h_m')/1000;time='Time from liftoff [s]';xpos='Axial station z from base [m]'
    _mission_overview(c,out,ascent,rv,dyn)
    fig,ax=plt.subplots(3,3,figsize=(15,11))
    spec=[('h_m',.001,'Altitude','km'),('speed_inertial_mps',.001,'Inertial speed','km/s'),('mass_kg',.001,'Attached mass','t'),
        ('pitch_deg',1,'Thrust elevation','deg'),('pitch_rate_degps',1,'Thrust elevation rate','deg/s'),('thrust_N',1e-6,'Net thrust','MN'),
        ('q_Pa',.001,'Dynamic pressure','kPa'),('proper_g',1,'Proper acceleration','g'),('alpha_deg',1,'Angle of attack','deg')]
    for a,(key,f,title,unit) in zip(ax.flat,spec):
        line(a,t,values(ascent,key)*f,title,time,unit)
        for e in events:
            if 'separation' in e['event']:a.axvline(e['t_s'],color='#999',ls=':',lw=.8)
    ax[2,0].axhline(c['mission']['q_limit_Pa']/1000,color=C[4],ls='--');ax[2,1].axhline(c['mission']['proper_g_limit'],color=C[4],ls='--')
    save(fig,out,'05_ascent_time_histories','Ascent flight history')
    fig,ax=plt.subplots(2,3,figsize=(15,8))
    for a,key,title,unit in zip(ax.flat,['throttle','qalpha_Pa_rad','integrated_thrust_dv_mps','drag_loss_mps','gravity_loss_projection_mps','east_wind_mps'],
         ['Thrust shaping / q relief','q × incidence','Integrated thrust / mass','Accumulated drag loss','Gravity projection loss','Assumed ascent wind'],
         ['fraction','Pa rad','m/s','m/s','m/s','m/s']):line(a,t,values(ascent,key),title,time,unit)
    save(fig,out,'06_propulsion_and_velocity_ledger','Propulsion, loads and velocity accounting')
    fig,ax=plt.subplots(1,3,figsize=(15,5.5));ss=np.arange(3)
    ax[0].bar(ss,[r['propellant_kg']/1000 for r in stages],label='Loaded propellant',color=C[0]);ax[0].bar(ss,[r['dry_kg']/1000 for r in stages],bottom=[r['propellant_kg']/1000 for r in stages],label='Dry budget',color=C[1]);ax[0].set(xticks=ss,xticklabels=['S1','S2','S3'],ylabel='Mass [t]',title='Stage mass ledger');ax[0].legend()
    ax[1].bar(ss,[r['dv_ideal_mps']/1000 for r in stages],color=C[:3]);ax[1].set(xticks=ss,xticklabels=['S1','S2','S3'],ylabel='Ideal velocity [km/s]',title='Sizing allocation at selected Isp')
    ax[2].bar(range(3),[r['mass_kg']/1000 for r in trade],color=[C[2] if r['feasible'] else C[4] for r in trade]);ax[2].set(xticks=range(3),xticklabels=[r['case'] for r in trade],ylabel='Liftoff mass [t]',title='Structural-factor sensitivity');ax[2].tick_params(axis='x',rotation=15)
    save(fig,out,'02_mass_and_staging','Sizing allocation and structural-factor trade','Red bar: constrained sizing did not meet all requirements; optimistic factors need substantiation')
    fig,ax=plt.subplots(2,1,figsize=(18,5));zz=np.linspace(0,c['tip_z_m'],1500);rr=radius(c['tip_z_m']-zz,c)
    for a in ax:
        a.plot(zz,rr,color='#243746');a.plot(zz,-rr,color='#243746');a.set(xlabel=xpos,ylabel='Radius [m]',ylim=(-3.2,4),xlim=(-1,c['tip_z_m']+1));a.set_aspect('equal',adjustable='box')
    for s in c['stages']:
        ax[0].axvspan(s['aft_z_m'],s['forward_z_m'],alpha=.12,color=C[int(s['stage'])-1]);ax[0].text((s['aft_z_m']+s['forward_z_m'])/2,2.7,f'S{s["stage"]}: {s["length_m"]:.2f} m',ha='center')
    for tank in c['tanks']:
        for a in [ax[1]]:
            lo=tank['z0_m'];hi=tank['z1_m'];r=tank['inner_radius_m'];hd=tank['head_depth_m'];color=C[1] if tank['fluid']=='RP1' else C[0] if tank['fluid']=='LOX' else C[2]
            a.add_patch(Ellipse((lo+hd,0),2*hd,2*r,fc=color,alpha=.4));a.add_patch(Ellipse((hi-hd,0),2*hd,2*r,fc=color,alpha=.4))
            a.add_patch(Rectangle((lo+hd,-r),max(0,hi-lo-2*hd),2*r,fc=color,alpha=.4))
            a.text((lo+hi)/2,0,f'{tank["name"]}\n{hi-lo:.2f} m',ha='center',va='center',fontsize=8)
    ax[0].set_title(f'External envelope: {c["tip_z_m"]:.2f} m × {c["diameter_m"]:.2f} m');ax[1].set_title('Tank meridian: computed from loaded mass, density, ullage, wall and insulation')
    save(fig,out,'03_dimensioned_vehicle','Vehicle meridian and tank dimensions','Exact conceptual coordinates: data/components.csv, data/tanks.csv and data/vehicle_profile.dxf | no production CAD claim')
    fig,ax=plt.subplots(1,3,figsize=(14,5.5))
    for sid,a in enumerate(ax,1):
        a.add_patch(Circle((0,0),c['diameter_m']/2,fill=False,color='#243746',lw=2))
        for r in layout:
            if r['stage']!=sid:continue
            a.add_patch(Circle((r['x_m'],r['y_m']),r['nozzle_radius_m'],fc=C[sid-1],alpha=.45));a.add_patch(Circle((r['x_m'],r['y_m']),r['gimbal_envelope_radius_m'],fill=False,color=C[4],ls='--',lw=1));a.text(r['x_m'],r['y_m'],str(r['engine']),ha='center',va='center')
        a.set(xlim=(-2.5,2.5),ylim=(-2.5,2.5),xlabel='x [m]',ylabel='y [m]',title=f'Stage {sid}: exit and gimbal envelope');a.set_aspect('equal')
    save(fig,out,'04_engine_installation','Engine installation feasibility','Dashed: assumed swept gimbal envelope | overlaps are unresolved packaging constraints, not approved placements')
    db=rv['history'];rt=values(db,'t_s')/3600
    fig=plt.figure(figsize=(14,7));ax=fig.add_subplot(121,projection='3d');b=fig.add_subplot(122)
    u=np.linspace(0,2*np.pi,45);v=np.linspace(0,np.pi,24);ax.plot_surface(RE/1e3*np.outer(np.cos(u),np.sin(v)),RE/1e3*np.outer(np.sin(u),np.sin(v)),RE/1e3*np.outer(np.ones_like(u),np.cos(v)),color='#aecbd8',alpha=.4,linewidth=0)
    for phase,color in [('two_period_coast',C[0]),('J2_phasing',C[1]),('prescribed_docking_approach',C[2])]:
        rs=[r for r in db if r['phase']==phase];ax.plot(values(rs,'rx_m')/1000,values(rs,'ry_m')/1000,values(rs,'rz_m')/1000,color=color,label=phase.replace('_',' '))
    ax.set(xlabel='ECI x [km]',ylabel='ECI y [km]',zlabel='ECI z [km]');ax.set_box_aspect([1,1,1]);ax.legend(fontsize=8,loc='upper left')
    line(b,rt,values(db,'altitude_m')/1000,'Actual-state orbital mission','Time from liftoff [h]','Altitude [km]');b.axhline(c['mission']['orbit_altitude_m']/1000,ls='--',color=C[3],label='BAS nominal');b.legend()
    save(fig,out,'07_orbit_and_phasing','J2 propagation and phasing orbit')
    fig,ax=plt.subplots(2,1,figsize=(13,8))
    for rows,label,color in [(ascent,'Ascent',C[0]),(db,'Orbital crew',C[1])]:
        lon=values(rows,'lon_deg');lat=values(rows,'lat_deg');lon=lon.copy();lon[np.r_[False,abs(np.diff(lon))>180]]=np.nan
        ax[0].plot(lon,lat,label=label,color=color,lw=1.2)
    ax[0].set(xlim=(-180,180),ylim=(-70,70),xlabel='Longitude [deg]',ylabel='Geocentric latitude [deg]',title='Earth-fixed track at common UTC epoch');ax[0].legend()
    line(ax[1],t,hm,'Ascent altitude along track',time,'Altitude [km]')
    save(fig,out,'08_ground_tracks','Ground tracks','Spherical Earth + GMST | no overflight clearance, debris dispersion or descent simulation')
    fig,ax=plt.subplots(2,2,figsize=(13,8))
    ax[0,0].semilogy(rt,np.maximum(values(db,'relative_distance_m'),.01)/1000);ax[0,0].set(xlabel='Time [h]',ylabel='Separation [km]',title='Crew–BAS distance')
    line(ax[0,1],rt,values(db,'relative_speed_mps'),'Relative speed in rotating RTN','Time [h]','m/s')
    for key,label in [('radial_m','Radial'),('alongtrack_m','Along track'),('crosstrack_m','Cross track')]:ax[1,0].plot(rt,values(db,key)/1000,label=label)
    ax[1,0].set(xlabel='Time [h]',ylabel='RTN distance [km]',title='Relative geometry');ax[1,0].legend()
    app=rv['approach'];tt=(values(app,'t_s')-app[0]['t_s'])/60
    line(ax[1,1],tt,values(app,'relative_distance_m'),'Prescribed final approach','Approach time [min]','Distance [m]')
    save(fig,out,'09_rendezvous_history','Phasing, arrival and prescribed docking','Coast and phasing are numerical J2 solutions; final approach is a commanded path, not autonomous docking validation')
    fig,ax=plt.subplots(1,2,figsize=(12,5));aa=rv['approach']
    line(ax[0],values(aa,'alongtrack_m'),values(aa,'radial_m'),'Final approach in RTN','Along track [m]','Radial [m]')
    ax[0].set_ylim(-10,10) # keep roundoff at 1e-10 m from appearing as physical radial oscillation
    line(ax[1],tt,values(aa,'commanded_accel_mps2')*c['mission']['payload_kg'],'Inverse-dynamics thrust demand','Time [min]','N');ax[1].axhline(c['mission']['approach_thrust_available_N'],ls='--',color=C[4],label='Assumed available');ax[1].legend()
    save(fig,out,'10_docking_feasibility','Terminal approach acceleration screen')
    aero=rv['aerodb'];polars=rv['polars'];fig,ax=plt.subplots(2,3,figsize=(15,8))
    for a,key,title,unit in zip(ax.flat,['Ca','Cd','Cn','Cnalpha_per_rad','Xcp_from_nose_m','Cl_roll'],['Axial coefficient','Drag coefficient','Normal coefficient at 3 deg','Normal-force slope','Center of pressure','Rolling coefficient (symmetric body)'],['Ca','Cd','Cn','1/rad','m from nose','Cl']):
        for sid in [1,2,3]:
            rows=[r for r in aero if r['stage']==sid];a.plot(values(rows,'Mach'),values(rows,key),label=f'S{sid}')
        a.set(xlabel='Mach number',ylabel=unit,title=title);a.legend(fontsize=8)
        a.axvspan(.8,1.3,alpha=.08,color=C[1])
    save(fig,out,'11_aerodynamic_mach_sweep','Aerodynamic database across the ascent Mach range','Engineering correlations at stated reference Reynolds number | beyond perfect-gas CFD: screening only | zero roll is a symmetry assumption')
    fig,ax=plt.subplots(1,2,figsize=(12,5.5))
    for M in [.3,.98,2,5,15]:
        unique=sorted(set(r['Mach'] for r in polars));nearest=min(unique,key=lambda q:abs(q-M));rows=[r for r in polars if r['stage']==1 and r['Mach']==nearest]
        line(ax[0],values(rows,'alpha_deg'),values(rows,'Cn'),'Normal-force polar','Angle of attack [deg]','Cn',f'M={nearest:.2f}')
    for sid in [1,2,3]:
        rr=[r for r in aero if r['stage']==sid];ax[1].plot(values(rr,'Mach'),values(rr,'Cd_friction'),label=f'S{sid} friction');
    ax[1].set(xlabel='Mach',ylabel='Drag coefficient',title='Skin-friction contribution by exposed length');ax[1].legend()
    save(fig,out,'12_aero_polars_and_viscous_estimates','Incidence and viscous corrections','Laminar/turbulent flat-plate correlations are not resolved viscous flow fields')
    fig,ax=plt.subplots(2,2,figsize=(13,8))
    for key in ['Ca','Cd','Cn']:ax[0,0].plot(t,values(ascent,key),label=key)
    ax[0,0].set(xlabel=time,ylabel='Coefficient',title='Coefficients used in flight EOM');ax[0,0].legend()
    line(ax[0,1],t,values(ascent,'drag_N')/1000,'Aerodynamic axial load',time,'kN')
    line(ax[1,0],t,values(ascent,'normal_N')/1000,'Aerodynamic normal load',time,'kN')
    line(ax[1,1],hm,values(ascent,'east_wind_mps'),'Assumed wind profile','Altitude [km]','East wind [m/s]')
    save(fig,out,'13_aerodynamic_loads','Aerodynamic forces and flight environment')
    fig,ax=plt.subplots(1,3,figsize=(15,5))
    for a,key,title,unit in zip(ax,['axial_N','shear_N','bending_Nm'],['Axial envelope','Shear envelope','Bending envelope'],['MN','MN','MN m']):line(a,values(envelope,'z_m'),values(envelope,key)/1e6,title,xpos,unit)
    save(fig,out,'15_section_load_envelopes','Flight section-cut force and moment envelopes','All ascent samples | quasi-static trim, mass inertia, drag and normal force | data/loads_by_case.csv retains simultaneous cases')
    fig,ax=plt.subplots(1,2,figsize=(15,6));names=[r['component'] for r in margins];ym=[float(r.get('yield_margin') or np.nan) for r in margins];um=[float(r.get('ultimate_margin') or np.nan) for r in margins]
    yy=np.arange(len(names));ax[0].barh(yy-.17,ym,.34,label='Yield');ax[0].barh(yy+.17,um,.34,label='Ultimate');ax[0].axvline(0,color=C[4]);ax[0].set(yticks=yy,yticklabels=names,xlabel='Margin of safety',title='Metallic shell and joint screens');ax[0].legend()
    x=np.arange(3);ax[1].bar(x-.18,[r['dry_budget_kg']/1000 for r in closure],.36,label='Dry budget');ax[1].bar(x+.18,[r['estimated_hardware_kg']/1000 for r in closure],.36,label='Estimated hardware');ax[1].set(xticks=x,xticklabels=['S1','S2','S3'],ylabel='Mass [t]',title='Hardware estimate versus sizing budget');ax[1].legend()
    save(fig,out,'16_structural_margins_and_mass','Structural margins and hardware mass closure','Negative margins remain visible | pressure + envelope bending/axial screen | no qualification, fatigue, local joints or fracture-control closure')
    fig,ax=plt.subplots(2,3,figsize=(14,8))
    for k,a in enumerate(ax.flat):
        rows=[r for r in dyn['shapes'] if r['mode']==k+1];line(a,values(rows,'z_m'),values(rows,'lateral_normalized'),f'Wet free-free mode {k+1}: {dyn["wet"]["frequency_Hz"][k]:.2f} Hz',xpos,'Normalized lateral shape')
    save(fig,out,'17_vehicle_mode_shapes','Homogenized Timoshenko beam modes','Rigidly attached liquids; rigid joints; bending and axial modes may be interleaved; mode labels in data/modal_frequencies.csv')
    fig,ax=plt.subplots(1,3,figsize=(15,5.5));sp=dyn['spectral'];f=values(sp,'frequency_Hz')
    ax[0].loglog(f,values(sp,'pressure_PSD_Pa2_Hz'),label='Assumed acoustic');ax[0].loglog(f,values(sp,'buffet_PSD_Pa2_Hz'),label='Assumed buffet');ax[0].set(xlabel='Frequency [Hz]',ylabel='Pa²/Hz',title='Prescribed pressure spectra');ax[0].legend()
    ax[1].loglog(f,values(sp,'acceleration_PSD_m2s4_Hz'));ax[1].set(xlabel='Hz',ylabel='(m/s²)²/Hz',title='One-mode acoustic response')
    for tank in sorted(set(r['tank'] for r in dyn['slosh'])):
        rows=[r for r in dyn['slosh'] if r['tank']==tank];ax[2].plot(values(rows,'fill_fraction'),values(rows,'frequency_Hz'),label=tank)
    ax[2].set(xlabel='Fill fraction',ylabel='Hz',title='First lateral slosh mode at 1 g');ax[2].legend(fontsize=7)
    save(fig,out,'18_acoustics_buffet_and_slosh','Acoustics, transonic buffet and tank slosh','Input spectra and baffle damping are assumptions; no acoustic-source or turbulent buffet simulation')
    fig,ax=plt.subplots(1,3,figsize=(15,5.5));pg=dyn['pogo'];ff=values(pg,'frequency_Hz')
    ax[0].plot(values(pg,'loop_real'),values(pg,'loop_imag'));ax[0].plot(-1,0,'x',color=C[4]);ax[0].set(xlabel='Re(L)',ylabel='Im(L)',title='Delay-aware POGO loop');ax[0].set_xlim(-3,3);ax[0].set_ylim(-3,3)
    ax[1].semilogx(ff,values(pg,'loop_gain_dB'));ax[1].axhline(0,color=C[4],ls='--');ax[1].set(xlabel='Hz',ylabel='dB',title='Open-loop gain')
    ax[2].loglog(ff,values(pg,'closed_loop_magnitude'));ax[2].set(xlabel='Hz',ylabel='|1/(1+L)|',title=dyn['pogostat']['status'].replace('_',' '))
    save(fig,out,'19_pogo_feedback','Longitudinal propulsion–structure feedback','Scalar assumed feed/structure transfer functions; exact-delay Nyquist screen; no identified engine/feed-system model')
    pogo_transient(out,dyn['pogo_transient'])
    fig,ax=plt.subplots(1,3,figsize=(15,5));gw=dyn['ground'];z=values(gw,'z_m')
    for a,key,title,unit in zip(ax,['wind_mps','bending_Nm','deflection_m'],['Pad wind profile','Cantilever bending','Static ground-wind deflection'],['m/s','N m','m']):line(a,values(gw,key),z,title,unit,'Height [m]')
    save(fig,out,'14_ground_wind_loads','Ground wind load screen','Assumed 20 m/s reference wind, 1.5 gust factor and rigid-base beam; no measured turbulence or vortex lock-in simulation')
    th,tps,spatial,cryo,bays,_=thermal;fig,ax=plt.subplots(2,2,figsize=(13,8))
    line(ax[0,0],values(th,'t_s'),values(th,'stagnation_flux_Wm2')/1000,'Ascent stagnation heating',time,'kW/m²')
    line(ax[0,1],values(th,'t_s'),values(th,'base_heat_proxy_Wm2')/1000,'Assumed plume/base exposure',time,'kW/m²')
    for key in ['surface_K','bondline_K','midwall_K']:ax[1,0].plot(values(tps,'t_s'),values(tps,key),label=key.replace('_',' '))
    ax[1,0].axhline(c['thermal']['bondline_limit_K'],ls='--',color=C[4]);ax[1,0].set(xlabel=time,ylabel='K',title='Implicit conduction TPS sizing');ax[1,0].legend()
    line(ax[1,1],values(spatial,'surface_angle_deg'),values(spatial,'heat_flux_Wm2')/1000,'Capsule exposure envelope','Angle from stagnation [deg]','kW/m²')
    save(fig,out,'20_ascent_thermal_and_TPS','Ascent aeroheating and thermal protection','Sutton–Graves + assumed spatial exposure; 1-D transient conduction + radiation; no descent, ablation or plume CFD')
    fig,ax=plt.subplots(1,3,figsize=(15,5.5));xx=np.arange(len(cryo));names=[r['tank'] for r in cryo]
    for a,key,title,unit in zip(ax,['heat_leak_W','pad_plus_flight_boiloff_kg','equivalent_vent_diameter_m'],['Cryogenic conductive heat leak','Pad + ascent boiloff screen','Assumed cryogen vent sizing'],['W','kg','m']):
        a.bar(xx,values(cryo,key),color=C[0]);a.set(xticks=xx,xticklabels=names,ylabel=unit,title=title);a.tick_params(axis='x',rotation=30)
    save(fig,out,'21_cryo_insulation_and_vents','Cryogenic insulation and vent sizing','Latent heat, ideal-gas orifice and fixed temperature assumptions; no two-phase vent/pressurization model')
    fig,ax=plt.subplots(1,2,figsize=(12,5.5))
    for event in sorted(set(r['event'] for r in dyn['separation'])):
        sr=[r for r in dyn['separation'] if r['event']==event];ax[0].plot(values(sr,'time_after_sep_s'),values(sr,'axial_clearance_m'),label=event);ax[1].plot(values(sr,'time_after_sep_s'),values(sr,'lateral_tip_sweep_m'),label=event)
    ax[0].axhline(c['dynamics']['separation_clearance_m'],color=C[4],ls='--');ax[0].set(xlabel='Time after release [s]',ylabel='Axial gap [m]',title='Cold separation push');ax[1].set(xlabel='Time [s]',ylabel='Lateral sweep [m]',title='Assumed tip-off sweep');ax[0].legend();ax[1].legend()
    save(fig,out,'22_stage_separation','Stage separation clearance screen','Collinear two-body impulse and assumed angular rate; no hot staging or six-degree-of-freedom re-contact model')
    _flow_figures(c,out,*flows)
    _elements(c,out)
    _iterations(out)
    return sorted(p.name for p in (out/'figures').glob('*.png'))


def _mission_overview(c,out,ascent,rv,dyn):
    fig,ax=plt.subplots(2,1,figsize=(14,7.5),gridspec_kw={'height_ratios':[1,2]})
    phases=[('Powered ascent',0,ascent[-1]['t_s']/3600,C[0])]
    for name,label,color in [('two_period_coast','Two nominal periods',C[1]),('J2_phasing','J2-targeted phasing',C[2]),('prescribed_docking_approach','Prescribed approach',C[3])]:
        rows=[r for r in rv['history'] if r['phase']==name];lo=rows[0]['t_s']/3600;hi=rows[-1]['t_s']/3600;phases.append((label,lo,hi-lo,color))
    for j,(name,start,duration,color) in enumerate(phases):ax[0].barh(j,duration,left=start,color=color);ax[0].text(start+duration+.05,j,f'{duration:.2f} h',va='center',fontsize=9)
    ax[0].set(yticks=range(4),yticklabels=[p[0] for p in phases],xlabel='Time from liftoff [h]',title='Mission sequence; phase durations are computed')
    tt=values(rv['history'],'t_s')/3600;ax[1].semilogy(tt,np.maximum(values(rv['history'],'relative_distance_m'),.01)/1000);ax[1].set(xlabel='Time from liftoff [h]',ylabel='Crew–BAS range [km]',title='Rendezvous range; terminal contact prescribed')
    save(fig,out,'01_mission_overview','Mission sequence and relative geometry','12,500 kg spacecraft | 400 km / 51.6 deg target | liquid–liquid–cryogenic | no descent analysis')
    fig,ax=plt.subplots(1,2,figsize=(14,6));rr=dyn['redundancy'];xx=np.arange(len(rr))
    ax[0].bar(xx-.18,[r['independent_function_failure_probability'] for r in rr],.36,label='Independent channels');ax[0].bar(xx+.18,[r['illustrative_common_cause_function_failure_probability'] for r in rr],.36,label='With assumed common cause');ax[0].set_yscale('log');ax[0].set(xticks=xx,xticklabels=[r['function'].replace('_','\n') for r in rr],ylabel='Illustrative function failure probability',title='Sensitivity to common-cause assumptions');ax[0].legend(fontsize=8)
    for sid in [1,2]:
        rows=[r for r in ascent if r['stage']==sid and r['throttle']>0];count=c['stages'][sid-1]['engine_count'];tw=values(rows,'thrust_N')*(count-1)/count/(values(rows,'mass_kg')*9.80665)
        ax[1].plot(values(rows,'t_s'),tw,label=f'S{sid}: one engine unavailable')
    ax[1].axhline(1,color=C[4],ls='--');ax[1].set(xlabel='Time [s]',ylabel='Instantaneous thrust / weight',title='One-engine-unavailable thrust screen');ax[1].legend(fontsize=8)
    save(fig,out,'30_redundancy_and_engine_thrust','Human-flight provisions: reduced-order screens','Illustrative reliability, not crew risk | engine-out curve holds nominal mass/throttle and is not a failed-engine trajectory or control assessment')


def _iterations(out):
    specs=[('mass_geometry_coupling.csv','iteration','max_hardware_over_budget_kg','Mass / geometry coupling','Over-budget hardware [kg]'),
        ('ascent_targeting_iterations.csv','evaluation','scaled_residual_norm','Ascent targeting','Scaled residual norm'),
        ('rendezvous_targeting_iterations.csv','evaluation','miss_m','J2 transfer shooting','Endpoint miss [m]')]
    fig,ax=plt.subplots(1,3,figsize=(15,5.5))
    for a,(name,x,y,title,unit) in zip(ax,specs):
        rows=table(out/'data'/name);a.semilogy(values(rows,x),np.maximum(values(rows,y),1e-8));a.set(xlabel=x.capitalize(),ylabel=unit,title=title)
    save(fig,out,'29_design_iteration_histories','Actual design-loop convergence','Residual histories record calculations performed; optimizer convergence does not close model-fidelity or qualification gates')


def pogo_transient(out,tr):
    fig,ax=plt.subplots(1,2,figsize=(13,5));tt=values(tr,'t_s');qq=values(tr,'structure_displacement_m')
    ax[0].plot(tt,qq*1000,label='Structure');ax[0].plot(tt,values(tr,'feed_equivalent_displacement_m')*1000,alpha=.7,label='Feed equivalent');ax[0].legend();ax[0].set(xlabel='Time [s]',ylabel='Equivalent displacement [mm]',title='Assumed POGO loop: free transient')
    ax[1].semilogy(tt,np.maximum(abs(qq),1e-12));ax[1].set(xlabel='Time [s]',ylabel='Absolute displacement [m]',title='Growth / decay diagnostic')
    save(fig,out,'32_pogo_oscillation_history','POGO transient from the supplied delayed model','0.1 mm test disturbance | stops at 10 mm linear-range screen | unstable assumed loop requires identification and mitigation')


def _flow_figures(c,out,cases,solutions,study,cp):
    if solutions:
        chosen=[min(solutions,key=lambda x:abs(x[1]['Mach']-m)) for m in [.3,.98,2,5]]
        fig,axes=plt.subplots(len(chosen),3,figsize=(17,3*len(chosen)))
        for axs,(label,row,res) in zip(axes,chosen):
            mesh=res['mesh'];rho=res['U'][...,0];u=res['U'][...,1]/rho;v=res['U'][...,2]/rho
            for a,field,title in zip(axs,[u/row['v_mps'],v/row['v_mps'],(res['p']-row['p_Pa'])/row['q_Pa']],['u/U∞','v/U∞','Cp field']):
                im=a.pcolormesh(mesh.X,mesh.R,field,shading='flat',cmap='RdBu_r' if title!='u/U∞' else 'viridis',rasterized=True)
                a.plot(mesh.x,mesh.R[:,0],color='black',lw=1);a.set(xlabel='Distance from nose [m]',ylabel='Radius [m]',title=f'{label}: {title} | '+('converged' if res['summary']['converged'] else 'iteration-limited'))
                a.set_xlim(-3,min(c['tip_z_m'],28));plt.colorbar(im,ax=a,shrink=.8)
        save(fig,out,'23_computed_Euler_fields','Computed axial velocity, radial velocity and pressure','Actual finite-volume iterates | inviscid, zero incidence | no viscous separation or laminar/turbulent field claims')
        label,row,res=chosen[1];mesh=res['mesh'];fig,ax=plt.subplots(2,1,figsize=(15,7))
        for a,limit in zip(ax,[None,(-3,12)]):
            a.plot(mesh.X,mesh.R,color='#8ab3c2',lw=.4);a.plot(mesh.X.T,mesh.R.T,color='#8ab3c2',lw=.4);a.plot(mesh.x,mesh.R[:,0],color='#243746',lw=2);a.set(xlabel='Nose-based coordinate [m]',ylabel='Radius [m]')
            if limit:a.set_xlim(*limit);a.set_ylim(0,5)
        save(fig,out,'24_FVM_grid','Body-fitted finite-volume ring mesh','Axisymmetric ring-volume discretization; exact metric check retained per case')
        fig,ax=plt.subplots(1,2,figsize=(14,5.5))
        for label,row,res in solutions:
            h=res['history'];ax[0].semilogy(values(h,'iteration'),values(h,'residual'),label=label);ax[1].plot(values(h,'iteration'),values(h,'CD_pressure_only'),label=label)
        ax[0].axhline(c['aerodynamics']['cfd_residual_tolerance'],color=C[4],ls='--');ax[0].set(xlabel='Iteration',ylabel='Scaled residual',title='Actual iterative residual');ax[1].set(xlabel='Iteration',ylabel='Pressure-only Cd',title='Force convergence');ax[0].legend(fontsize=7);ax[1].legend(fontsize=7)
        save(fig,out,'25_CFD_convergence','Euler residual and force convergence','An iteration cap is not convergence | nonconverged cases are excluded from validated-load claims')
        fig,ax=plt.subplots(1,2,figsize=(14,5.5))
        for label,row,res in chosen:
            surface=res['surface'];ax[0].plot(values(surface,'x_from_nose_m'),values(surface,'Cp'),label=label)
        ax[0].set(xlabel='Distance aft of nose [m]',ylabel='Cp',title='Computed zero-incidence surface pressure');ax[0].legend(fontsize=8)
        for m in [.3,.98,2,5]:
            for side in ['windward','leeward']:
                rr=[r for r in cp if r['Mach']==m and r['side']==side];ax[1].plot(values(rr,'x_over_capsule_length'),values(rr,'Cp'),ls='-' if side=='windward' else '--',label=f'{m:g} {side}')
        ax[1].set(xlabel='Normalized capsule coordinate',ylabel='Cp proxy',title='Analytical incidence pressure estimate');ax[1].legend(fontsize=7)
        save(fig,out,'26_surface_pressure_and_incidence','Pressure distributions with separate provenance','Left: actual Euler iterates | right: analytic exposure proxy; not computed windward/leeward viscous CFD')
        fig,axes=plt.subplots(2,2,figsize=(14,8))
        for a,(label,row,res) in zip(axes.flat,chosen):
            mesh=res['mesh'];rho=res['U'][...,0];grad=np.hypot(*np.gradient(rho));field=np.log10(grad/max(row['rho_kgm3'],1e-20)+1e-5)
            im=a.pcolormesh(mesh.X,mesh.R,field,cmap='Greys',shading='flat');a.contour(mesh.xc,mesh.rc,res['Mach'],levels=[1],colors=[C[1]],linewidths=.8)
            a.set(xlim=(-3,25),ylim=(0,10),xlabel='Distance from nose [m]',ylabel='Radius [m]',title=label+' | density-gradient diagnostic');plt.colorbar(im,ax=a,label='log10 normalized grid gradient',shrink=.8)
        save(fig,out,'31_shock_and_sonic_diagnostics','Computed pressure-wave and sonic diagnostics','Density gradient is a numerical schlieren diagnostic, not a physical optical signal | orange: local Mach 1 contour when present')
    if study:
        fig,ax=plt.subplots(1,2,figsize=(12,5));cells=np.array([r['Nx']*r['Nr'] for r in study]);ax[0].plot(cells,values(study,'CD_pressure_only'),'o-');ax[0].set(xlabel='Cell count',ylabel='Pressure-only Cd',title='Observed grid sensitivity')
        ax[1].semilogy(cells,values(study,'residual'),'o-');ax[1].set(xlabel='Cell count',ylabel='Final residual',title='Convergence status of grid study')
        save(fig,out,'27_CFD_grid_sensitivity','Spatial discretization sensitivity','No grid-independent coefficient or GCI claimed unless iterative and asymptotic convergence are demonstrated')


def _elements(c,out):
    fig,ax=plt.subplots(1,3,figsize=(15,5.5));a=ax[0]
    for i in np.arange(-1,2,.3):
        a.plot([0,2],[i,i],color=C[0],lw=1)
        a.plot([0,2],[i,i+2*np.sqrt(3)],color=C[0],lw=1);a.plot([0,2],[i,i-2*np.sqrt(3)],color=C[0],lw=1)
    a.set(xlim=(0,2),ylim=(0,2),aspect='equal',title='Isogrid topology (schematic)',xlabel='Unwrapped skin [m]',ylabel='Axial skin coordinate [m]')
    a=ax[1];a.add_patch(Rectangle((0,0),2,4,fill=False,lw=2));
    for zz in [.8,2,3.2]:a.plot([0,.35],[zz,zz],color=C[1],lw=5);a.plot([1.65,2],[zz,zz],color=C[1],lw=5)
    a.set(xlim=(-.2,2.2),ylim=(-.2,4.2),aspect='equal',title='Ring baffle meridian (schematic)',xlabel='Tank section',ylabel='Height')
    a=ax[2];a.add_patch(Circle((0,0),2,fill=False,lw=2));
    for th in np.arange(0,2*np.pi,2*np.pi/10):a.plot([0,1.7*np.cos(th)],[0,1.7*np.sin(th)],color=C[0],lw=3)
    a.set(xlim=(-2.2,2.2),ylim=(-2.2,2.2),aspect='equal',title='Thrust-frame load paths (schematic)',xlabel='x',ylabel='y')
    save(fig,out,'28_structural_elements','Structural concepts and load paths','Schematics only; rib geometry, attachments, feedline cutouts and buckling modes need detailed finite-element design')


def animations(c,out,ascent,rv,flows,video=True):
    style();folder=out/'gifs';folder.mkdir(exist_ok=True);videos=out/'videos';videos.mkdir(exist_ok=True);status=[]
    def finish(anim,stem):
        anim.save(folder/f'{stem}.gif',writer=animation.PillowWriter(fps=12),dpi=85)
        if video and animation.writers.is_available('ffmpeg'):
            anim.save(videos/f'{stem}.mp4',writer=animation.FFMpegWriter(fps=12,bitrate=1800),dpi=110);status.append(dict(animation=stem,gif=True,mp4=True))
        else:status.append(dict(animation=stem,gif=True,mp4=False))
    fig,ax=plt.subplots(1,2,figsize=(11,5));t=values(ascent,'t_s');h=values(ascent,'h_m')/1000;v=values(ascent,'speed_inertial_mps')/1000
    for a,data,title,unit in [(ax[0],h,'Ascent altitude','km'),(ax[1],v,'Inertial speed','km/s')]:a.plot(t,data,color='#d8e1e7');a.set(xlim=(0,t[-1]),ylim=(0,max(data)*1.1),title=title,xlabel='Time [s]',ylabel=unit)
    l1,=ax[0].plot([],[]);p1,=ax[0].plot([],[],'o');l2,=ax[1].plot([],[]);p2,=ax[1].plot([],[],'o');caption=fig.suptitle('TARA_LV ascent')
    indices=np.linspace(0,len(t)-1,90).astype(int)
    def update(k):
        i=indices[k];l1.set_data(t[:i+1],h[:i+1]);p1.set_data([t[i]],[h[i]]);l2.set_data(t[:i+1],v[:i+1]);p2.set_data([t[i]],[v[i]]);caption.set_text(f'TARA_LV | t={t[i]:.0f} s | Stage {ascent[i]["stage"]}')
    finish(animation.FuncAnimation(fig,update,frames=len(indices),interval=80),'ascent_time_history');plt.close(fig)
    app=rv['approach'];fig,ax=plt.subplots(figsize=(9,4));ax.set(xlim=(-1100,80),ylim=(-100,100),xlabel='Along-track separation [m]',ylabel='Radial separation [m]',title='TARA_LV prescribed final approach');ax.plot(0,0,'s',ms=12,color=C[1],label='BAS');dot,=ax.plot([],[],'o',ms=9,label='Crew');ax.legend();label=ax.text(.02,.9,'',transform=ax.transAxes)
    ids=np.linspace(0,len(app)-1,80).astype(int)
    def up(k):
        r=app[ids[k]];dot.set_data([r['alongtrack_m']],[r['radial_m']]);label.set_text(f'Range {r["relative_distance_m"]:.1f} m | prescribed path')
    finish(animation.FuncAnimation(fig,up,frames=len(ids),interval=80),'prescribed_docking_approach');plt.close(fig)
    # 3-D state replay with explicitly nonuniform compression by mission phase.
    from .orbits import propagate
    fig=plt.figure(figsize=(11,6));ax=fig.add_subplot(121,projection='3d');panel=fig.add_subplot(122)
    u=np.linspace(0,2*np.pi,35);v=np.linspace(0,np.pi,20);er=RE/1000
    ax.plot_surface(er*np.outer(np.cos(u),np.sin(v)),er*np.outer(np.sin(u),np.sin(v)),er*np.outer(np.ones_like(u),np.cos(v)),color='#89b7cd',alpha=.45,linewidth=0)
    aindices=np.unique(np.linspace(0,len(ascent)-1,25).astype(int));arows=[ascent[k] for k in aindices]
    unique={float(r['t_s']):r for r in arows};at=np.array(sorted(unique));arows=[unique[t] for t in at]
    bas_a=propagate(np.array(c['bas_initial_state']),at)
    frames=[(r,b[:3],'ascent') for r,b in zip(arows,bas_a)]
    for phase,count in [('two_period_coast',20),('J2_phasing',45),('prescribed_docking_approach',20)]:
        rr=[r for r in rv['history'] if r['phase']==phase]
        for i in np.linspace(0,len(rr)-1,count).astype(int):
            r=rr[i];frames.append((r,np.array([r['BAS_rx_m'],r['BAS_ry_m'],r['BAS_rz_m']]),phase))
    pos=np.array([[r['rx_m'],r['ry_m'],r['rz_m']] for r,b,p in frames])/1000
    dense=ascent+rv['history'];dense_t=values(dense,'t_s');dense_pos=np.c_[values(dense,'rx_m'),values(dense,'ry_m'),values(dense,'rz_m')]/1000
    ax.plot(dense_pos[:,0],dense_pos[:,1],dense_pos[:,2],color='#c8d8df',lw=.8);crew,=ax.plot([],[],[],'o',color=C[1],label='Crew',ms=6);station,=ax.plot([],[],[],'s',color=C[0],label='BAS',ms=5);trail,=ax.plot([],[],[],color=C[1],lw=1)
    ax.set(xlim=(-7500,7500),ylim=(-7500,7500),zlim=(-7500,7500),xlabel='ECI x [km]',ylabel='ECI y [km]',zlabel='ECI z [km]');ax.set_box_aspect([1,1,1]);ax.legend(loc='upper left',fontsize=8)
    panel.axis('off');clock=panel.text(.05,.8,'',transform=panel.transAxes,fontsize=17,va='top');panel.text(.05,.15,'Actual integrated state replay\nTime compression varies by phase.\nTerminal approach is prescribed.\nNo descent simulation.',transform=panel.transAxes,color='#63778b',fontsize=10)
    fig.suptitle('TARA_LV | Launch, two-period coast and BAS rendezvous',fontsize=15,fontweight='bold')
    def stateframe(k):
        r,b,phase=frames[k];p=pos[k];crew.set_data_3d([p[0]],[p[1]],[p[2]]);bb=b/1000;station.set_data_3d([bb[0]],[bb[1]],[bb[2]]);end=np.searchsorted(dense_t,r['t_s'],side='right');trail.set_data_3d(dense_pos[:end,0],dense_pos[:end,1],dense_pos[:end,2]);height=r.get('altitude_m',r.get('h_m',0))/1000
        clock.set_text(f'{phase.replace("_"," ").upper()}\n\nElapsed: {r["t_s"]/3600:.2f} h\nAltitude: {height:.1f} km')
    finish(animation.FuncAnimation(fig,stateframe,frames=len(frames),interval=80),'mission_state_replay');plt.close(fig)
    if flows[1]:
        fig,ax=plt.subplots(figsize=(10,4));label,row,res=flows[1][0];mesh=res['mesh'];im=ax.pcolormesh(mesh.X,mesh.R,res['Mach'],shading='flat',cmap='turbo',vmin=0,vmax=6);plt.colorbar(im,ax=ax,label='Local Mach');ax.set(xlim=(-3,30),ylim=(0,10),xlabel='Distance from nose [m]',ylabel='Radius [m]');title=ax.set_title(label)
        # Each frame is a solved condition, not physical flow-time evolution.
        def flowframe(k):
            name,r,s=flows[1][k//12];im.set_array(s['Mach']);title.set_text(f'{name} | '+('converged' if s['summary']['converged'] else 'iteration-limited')+' | condition sweep')
        finish(animation.FuncAnimation(fig,flowframe,frames=len(flows[1])*12,interval=80),'computed_Euler_condition_sweep');plt.close(fig)
    return status
