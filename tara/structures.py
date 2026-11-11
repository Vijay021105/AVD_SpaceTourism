"""Flight section-cut loads, metallic shell/joint screens and dry-mass audit."""
import numpy as np
from .io import G0,write_csv


def hardware_mass(c):
    """Same item estimates used by staging dry floors and the final mass audit."""
    st=c['structures'];mat=c['materials'][0];R=c['diameter_m']/2;items=[]
    for s in c['stages']:
        sid=int(s['stage']);ts=[t for t in c['tanks'] if t['stage']==sid]
        for t in ts:
            items.extend([dict(stage=sid,component=t['name']+'_shell',mass_kg=t['shell_mass_kg']),
                dict(stage=sid,component=t['name']+'_insulation',mass_kg=t['insulation_mass_kg']),
                dict(stage=sid,component=t['name']+'_baffles',mass_kg=st['baffle_count']*2*np.pi*t['inner_radius_m']*st['baffle_width_m']*st['baffle_thickness_m']*mat['density_kgm3'])])
        shell=sum(t['shell_mass_kg'] for t in ts)
        items.extend([dict(stage=sid,component='stiffeners_isogrid_allowance',mass_kg=shell*st['stiffener_mass_fraction']),
            dict(stage=sid,component='engines',mass_kg=s['engine_count']*s['engine']['engine_mass_kg']),
            dict(stage=sid,component='thrust_frame',mass_kg=s['engine_count']*s['engine']['thrust_vac_N']*st['engine_frame_mass_per_N']),
            dict(stage=sid,component='interstage_and_engine_skirt',mass_kg=2*np.pi*R*(s['interstage_m']+s['engine_bay_m'])*st['effective_skin_m']*mat['density_kgm3']),
            dict(stage=sid,component='joint_rings',mass_kg=2*2*np.pi*R*st['joint_mass_per_m']),
            dict(stage=sid,component='avionics_feed_pressurization_allowance',mass_kg=s['fixed_dry_allowance_kg'])])
    return items


def mass_distribution(c,row,z):
    """Nodal lumped masses with exact component sums; attached stages only."""
    mass=np.zeros_like(z);stage=int(row['stage'])
    def put(total,lo,hi):
        weights=((z>=lo)&(z<=hi)).astype(float)
        if not weights.any():weights[np.argmin(abs(z-(lo+hi)/2))]=1
        mass[:]+=total*weights/weights.sum()
    for s in c['stages']:
        sid=int(s['stage'])
        if sid<stage:continue
        e=s['engine'];eng=s['engine_count']*e['engine_mass_kg']
        put(eng,s['aft_z_m'],s['aft_z_m']+s['engine_bay_m'])
        put(max(0,s['dry_kg']-eng),s['aft_z_m'],s['forward_z_m'])
        frac=(1-(1-s['reserve_fraction'])*row['active_prop_fraction']) if sid==stage else 1.
        for t in c['tanks']:
            if t['stage']==sid:
                # Settled liquid approximated by a uniform column, shrinking from forward end.
                put(t['loaded_mass_kg']*frac,t['z0_m'],t['z0_m']+max(.01,(t['z1_m']-t['z0_m'])*frac))
    put(c['mission']['payload_kg'],c['stages'][2]['forward_z_m'],c['tip_z_m'])
    if abs(mass.sum()-row['mass_kg'])/row['mass_kg']>1e-7:raise ValueError('Distributed mass and ascent mass disagree')
    return mass


def run(c,ascent,out):
    st=c['structures'];mat=c['materials'][0];R=c['diameter_m']/2;L=c['tip_z_m'];z=np.linspace(0,L,301)
    # All flight samples enter the envelope; each stage has its own surviving stack.
    Nx=np.zeros_like(z);Vx=np.zeros_like(z);Mx=np.zeros_like(z);case_rows=[];selected=[]
    critical={int(np.argmax([r[k] for r in ascent])) for k in ['q_Pa','qalpha_Pa_rad','proper_g','thrust_N']}
    critical.update(i for i,r in enumerate(ascent) if i==len(ascent)-1 or (i<len(ascent)-1 and ascent[i+1]['stage']!=r['stage']))
    inertia_error=0.
    for n,row in enumerate(ascent):
        mass=mass_distribution(c,row,z);s=c['stages'][int(row['stage'])-1];base=s['aft_z_m'];mask=z>=base
        zcg=np.sum(z*mass)/mass.sum();aero=np.exp(-((z-(L-row['Xcp_from_nose_m']))/(.07*L))**2)*mask
        aero/=max(aero.sum(),1e-20)
        ext=row['normal_N']*aero
        arm=base-zcg
        tvc=-np.sum(ext*(z-zcg))/arm if abs(arm)>1e-9 else 0.
        external=ext.copy();external[np.argmin(abs(z-base))]+=tvc
        lateral_acc=external.sum()/mass.sum();net=external-mass*lateral_acc
        shear=np.cumsum(net[::-1])[::-1]
        moment=np.array([np.sum(net[k:]*(z[k:]-zz)) for k,zz in enumerate(z)])
        drag_above=np.cumsum((row['drag_N']*aero)[::-1])[::-1]
        axial=np.cumsum(mass[::-1])[::-1]*row['proper_axial_mps2']+drag_above
        axial[~mask]=0;shear[~mask]=0;moment[~mask]=0
        Nx=np.maximum(Nx,abs(axial));Vx=np.maximum(Vx,abs(shear));Mx=np.maximum(Mx,abs(moment))
        inertia_error=max(inertia_error,abs(net.sum())/max(abs(external).sum(),1))
        if n in critical:
            case=f'LC_{n:04d}_S{row["stage"]}'
            for zz,nn,vv,mm in zip(z,axial,shear,moment):selected.append(dict(case=case,t_s=row['t_s'],stage=row['stage'],z_m=zz,axial_N=nn,shear_N=vv,bending_Nm=mm))
            case_rows.append(dict(case=case,t_s=row['t_s'],stage=row['stage'],q_Pa=row['q_Pa'],proper_g=row['proper_g'],
                TVC_lateral_trim_N=tvc,max_axial_N=max(abs(axial)),max_shear_N=max(abs(shear)),max_bending_Nm=max(abs(moment)),
                lateral_equilibrium_relative_error=abs(net.sum())/max(abs(external).sum(),1)))
    envelope=[dict(z_m=zz,axial_N=nn,shear_N=vv,bending_Nm=mm) for zz,nn,vv,mm in zip(z,Nx,Vx,Mx)]
    margins=[]
    for s in c['stages']:
        sid=int(s['stage']);ts=[t for t in c['tanks'] if t['stage']==sid]
        for t in ts:
            mask=(z>=t['z0_m'])&(z<=t['z1_m']);n=max(Nx[mask]);moment=max(Mx[mask]);v=max(Vx[mask]);r=t['inner_radius_m'];th=t['wall_m']
            A=2*np.pi*r*th;I=np.pi*r**3*th
            hoop=t['design_pressure_Pa']*r/th
            sigma=hoop/2-n/A-moment*r/I;sigplus=hoop/2-n/A+moment*r/I
            tau=v/(np.pi*r*th)
            vm=max(np.sqrt(hoop**2+sigma**2-hoop*sigma+3*tau*tau),np.sqrt(hoop**2+sigplus**2-hoop*sigplus+3*tau*tau))*st['load_uncertainty']
            comp=max(0,-sigma)*st['load_uncertainty']
            buck=st['buckling_knockdown']*mat['E_Pa']/np.sqrt(3*(1-mat['nu']**2))*th/r
            margins.append(dict(component=t['name'],stage=sid,wall_m=th,hoop_stress_Pa=hoop,
                limit_von_mises_Pa=vm,yield_margin=mat['yield_Pa']/(st['yield_fos']*vm)-1,
                ultimate_margin=mat['ultimate_Pa']/(st['ultimate_fos']*vm)-1,
                buckling_margin=(buck/(st['ultimate_fos']*comp)-1) if comp>0 else '',
                buckling_status='compression_screen' if comp>0 else 'inactive_no_compression',
                note='Envelope components may occur at different times; conservative combination; no pressure buckling credit; no fatigue or fracture assessment'))
        idx=np.argmin(abs(z-s['forward_z_m']));bolt=(Nx[idx]+abs(Mx[idx])/st['joint_lever_arm_m'])/st['joint_bolts']
        margins.append(dict(component=f'S{sid}_joint',stage=sid,limit_bolt_force_N=bolt,
            ultimate_margin=st['bolt_allowable_N']/(max(bolt,1)*st['ultimate_fos']*st['load_uncertainty'])-1,
            note='Equal bolt load plus simple moment couple; local fitting/prying/contact not modeled'))
    masslist=hardware_mass(c)
    closure=[]
    for s in c['stages']:
        subtotal=sum(r['mass_kg'] for r in masslist if r['stage']==s['stage'])
        closure.append(dict(stage=s['stage'],dry_budget_kg=s['dry_kg'],estimated_hardware_kg=subtotal,
            unallocated_dry_budget_kg=s['dry_kg']-subtotal,closure_pass=bool(s['dry_kg']>=subtotal-.1),
            scope='Conceptual item mass screen; no missing subsystem can be assumed massless'))
    write_csv(out/'data/loads_envelope.csv',envelope);write_csv(out/'data/load_cases.csv',case_rows)
    write_csv(out/'data/loads_by_case.csv',selected);write_csv(out/'data/structural_margins.csv',margins)
    write_csv(out/'data/hardware_mass_breakdown.csv',masslist);write_csv(out/'data/dry_mass_closure.csv',closure)
    return envelope,margins,masslist,closure,dict(lateral_force_equilibrium_relative_error=inertia_error)
