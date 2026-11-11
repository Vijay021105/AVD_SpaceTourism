"""Documented engineering aerodynamic database, separate from Euler pressure CFD.
No separated-flow or transition field solution is implied by these correlations.
"""
import numpy as np
from .io import write_csv


def coefficients(c,M,alpha=0.,stage=1,rho=1.225,speed=100.,mu=1.8e-5):
    a=c['aerodynamics'];D=c['diameter_m'];length=c['tip_z_m']-c['stage_aft_z_m'][int(stage)-1]
    Re=max(1e3,rho*max(speed,1)*length/mu)
    cf_lam=1.328/np.sqrt(Re)
    cf_turb=.455/np.log10(Re)**2.58/(1+.144*M*M)**.65
    lam_fraction=min(1,a['transition_Re']/Re)
    cf=lam_fraction*cf_lam+(1-lam_fraction)*cf_turb
    friction=cf*(np.pi*D*length)/(np.pi*D*D/4)
    peak=a['wave_peak_cd']*np.exp(-((M-1.15)/.29)**2)
    wave=peak+np.where(M>1.3,.23*(1-np.exp(-(M-1.3)/.35))/max(M/2,1),0.)
    ca=a['base_cd']+a['form_cd']+friction+wave
    cna=a['cnalpha_sub_per_rad']*(1+.45*np.exp(-((M-1.05)/.3)**2))
    cna=cna/(1+.018*max(M-3,0))
    cn=cna*alpha+1.2*np.sin(alpha)*abs(np.sin(alpha))
    cd=ca*np.cos(alpha)+cn*np.sin(alpha)
    cp_nose=c['geometry']['capsule_length_m']*.55
    return dict(Ca=float(ca),Cd=float(cd),Cn=float(cn),Cnalpha_per_rad=float(cna),
        Xcp_from_nose_m=float(cp_nose),Cl_roll=0.,Re_L=Re,Cf_laminar=cf_lam,Cf_turbulent=cf_turb,
        Cd_friction=float(friction),Cd_wave=float(wave),Cd_base=a['base_cd'],
        model='engineering correlation; zero roll from symmetry; no separation resolved')


def database(c,out):
    mm=np.unique(np.r_[np.linspace(.05,.75,10),np.linspace(.8,1.4,25),np.linspace(1.5,6,25),np.linspace(7,c['aerodynamics']['mach_max'],24)])
    rows=[];polars=[]
    for stage in [1,2,3]:
        for m in mm:
            r=coefficients(c,m,np.deg2rad(c['aerodynamics']['alpha_reference_deg']),stage)
            rows.append(dict(stage=stage,Mach=m,alpha_deg=c['aerodynamics']['alpha_reference_deg'],**r))
            for alpha in [-6,-3,0,3,6]:
                polars.append(dict(stage=stage,Mach=m,alpha_deg=alpha,**coefficients(c,m,np.deg2rad(alpha),stage)))
    write_csv(out/'data/aerodynamic_database.csv',rows);write_csv(out/'data/aero_polars.csv',polars)
    return rows,polars
