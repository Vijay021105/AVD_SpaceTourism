"""ISA to 84.852 geopotential km, explicit thermospheric extrapolation and wind."""
import numpy as np
from .io import G0,RE

H=np.array([0,11000,20000,32000,47000,51000,71000,84852.])
L=np.array([-.0065,0,.001,.0028,0,-.0028,-.002,0.])
T=np.empty(8);P=np.empty(8);T[0]=288.15;P[0]=101325
for j in range(1,8):
    dh=H[j]-H[j-1];T[j]=T[j-1]+L[j-1]*dh
    P[j]=P[j-1]*(T[j]/T[j-1])**(-G0/(287.05287*L[j-1])) if L[j-1] else P[j-1]*np.exp(-G0*dh/(287.05287*T[j-1]))


def atmosphere(alt):
    alt=max(-500.,float(alt));h=RE*alt/(RE+alt);j=int(np.clip(np.searchsorted(H,h)-1,0,7))
    if h<=H[-1]:
        t=T[j]+L[j]*(h-H[j]);p=P[j]*(t/T[j])**(-G0/(287.05287*L[j])) if L[j] else P[j]*np.exp(-G0*(h-H[j])/(287.05287*T[j]))
        rho=p/(287.05287*t)
    else:
        # Log-density interpolation: explicitly a screening atmosphere, not NRLMSISE.
        hh=[86000,100000,120000,150000,200000,300000,400000,600000,1000000]
        rr=[6.96e-6,5.6e-7,2.2e-8,2.1e-9,2.8e-10,2.4e-11,3.7e-12,1.1e-13,3e-15]
        rho=np.exp(np.interp(alt,hh,np.log(rr)));t=np.interp(alt,[86000,120000,200000,400000],[187,360,850,1000]);p=rho*287.05287*t
    mu=1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
    return dict(rho_kgm3=rho,T_K=t,p_Pa=p,a_mps=np.sqrt(1.4*287.05287*t),mu_Pas=mu,
                atmosphere_valid=alt<86000,mean_free_path_m=mu/max(p,1e-30)*np.sqrt(np.pi*287.05287*t/2))


def wind(c,alt):
    h=[r['altitude_m'] for r in c['wind']]
    return np.array([np.interp(alt,h,[r[k] for r in c['wind']]) for k in ['east_mps','north_mps']])
