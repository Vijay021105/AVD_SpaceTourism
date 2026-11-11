"""ECI Cartesian astrodynamics with J2 and one UTC/GMST convention."""
from datetime import datetime,timezone
import numpy as np
from scipy.integrate import solve_ivp
from .io import MU,RE,J2,OMEGA


def gravity(r,j2=True):
    r=np.asarray(r);norm=np.linalg.norm(r,axis=-1);a=-MU*r/norm[...,None]**3
    if j2:
        zz=(r[...,2]/norm)**2;fac=1.5*J2*MU*RE**2/norm**5
        a+=fac[...,None]*r*np.stack([5*zz-1,5*zz-1,5*zz-3],axis=-1)
    return a


def propagate(y,times,j2=True,rtol=2e-10):
    times=np.asarray(times,dtype=float)
    if len(times)<2 or np.any(np.diff(times)<=0):raise ValueError('Propagation times must strictly increase')
    sol=solve_ivp(lambda t,y:np.r_[y[3:],gravity(y[:3],j2)],(times[0],times[-1]),y,
        t_eval=times,rtol=rtol,atol=1e-5,max_step=60,method='DOP853')
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y.T


def gmst(epoch):
    dt=datetime.fromisoformat(epoch.replace('Z','+00:00')).astimezone(timezone.utc)
    jd=dt.timestamp()/86400+2440587.5;d=jd-2451545
    return np.deg2rad((280.46061837+360.98564736629*d)%360)


def rotation_z(a):
    co,si=np.cos(a),np.sin(a)
    return np.array([[co,-si,0],[si,co,0],[0,0,1.]])


def frame(r):
    rhat=r/np.linalg.norm(r);east=np.cross([0,0,1],rhat);east/=np.linalg.norm(east);north=np.cross(rhat,east)
    return rhat,east,north


def rtn(y):
    r,v=y[:3],y[3:];R=r/np.linalg.norm(r);N=np.cross(r,v);N/=np.linalg.norm(N);T=np.cross(N,R)
    return np.c_[R,T,N]


def elements(y):
    r,v=np.asarray(y[:3]),np.asarray(y[3:6]);R=np.linalg.norm(r);h=np.cross(r,v);H=np.linalg.norm(h)
    ev=np.cross(v,h)/MU-r/R;e=np.linalg.norm(ev);energy=np.dot(v,v)/2-MU/R
    a=-MU/(2*energy);inc=np.arccos(np.clip(h[2]/H,-1,1));node=np.cross([0,0,1],h)
    raan=np.arctan2(node[1],node[0])%(2*np.pi)
    n=node/max(np.linalg.norm(node),1e-30);p=np.cross(h/H,n)
    u=np.arctan2(np.dot(r,p),np.dot(r,n))%(2*np.pi)
    return dict(a_m=a,eccentricity=e,inclination_deg=np.rad2deg(inc),RAAN_deg=np.rad2deg(raan),
        argument_latitude_deg=np.rad2deg(u),perigee_altitude_m=a*(1-e)-RE,apogee_altitude_m=a*(1+e)-RE,
        specific_energy_Jkg=energy)


def setup(c):
    m=c['mission'];lat=np.deg2rad(m['site_lat_deg']);lon=np.deg2rad(m['site_lon_deg'])+gmst(m['epoch_utc'])
    # Spherical launch site; flattening omitted consistently from this conceptual model.
    r=(RE+m['site_height_m'])*np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    rr,east,north=frame(r);inc=np.deg2rad(m['inclination_deg'])
    az=np.pi-np.arcsin(np.cos(inc)/np.cos(lat));tangent=np.sin(az)*east+np.cos(az)*north
    normal=np.cross(rr,tangent);normal/=np.linalg.norm(normal)
    node=np.cross([0,0,1],normal);node/=np.linalg.norm(node);other=np.cross(normal,node)
    u=np.deg2rad(m['bas_u0_deg']);R=RE+m['orbit_altitude_m']
    bas=np.r_[R*(np.cos(u)*node+np.sin(u)*other),np.sqrt(MU/R)*(-np.sin(u)*node+np.cos(u)*other)]
    c['plane_normal']=normal;c['site_r_eci']=r;c['bas_initial_state']=bas
    c['node_precession_radps']=-1.5*J2*np.sqrt(MU/R**3)*(RE/R)**2*np.cos(inc)
    c['launch_azimuth_deg']=np.rad2deg(az)
    return r,np.cross([0,0,OMEGA],r),bas


def groundtrack(c,t,y):
    pos=rotation_z(-(gmst(c['mission']['epoch_utc'])+OMEGA*t))@np.asarray(y[:3])
    return np.rad2deg(np.arctan2(pos[1],pos[0])),np.rad2deg(np.arctan2(pos[2],np.hypot(pos[0],pos[1])))
