"""One geometry generated from mass, density, ullage and exterior envelope.
z: aft to nose. CFD x: nose to aft. Tank heads: two 2:1 ellipsoids.
"""
import numpy as np
from scipy.optimize import minimize
from .io import G0,write_csv


def tank_dimensions(volume,r,ratio=.5):
    h=ratio*r; heads=4*np.pi*r*r*h/3
    if volume<heads:
        r=(3*volume/(4*np.pi*ratio))**(1/3);h=ratio*r;barrel=0.
    else:barrel=(volume-heads)/(np.pi*r*r)
    return r,h,barrel


def build(c,out=None):
    g=c['geometry'];st=c['structures'];mat=c['materials'][0]
    fluids={r['fluid']:r for r in c['fluids']}; D=c['mission']['diameter_m']
    z=0.;parts=[];tanks=[];stations=[]
    for j,s in enumerate(c['stages']):
        sid=j+1;s['aft_z_m']=z;stations.append(z)
        parts.append(dict(component=f'S{sid}_engine_bay',stage=sid,z0_m=z,z1_m=z+s['engine_bay_m'],radius_m=D/2,role='thrust frame and engine clearance'))
        z+=s['engine_bay_m']
        e=s['engine']; fuel_mass=s['propellant_kg']/(1+e['OF'])
        for fluid,mass in [(e['fuel'],fuel_mass),('LOX',fuel_mass*e['OF'])]:
            f=fluids[fluid];ins=g['rp1_insulation_m'] if fluid=='RP1' else g['cryo_insulation_m']
            bare_r=D/2-g['skin_m']-g['radial_clearance_m']-ins
            # Conservative pressure sizing includes full-column at 1.8 g.
            wall=st['shell_min_m']
            for _ in range(12):
                volume=mass/f['density_kgm3']/(1-g['ullage_fraction'])
                ri,head,barrel=tank_dimensions(volume,bare_r-wall,g['dome_depth_over_radius'])
                length=barrel+2*head
                pressure=f['pressure_abs_Pa']+f['density_kgm3']*1.8*G0*length
                wall_new=max(st['shell_min_m'],st['pressure_fos']*pressure*ri/mat['yield_Pa'])
                if abs(wall_new-wall)<1e-9:break
                wall=.5*(wall+wall_new)
            wall=wall_new
            ri,head,barrel=tank_dimensions(volume,bare_r-wall,g['dome_depth_over_radius']);length=barrel+2*head
            # Oblate spheroid area + barrel; two half heads form one spheroid.
            ecc=np.sqrt(max(0,1-(head/ri)**2))
            head_area=2*np.pi*ri*ri*(1+(1-ecc*ecc)/max(ecc,1e-12)*np.arctanh(min(ecc,1-1e-12)))
            area=2*np.pi*ri*barrel+head_area
            r=dict(name=f'S{sid}_{fluid}',stage=sid,fluid=fluid,loaded_mass_kg=mass,liquid_density_kgm3=f['density_kgm3'],
                volume_m3=volume,liquid_volume_m3=mass/f['density_kgm3'],ullage_fraction=g['ullage_fraction'],
                inner_radius_m=ri,wall_m=wall,bare_outer_radius_m=ri+wall,insulation_m=ins,
                installed_outer_radius_m=ri+wall+ins,z0_m=z,z1_m=z+length,
                head_depth_m=head,barrel_m=barrel,area_m2=area,shell_mass_kg=area*wall*mat['density_kgm3'],
                insulation_mass_kg=area*ins*c['thermal']['cryo_foam_density_kgm3'],
                design_pressure_Pa=pressure,pressure_abs_Pa=f['pressure_abs_Pa'],temperature_K=f['temperature_K'],
                capacity_mass_kg=volume*(1-g['ullage_fraction'])*f['density_kgm3'])
            tanks.append(r);parts.append(dict(component=r['name'],stage=sid,z0_m=z,z1_m=z+length,radius_m=ri+wall+ins,role='propellant tank'))
            z+=length+g['tank_gap_m']
        parts.append(dict(component=f'S{sid}_interstage',stage=sid,z0_m=z,z1_m=z+s['interstage_m'],radius_m=D/2,role='load carrying skirt; rigid-joint screen'))
        z+=s['interstage_m'];s['forward_z_m']=z;s['length_m']=z-s['aft_z_m']
    c['capsule_base_z_m']=z+g['trunk_length_m'];
    parts.append(dict(component='trunk',stage=4,z0_m=z,z1_m=c['capsule_base_z_m'],radius_m=g['capsule_base_radius_m'],role='spacecraft trunk'))
    z=c['capsule_base_z_m']
    parts.append(dict(component='crew_capsule',stage=4,z0_m=z,z1_m=z+g['capsule_length_m'],radius_m=g['capsule_base_radius_m'],role='conceptual crew outer mold line'))
    c['tip_z_m']=z+g['capsule_length_m'];c['stage_aft_z_m']=stations;c['diameter_m']=D
    c['tanks']=tanks;c['parts']=parts
    if out:
        write_csv(out/'data/tanks.csv',tanks);write_csv(out/'data/components.csv',parts)
        zz=np.linspace(0,c['tip_z_m'],1000)
        write_csv(out/'data/outer_profile.csv',[dict(z_m=x,radius_m=y) for x,y in zip(zz,radius(c['tip_z_m']-zz,c))])
        # Portable DXF meridian as one closed polyline. Units: metres ($INSUNITS=6).
        pts=list(zip(radius(c['tip_z_m']-zz,c),zz))+list(zip(-radius(c['tip_z_m']-zz[::-1],c),zz[::-1]))
        text='0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n6\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n0\nLWPOLYLINE\n8\nOUTLINE\n90\n'+str(len(pts))+'\n70\n1\n'
        text+=''.join(f'10\n{x:.8f}\n20\n{y:.8f}\n' for x,y in pts)+'0\nENDSEC\n0\nEOF\n'
        (out/'data/vehicle_profile.dxf').write_text(text)
        # Closed axisymmetric visualization mesh in metres (stage splits remain
        # in components.csv). Supplier engine internals are deliberately absent.
        zmesh=np.linspace(0,c['tip_z_m']-.001,180);rmesh=radius(c['tip_z_m']-zmesh,c);nt=48
        vertices=[(rr*np.cos(th),rr*np.sin(th),zz) for zz,rr in zip(zmesh,rmesh) for th in np.arange(nt)*2*np.pi/nt]
        base=len(vertices);vertices.append((0,0,0));tip=len(vertices);vertices.append((0,0,c['tip_z_m']));faces=[]
        for j in range(len(zmesh)-1):
            for k in range(nt):
                a=j*nt+k;b=j*nt+(k+1)%nt;d=(j+1)*nt+k;e=(j+1)*nt+(k+1)%nt
                faces.extend([(a,b,d),(b,e,d)])
        for k in range(nt):faces.extend([(base,(k+1)%nt,k),((len(zmesh)-1)*nt+k,(len(zmesh)-1)*nt+(k+1)%nt,tip)])
        obj='# TARA_LV concept exterior, metres; z from aft; no installation qualification\n'
        obj+=''.join('v %.7f %.7f %.7f\n'%v for v in vertices)
        obj+=''.join('f %d %d %d\n'%(a+1,b+1,d+1) for a,b,d in faces)
        (out/'data/vehicle_exterior.obj').write_text(obj)
    return tanks,parts


def radius(x,c):
    z=c['tip_z_m']-np.asarray(x,dtype=float);r=np.full_like(z,c['diameter_m']/2)
    base=c['capsule_base_z_m'];capr=c['geometry']['capsule_base_radius_m'];L=c['geometry']['capsule_length_m']
    forward=c['stages'][2]['forward_z_m'];adapter=c['stages'][2]['interstage_m']
    blend=np.clip((z-(forward-adapter))/adapter,0,1)
    r=np.where(z>=forward-adapter,c['diameter_m']/2*(1-blend)+capr*blend,r)
    # Cubic radius-squared contour: stated nose curvature, base radius and
    # zero base slope are all satisfied. Thermal R_n is thus the actual OML R_n.
    nose=c['geometry']['nose_radius_m'];s=np.clip(c['tip_z_m']-z,0,L)
    aa=(3*capr**2-4*nose*L)/L**2;bb=(-2*capr**2+2*nose*L)/L**3
    cap=np.sqrt(np.maximum(0,2*nose*s+aa*s*s+bb*s*s*s))
    r=np.where(z>=base,cap,r)
    return np.where(z>=c['tip_z_m'],0,r)


def engine_layout(c,out):
    e=c['stages'][0]['engine'];d=e['nozzle_exit_diameter_m'];gap=c['geometry']['engine_clearance_m']
    # Feasibility search of 10 circles. This is one local packing solution, not a proof of impossibility.
    theta=np.arange(9)*2*np.pi/9
    xy=np.vstack([[0,0],np.c_[1.5*np.cos(theta),1.5*np.sin(theta)]])
    def con(v):
        p=v[:-1].reshape(10,2);R=v[-1]
        return np.r_[R-d/2-gap-np.linalg.norm(p,axis=1),
            [np.linalg.norm(p[a]-p[b])-d-gap for a in range(10) for b in range(a)]]
    sol=minimize(lambda v:v[-1],np.r_[xy.ravel(),2.2],constraints=[dict(type='ineq',fun=con)],
                 bounds=[(-3,3)]*20+[(1,4)],method='SLSQP',options=dict(maxiter=500,ftol=1e-9))
    xy=sol.x[:-1].reshape(10,2);R=sol.x[-1]
    dynamic=d/2+e['length_m']*np.sin(np.deg2rad(e['gimbal_deg']))+gap
    rows=[dict(stage=1,engine=k+1,x_m=x,y_m=y,nozzle_radius_m=d/2,gimbal_envelope_radius_m=dynamic) for k,(x,y) in enumerate(xy)]
    for sid,xs in [(2,[-.7,.7]),(3,[0])]:
        ee=c['stages'][sid-1]['engine']
        rows.extend(dict(stage=sid,engine=k+1,x_m=x,y_m=0,nozzle_radius_m=ee['nozzle_exit_diameter_m']/2,
            gimbal_envelope_radius_m=ee['nozzle_exit_diameter_m']/2+ee['length_m']*np.sin(np.deg2rad(ee['gimbal_deg']))) for k,x in enumerate(xs))
    write_csv(out/'data/engine_layout.csv',rows)
    radial=max(np.linalg.norm(xy,axis=1)+dynamic)
    pair=min(np.linalg.norm(xy[a]-xy[b])-2*dynamic for a in range(10) for b in range(a))
    return rows,dict(static_packing_diameter_m=2*R,static_feasible_in_search=bool(2*R<=c['diameter_m']),
        packing_optimizer_success=bool(sol.success),gimbal_envelope_diameter_m=2*radial,
        pairwise_gimbal_clearance_m=pair,gimbal_installation_pass=bool(2*radial<=c['diameter_m'] and pair>=0),
        note='Exit diameter length and gimbal assumptions need supplier data; local packing is not a global impossibility proof')


class RingMesh:
    """Exact axisymmetric ring volumes and vector face areas; metric-preserving."""
    def __init__(self,c,stage=1,nx=140,nr=32,perturb=None):
        if not 1<=stage<=3 or nx<12 or nr<2:raise ValueError('Invalid ring mesh size or stage')
        L=c['tip_z_m']-c['stage_aft_z_m'][stage-1];up=c['cfd']['upstream_m'];outer=c['cfd']['outer_radius_m']
        nup=min(nx-2,max(6,int(nx*up/(L+up))))
        self.x=np.r_[np.linspace(-up,0,nup+1)[:-1],np.linspace(0,L,nx-nup+1)]
        rb=radius(self.x,c)
        if perturb is not None:rb+=np.interp(self.x,perturb[0],perturb[1],left=0,right=0)
        if np.max(rb)>=outer or np.min(rb)<0:raise ValueError('Mesh body outside domain')
        self.X=np.broadcast_to(self.x[:,None],(nx+1,nr+1)).copy()
        eta=(np.exp(2*np.linspace(0,1,nr+1))-1)/(np.exp(2)-1)
        self.R=rb[:,None]+(outer-rb[:,None])*eta
        self.xc=(self.X[:-1,:-1]+self.X[1:,:-1]+self.X[:-1,1:]+self.X[1:,1:])/4
        a,b=self.R[:-1,:-1],self.R[1:,:-1];A,B=self.R[:-1,1:],self.R[1:,1:];dx=np.diff(self.x)[:,None]
        self.V=np.pi*dx*((A*A+A*B+B*B)-(a*a+a*b+b*b))/3
        self.area2=dx*((A-a)+(B-b))/2;self.rc=(a+b+A+B)/4
        self.Sx=np.stack([np.pi*(self.R[:,1:]**2-self.R[:,:-1]**2),np.zeros((nx+1,nr))],axis=-1)
        rm=(self.R[:-1]+self.R[1:])/2
        self.Sr=2*np.pi*rm[:,:,None]*np.stack([-np.diff(self.R,axis=0),np.broadcast_to(dx,(nx,nr+1))],axis=-1)
        self.nx,self.nr=nx,nr
        if np.any(self.V<=0):raise ValueError('Nonpositive cell volume')
    def check(self):
        balance=self.Sx[1:]-self.Sx[:-1]+self.Sr[:,1:]-self.Sr[:,:-1]
        expected=np.stack([np.zeros_like(self.area2),2*np.pi*self.area2],axis=-1)
        return float(np.max(abs(balance-expected))/np.max(abs(expected)))
