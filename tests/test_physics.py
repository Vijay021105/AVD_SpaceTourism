"""Independent conservation, exact-solution and limiting-case checks."""
from pathlib import Path
import sys
import copy
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tara.io import load_inputs,MU,RE,G0
from tara.orbits import propagate,elements,gravity,setup
from tara.sizing import ledger,optimize
from tara.geometry import build,RingMesh,radius
from tara.cfd import state,rhs,flux
from tara.thermal import conduction,orifice_flux

ROOT=Path(__file__).resolve().parents[1]


class PhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=load_inputs(ROOT/'inputs');rows,_,_=optimize(cls.c)
        for s,r in zip(cls.c['stages'],rows):s.update(r)
        build(cls.c);cls.c['cfd']=dict(upstream_m=10,outer_radius_m=16)

    def test_two_body_exact_circular_period(self):
        r=RE+400000;v=np.sqrt(MU/r);T=2*np.pi*np.sqrt(r**3/MU);y=np.array([r,0,0,0,v,0])
        end=propagate(y,[0,T],False)[-1]
        self.assertLess(np.linalg.norm(end[:3]-y[:3]),.05)
        self.assertLess(np.linalg.norm(end[3:]-y[3:]),1e-4)

    def test_j2_regresses_node_at_expected_rate(self):
        c=copy.deepcopy(self.c);_,_,y=setup(c);ts=np.linspace(0,86400,500);ys=propagate(y,ts)
        raan=np.unwrap(np.deg2rad([elements(yy)['RAAN_deg'] for yy in ys]));slope=np.polyfit(ts,raan,1)[0]
        self.assertLess(abs(slope/c['node_precession_radps']-1),.02)

    def test_gravity_points_inward(self):
        for r in [[7e6,0,0],[0,0,7e6],[3e6,3e6,6e6]]:self.assertLess(np.dot(r,gravity(np.array(r))),0)

    def test_tank_volumes_and_envelope(self):
        for t in self.c['tanks']:
            volume=np.pi*t['inner_radius_m']**2*t['barrel_m']+4*np.pi*t['inner_radius_m']**2*t['head_depth_m']/3
            self.assertAlmostEqual(volume,t['volume_m3'],places=7)
            self.assertAlmostEqual(t['capacity_mass_kg'],t['loaded_mass_kg'],places=6)
            self.assertLessEqual(t['installed_outer_radius_m'],2)

    def test_ring_metrics_and_static_pressure_balance(self):
        mesh=RingMesh(self.c,1,60,16);self.assertLess(mesh.check(),1e-12)
        inf=state(1.225,0,0,101325);U=np.broadcast_to(inf,(60,16,4)).copy()
        rate,_=rhs(U,mesh,inf,1.4)
        self.assertLess(np.max(abs(rate)),1e-7)

    def test_outer_nose_curvature_matches_thermal_input(self):
        distance=1e-5
        rr=float(radius(distance,self.c))
        inferred_radius=rr*rr/(2*distance)
        self.assertLess(abs(inferred_radius/self.c['geometry']['nose_radius_m']-1),1e-4)

    def test_equal_state_flux_is_physical_flux(self):
        U=state(1,200,0,100000);S=np.array([2.,0.]);F,_=flux(U,U,S,1.4)
        expected=np.array([400,2*(200**2+100000),0,2*(U[3]+100000)*200])
        np.testing.assert_allclose(F,expected,rtol=1e-13)

    def test_mass_ratios_include_reserve(self):
        prop=np.array([400000,140000,25000.]);r=ledger(self.c,prop)
        for row,s in zip(r,self.c['stages']):
            expected=G0*s['isp_sizing_s']*np.log((row['upper_stack_kg']+row['dry_kg']+row['propellant_kg'])/(row['upper_stack_kg']+row['dry_kg']+row['reserve_kg']))
            self.assertAlmostEqual(row['dv_ideal_mps'],expected,places=9)
        self.assertAlmostEqual(r[0]['cutoff_mass_kg']-r[0]['dry_kg']-r[0]['reserve_kg'],r[1]['ignition_mass_kg'])

    def test_zero_heat_ambient_equilibrium(self):
        c=copy.deepcopy(self.c);c['thermal']['initial_K']=300;c['thermal']['ambient_K']=300
        back,_,T=conduction(c,np.linspace(0,1000,100),np.zeros(100),.03)
        np.testing.assert_allclose(T,300,atol=1e-8)

    def test_vent_no_reverse_flow_and_choking(self):
        self.assertEqual(orifice_flux(100000,110000,300),0)
        choked=orifice_flux(200000,0,300)
        self.assertAlmostEqual(choked,orifice_flux(200000,50000,300),places=10)
        self.assertGreater(choked,orifice_flux(200000,180000,300))


if __name__=='__main__':unittest.main()
