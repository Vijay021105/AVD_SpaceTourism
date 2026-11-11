"""Write the reviewed nominal input deck. Run explicitly; never invoked by the driver.
Existing CSV edits are overwritten only when --overwrite is supplied.
"""
from pathlib import Path
import argparse
from tara.io import write_csv


def regenerate(folder):
    def scalars(name, items):
        write_csv(folder/f'{name}.csv', [dict(parameter=k,value=v,unit=u,source=s) for k,v,u,s in items])
    A='Conceptual design assumption; replace with qualified data'
    U='User requirement 2026-11-11'
    B='Suresh & Sivan 2015; method not vehicle-specific numerical value'
    scalars('mission',[
        ('payload_kg',12500,'kg',U),('diameter_m',4,'m',U),('orbit_altitude_m',400000,'m',U),
        ('inclination_deg',51.6,'deg',U),('total_dv_budget_mps',11400,'m/s',U),
        ('spacecraft_dv_allocation_mps',300,'m/s',A),('coast_periods',2,'period',U),
        ('phasing_revolutions',20,'count',A),('phasing_fraction',.25,'turn','Non-integer arrival geometry avoids near-return shooting singularity'),('bas_u0_deg',256,'deg','Inherited orbit study phase; common liftoff epoch'),
        ('epoch_utc','2027-01-01T00:00:00Z','UTC',A),('site_lat_deg',8.3943,'deg','Inherited Kulasekarapattinam vicinity; no range approval'),
        ('site_lon_deg',78.0516,'deg','Inherited conceptual site'),('site_height_m',26.4,'m',A),
        ('q_limit_Pa',45000,'Pa',A),('qalpha_limit_Pa_rad',3000,'Pa rad',A),
        ('proper_g_limit',3.5,'g',A),('liftoff_tw_min',1.15,'ratio',A),('stage2_tw_min',.65,'ratio',A),
        ('stage3_tw_min',.45,'ratio',A),('standoff_m',1000,'m',A),('approach_duration_s',2400,'s',A),
        ('spacecraft_isp_s',320,'s',A),('approach_thrust_available_N',400,'N',A),
        ('orbit_alt_tolerance_m',2000,'m',A),('orbit_speed_tolerance_mps',2,'m/s',A),
        ('orbit_incl_tolerance_deg',.1,'deg',A)])
    write_csv(folder/'engines.csv',[
        dict(engine_id='M1D_standard',thrust_sl_N=845000,thrust_vac_N=845000*311/282,isp_sl_s=282,isp_vac_s=311,
             min_throttle=.57,engine_mass_kg=470,OF=2.35,fuel='RP1',oxidizer='LOX',nozzle_exit_diameter_m=1.05,
             length_m=2.9,gimbal_deg=5,source='SpaceX guide 2025 for 845 kN and throttle; Isp mass O/F dimensions provisional inherited assumptions'),
        dict(engine_id='CE32_project_CE20_proxy',thrust_sl_N=0,thrust_vac_N=215746.3,isp_sl_s=0,isp_vac_s=442.5,
             min_throttle=1,engine_mass_kg=588,OF=5.05,fuel='LH2',oxidizer='LOX',nozzle_exit_diameter_m=1.5,
             length_m=3.2,gimbal_deg=4,source='ISRO CE20 19-22 tf and 442.5 s; 22 tf selected; other fields provisional')])
    write_csv(folder/'stages.csv',[dict(stage=j+1,engine_id=('M1D_standard' if j<2 else 'CE32_project_CE20_proxy'),engine_count=[10,2,1][j],
          epsilon=[.065,.06,.11][j],epsilon_low=[.06,.05,.09][j],epsilon_high=[.09,.09,.16][j],
          fixed_dry_allowance_kg=[800,400,250][j],reserve_fraction=.01,initial_prop_kg=[340000,110000,24000][j],
          isp_sizing_s=[295,311,442.5][j],engine_bay_m=[3.5,3.4,3.7][j],interstage_m=[2.5,2.3,.8][j],
          source='Selected optimistic structural factors; bounded sensitivity; Ch5 methodology. Not measured human-rated fractions') for j in range(3)])
    write_csv(folder/'fluids.csv',[
        dict(fluid='RP1',density_kgm3=810,temperature_K=293,pressure_abs_Pa=250000,latent_Jkg=0),
        dict(fluid='LOX',density_kgm3=1141,temperature_K=90,pressure_abs_Pa=300000,latent_Jkg=213000),
        dict(fluid='LH2',density_kgm3=70.85,temperature_K=20.3,pressure_abs_Pa=250000,latent_Jkg=446000)])
    scalars('geometry',[
        ('skin_m',.006,'m',A),('radial_clearance_m',.015,'m',A),('tank_gap_m',.35,'m',A),
        ('ullage_fraction',.03,'fraction',A),('dome_depth_over_radius',.5,'ratio','2:1 ellipsoidal heads'),
        ('cryo_insulation_m',.05,'m',A),('rp1_insulation_m',.02,'m',A),
        ('capsule_length_m',4.4,'m','Concept capsule envelope; not Dragon production geometry'),
        ('trunk_length_m',3.7,'m',A),('capsule_base_radius_m',1.9,'m',A),
        ('nose_radius_m',.25,'m',A),('engine_clearance_m',.025,'m',A)])
    scalars('aerodynamics',[
        ('base_cd',.12,'coefficient',A),('form_cd',.08,'coefficient',A),('wave_peak_cd',.45,'coefficient',A),
        ('cnalpha_sub_per_rad',2,'1/rad','Slender-body forebody approximation'),('transition_Re',5e5,'Re',A),
        ('alpha_reference_deg',3,'deg',A),('uncertainty_fraction',.25,'fraction',A),
        ('mach_max',30,'Mach','Database screening extends beyond perfect-gas CFD validity'),
        ('cfd_max_total_temp_K',1500,'K','Calorically perfect-gas model applicability screen'),
        ('cfd_residual_tolerance',1e-5,'relative',A),('cfd_drag_change_tolerance',.002,'relative',A)])
    write_csv(folder/'materials.csv',[
        dict(material='Al2219_weld_proxy',E_Pa=73e9,nu=.33,density_kgm3=2840,yield_Pa=150e6,ultimate_Pa=250e6,k_WmK=120,cp_JkgK=900,source=A),
        dict(material='TPS_silica_proxy',E_Pa=1e7,nu=.2,density_kgm3=150,yield_Pa=1e6,ultimate_Pa=1.5e6,k_WmK=.08,cp_JkgK=1000,source=A)])
    scalars('structures',[
        ('yield_fos',1.25,'factor','Selected metallic protoflight screen; NASA-STD-5001B Table 1 applicability must be reviewed'),
        ('ultimate_fos',1.4,'factor','Selected metallic screen; not complete human-rating qualification'),
        ('pressure_fos',1.5,'factor',A),('buckling_knockdown',.3,'factor',A),('load_uncertainty',1.15,'factor',A),
        ('shell_min_m',.003,'m',A),('effective_skin_m',.006,'m',A),('stiffener_mass_fraction',.35,'fraction',A),
        ('engine_frame_mass_per_N',.00012,'kg/N',A),('joint_mass_per_m',20,'kg/m circumference',A),
        ('bolt_allowable_N',45000,'N/bolt',A),('joint_bolts',96,'count',A),('joint_lever_arm_m',1.8,'m',A),
        ('baffle_count',3,'per tank',A),('baffle_width_m',.15,'m',A),('baffle_thickness_m',.003,'m',A),
        ('isogrid_pitch_m',.3,'m',A),('isogrid_rib_height_m',.025,'m',A),('isogrid_rib_width_m',.004,'m',A)])
    scalars('thermal',[
        ('sutton_graves_k',1.7415e-4,'SI','Sutton-Graves earth-air correlation; screening'),
        ('emissivity',.85,'fraction',A),('initial_K',293,'K',A),('bondline_limit_K',423,'K',A),
        ('tps_min_m',.005,'m',A),('tps_max_m',.2,'m',A),('cryo_foam_k_WmK',.025,'W/m/K',A),
        ('cryo_foam_density_kgm3',45,'kg/m3',A),('pad_hold_s',3600,'s',A),('ambient_K',300,'K',A),
        ('base_heat_fraction',.003,'fraction of chamber chemical power',A),('combustion_energy_Jkg',12e6,'J/kg mixture',A),
        ('vent_Cd',.65,'coefficient',A),('vent_diff_pressure_Pa',3000,'Pa',A),('bay_temperature_K',300,'K',A),
        ('LH2_vent_R_JkgK',4124,'J/kg/K',A),('LOX_vent_R_JkgK',259.8,'J/kg/K',A)])
    scalars('dynamics',[
        ('modal_damping',.02,'ratio',A),('OASPL_dB',145,'dB re 20uPa',A),('buffet_rms_q_fraction',.005,'ratio',A),
        ('pogo_feed_Hz',8,'Hz',A),('pogo_struct_Hz',10,'Hz',A),('pogo_feed_damping',.04,'ratio',A),
        ('pogo_struct_damping',.02,'ratio',A),('pogo_loop_gain',.15,'dimensionless',A),
        ('pogo_delay_s',.03,'s',A),('ground_wind_mps',20,'m/s',A),('gust_factor',1.5,'ratio',A),
        ('separation_push_N',20000,'N',A),('separation_push_duration_s',.15,'s',A),('separation_clearance_m',2,'m',A),
        ('separation_tipoff_degps',.5,'deg/s',A)])
    scalars('guidance',[
        ('vertical_time_s',12,'s',A),('pitch_fraction_1',.3,'fraction',A),('pitch_fraction_2',.65,'fraction',A),
        ('pitch1_deg',75.820312773394,'deg above horizon',A),('pitch2_deg',63.137495236892,'deg above horizon',A),('pitch3_deg',39.532002192597,'deg above horizon',A),
        ('stage2_end_deg',27.844135309831,'deg above horizon',A),
        ('stage3_start_deg',32.362898366654,'deg above horizon',A),('stage3_end_deg',10.885174407469,'deg above horizon',A),
        ('coast12_s',3,'s',A),('coast23_s',134.90339451515,'s',A),('plane_trim_deg',0.10560784273312,'deg','Numerical yaw-plane trim; BAS orbit remains fixed'),('stage3_use_fraction',0.78977012939938,'fraction',A),
        ('q_throttle_on_fraction',.75,'fraction',A),('q_controller_gain',.65,'gain',A)])
    write_csv(folder/'wind.csv',[dict(altitude_m=h,east_mps=e,north_mps=n,source='Assumed design wind; not launch weather') for h,e,n in
        [(0,3,1),(1000,5,2),(3000,10,3),(6000,16,5),(10000,25,8),(15000,30,5),(20000,15,2),(40000,5,0),(80000,0,0),(150000,0,0)]])
    write_csv(folder/'payload.csv',[dict(component=k,mass_kg=m,source=A) for k,m in [
        ('crew_cabin_structure_and_TPS',4900),('crew_seats_life_support',2200),('trunk_power_avionics',1700),
        ('crew_and_cargo',800),('rendezvous_propulsion_dry',700),('spacecraft_propellant',1200),('adapter_and_margin',1000)]])
    write_csv(folder/'human_rating.csv',[
        dict(function='flight_computers',channels=3,required_healthy=2,unit_failure_probability=.001,common_cause_beta=.05,mass_allowance_kg=120,status='illustrative independent-channel screen'),
        dict(function='IMU_navigation',channels=3,required_healthy=2,unit_failure_probability=.002,common_cause_beta=.05,mass_allowance_kg=60,status='illustrative independent-channel screen'),
        dict(function='power_buses',channels=2,required_healthy=1,unit_failure_probability=.001,common_cause_beta=.1,mass_allowance_kg=200,status='illustrative independent-channel screen'),
        dict(function='separation_initiation',channels=2,required_healthy=1,unit_failure_probability=.005,common_cause_beta=.1,mass_allowance_kg=60,status='illustrative independent-channel screen')])


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,default=Path(__file__).parent/'inputs');p.add_argument('--overwrite',action='store_true');a=p.parse_args()
    a.inputs.mkdir(parents=True,exist_ok=True)
    if list(a.inputs.glob('*.csv')) and not a.overwrite:p.error('Existing deck: use --overwrite explicitly')
    regenerate(a.inputs)
    print('Nominal CSV deck written:',a.inputs)
