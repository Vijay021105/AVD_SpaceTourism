# TARA_LV computed mission report

This run is a conceptual analysis of a 12,500 kg crew spacecraft launched to the BAS orbital plane, waiting two nominal 400 km periods, then executing a J2-targeted phasing transfer and a prescribed terminal docking approach. Descent and re-entry are outside this run. **Readiness: CONCEPTUAL_RESEARCH_NOT_HUMAN_RATED.**

## Computed baseline

| Quantity | Result |
|---|---:|
| Liftoff mass | 728.535 t |
| Overall conceptual length | 93.505 m |
| Maximum exterior diameter | 4.000 m |
| Launcher sizing ideal delta-v | 11100.00 m/s |
| Spacecraft allocation | 300.00 m/s |
| Injection elapsed time | 1008.658 s |
| Injection altitude | 400.0000 km |
| Injection inertial speed | 7668.5582 m/s |
| Injection radial speed | -0.0000 m/s |
| Injection inclination | 51.600000 deg |
| Osculating perigee / apogee | 400.00 / 400.00 km |
| Peak dynamic pressure | 32.904 kPa |
| Maximum proper acceleration | 3.501 g |
| Doubled-resolution position difference | 479.412 m |
| Doubled-resolution velocity difference | 0.7092 m/s |
| Nominal orbital period | 92.5604 min |
| Two-period coast | 3.0853 h |
| BAS phase at departure | 112.4457 deg |
| Phasing duration | 30.7573 h |
| Spacecraft maneuver delta-v | 161.898 m/s |
| Spacecraft delta-v reserve | 138.102 m/s |
| Transfer endpoint miss at 1 km standoff | 0.0000 m |
| Minimum phasing altitude | 242.096 km |
| Prescribed docking elapsed time | 34.789 h |
| Ascent stagnation heat-flux peak | 60.321 kW/m2 |
| Ascent TPS thickness | 5.124 mm |

The 11,400 m/s number is a **sizing budget**, split into launcher ideal capability and spacecraft allocation. It is not orbital speed and is not added to the gravity and drag losses a second time. Actual integrated thrust/mass, gravity projection and drag histories are exported separately. Upper-stage unused propellant is retained on that stage; it is not spacecraft docking propellant.

The BAS initial state is fixed at the common liftoff epoch. The achieved ascent state initializes the crew coast. The transfer is shot to the propagated BAS standoff state under J2. The terminal distance reaches zero because a smooth approach path is prescribed; that endpoint is **not evidence of autonomous docking or contact qualification**.

## Stage ledger

| Stage | Engine | Count | Propellant kg | Dry kg | Ignition T/W | Ideal delta-v m/s |
|---|---|---:|---:|---:|---:|---:|
| 1 | M1D_standard | 10 | 492174.44 | 38597.20 | 1.183 | 3196.92 |
| 2 | M1D_standard | 2 | 139565.54 | 9308.44 | 0.961 | 3658.41 |
| 3 | CE32_project_CE20_proxy | 1 | 30814.69 | 5574.20 | 0.450 | 4244.67 |

## Acceptance and open work

| Check | Status | Computed value | Criterion |
|---|---|---|---|
| constrained_sizing | **PASS** | 11100.00000000606 | launcher ideal delta-v >= 11100.0 m/s and selected ignition T/W |
| dry_hardware_floor | **PASS** | 0.8006985487227212 | estimated hardware <= dry budgets; missing hardware still requires design |
| tank_capacity | **PASS** | 0.0 | mass-volume equality <0.01 kg |
| diameter | **PASS** | 3.958 | installed tanks inside 4 m envelope |
| engine_installation | **OPEN** | 4.748593038111737 | 4 m radial and pairwise gimbal clearance; assumed supplier dimensions |
| orbit_injection | **PASS** | 399999.99999988545 | altitude ±2 km, speed ±2 m/s, radial speed <2 m/s, inclination ±0.1 deg |
| ascent_integration_refinement | **PASS** | 479.4115440928071 | doubled fraction resolution: <1000 m and <2 m/s |
| dynamic_pressure | **PASS** | 32903.64892826528 | <= 45000.0 Pa |
| q_alpha | **PASS** | 3000.000000000082 | <= 3000.0 Pa rad |
| crew_acceleration | **PASS** | 3.5006116704605668 | <= 3.5 g; physiological qualification not established |
| phasing_shooting | **PASS** | 2.184782717992487e-06 | <10 m endpoint miss before velocity-match impulse |
| phasing_orbit_clearance | **PASS** | 242096.32882029284 | screening perigee >180 km |
| spacecraft_delta_v | **PASS** | 161.89793534254235 | within allocated 300 m/s; no double counting launcher ideal reserve |
| spacecraft_propellant | **PASS** | 628.5300888816084 | required propellant <= included payload propellant |
| approach_thrust | **PASS** | 17.646708341798625 | prescribed-path demand <= assumed available thrust |
| structural_screens | **OPEN** | -0.6775062093333608 | all calculated metallic and buckling margins >=0; qualification work remains |
| TPS_ascent | **PASS** | 422.9952987887043 | bondline <= assumed 423 K for ascent only |
| cryo_boiloff_reserve | **OPEN** | 0.02656666678495702 | pad+ascent boiloff <= selected 1% liquid reserve; operational loading must account for losses |
| POGO_assumed_loop | **OPEN** | RHP_INSTABILITY | no unstable roots in supplied scalar delay loop; not flight validation |
| Euler_iterations | **PASS** | 6 | 6 admissible cases must meet residual + force criteria |
| human_rating | **OPEN** | NOT_ESTABLISHED | hazard analysis, abort/escape, FDIR, ECLSS, qualification, fracture control, contact dynamics and test evidence required |
| viscous_separated_CFD | **OPEN** | NOT_IMPLEMENTED | 3-D transition/RANS/URANS, base/plume and incidence validation required |
| aeroelastic_stability | **OPEN** | NOT_ESTABLISHED | two-way aerodynamic/structural/control coupling and flutter/divergence analysis required |

## CFD evidence

The trajectory uses the engineering database. The Euler snapshots are an independent diagnostic and are **not fed back as a validated drag table**. Slip-wall Euler resolves pressure-wave structure but cannot resolve boundary layers, viscous separation, base recirculation, turbulence, aeroacoustic generation or heat transfer. Axisymmetry at zero incidence cannot identify Cn, Cnalpha or windward/leeward pressure differences. Those outputs use separately labelled engineering estimates.

| Case | Status | Reason |
|---|---|---|
| M00.30_S1 | CONVERGED |  |
| M00.80_S1 | CONVERGED |  |
| M00.98_S1 | CONVERGED |  |
| M01.20_S1 | CONVERGED |  |
| M02.00_S1 | CONVERGED |  |
| M03.00_S1 | CONVERGED |  |
| M05.00_S1 | EXCLUDED | perfect_gas_temperature_limit |
| M08.00_S3 | EXCLUDED | perfect_gas_temperature_limit;rarefied_flow;upper_atmosphere_extrapolated |
| M15.00_S3 | EXCLUDED | Mach_not_reached;perfect_gas_temperature_limit;rarefied_flow;upper_atmosphere_extrapolated |
| M25.00_S3 | EXCLUDED | Mach_not_reached;perfect_gas_temperature_limit;rarefied_flow;upper_atmosphere_extrapolated |

## Book-to-implementation map

The supplied Suresh and Sivan volume is the organizing reference. Its chapter numbering, rather than the unrelated Venus-mission README numbering, is used here. Numerical assumptions are not attributed to the book unless it actually supplies them.

| Topic | Book chapter | Implementation | Scope |
|---|---|---|---|
| Mission definition and coupling | 2, 7, 8 | inputs/mission.csv; tara/pipeline.py | Implemented dependency chain; open gates retained |
| Orbit, J2 and injection | 3, 4, 6 | tara/orbits.py; tara/rendezvous.py | Numerical ECI J2; spherical surface and GMST |
| Staging and mass allocation | 5 | tara/sizing.py | Constrained numerical optimization plus hardware dry floors |
| Propulsion | 9 | inputs/engines.csv; tara/ascent.py | Existing-engine performance proxy; no engine internal design |
| Aerodynamics | 10 | tara/aerodynamics.py; tara/cfd.py | Engineering database plus inviscid axisymmetric Euler fields |
| Structures and materials | 11 | tara/structures.py; tara/beam.py | Section cuts, metallic screening, beam modes; detailed shell FEM open |
| Thermal | 12 | tara/thermal.py | Ascent heating, conduction TPS, heat leak, ideal gas vent sizing |
| Stage auxiliary systems | 13 | tara/dynamics.py | Separation impulse and tip-off screens |
| Guidance and redundancy | 14 | tara/ascent.py; inputs/human_rating.csv | Point-mass targeting and illustrative channel reliability; full GNC open |
| Re-entry | 15 | Not included | Excluded by requested scope |

## Figures

- [01_mission_overview.png](figures/01_mission_overview.png)
- [02_mass_and_staging.png](figures/02_mass_and_staging.png)
- [03_dimensioned_vehicle.png](figures/03_dimensioned_vehicle.png)
- [04_engine_installation.png](figures/04_engine_installation.png)
- [05_ascent_time_histories.png](figures/05_ascent_time_histories.png)
- [06_propulsion_and_velocity_ledger.png](figures/06_propulsion_and_velocity_ledger.png)
- [07_orbit_and_phasing.png](figures/07_orbit_and_phasing.png)
- [08_ground_tracks.png](figures/08_ground_tracks.png)
- [09_rendezvous_history.png](figures/09_rendezvous_history.png)
- [10_docking_feasibility.png](figures/10_docking_feasibility.png)
- [11_aerodynamic_mach_sweep.png](figures/11_aerodynamic_mach_sweep.png)
- [12_aero_polars_and_viscous_estimates.png](figures/12_aero_polars_and_viscous_estimates.png)
- [13_aerodynamic_loads.png](figures/13_aerodynamic_loads.png)
- [14_ground_wind_loads.png](figures/14_ground_wind_loads.png)
- [15_section_load_envelopes.png](figures/15_section_load_envelopes.png)
- [16_structural_margins_and_mass.png](figures/16_structural_margins_and_mass.png)
- [17_vehicle_mode_shapes.png](figures/17_vehicle_mode_shapes.png)
- [18_acoustics_buffet_and_slosh.png](figures/18_acoustics_buffet_and_slosh.png)
- [19_pogo_feedback.png](figures/19_pogo_feedback.png)
- [20_ascent_thermal_and_TPS.png](figures/20_ascent_thermal_and_TPS.png)
- [21_cryo_insulation_and_vents.png](figures/21_cryo_insulation_and_vents.png)
- [22_stage_separation.png](figures/22_stage_separation.png)
- [23_computed_Euler_fields.png](figures/23_computed_Euler_fields.png)
- [24_FVM_grid.png](figures/24_FVM_grid.png)
- [25_CFD_convergence.png](figures/25_CFD_convergence.png)
- [26_surface_pressure_and_incidence.png](figures/26_surface_pressure_and_incidence.png)
- [27_CFD_grid_sensitivity.png](figures/27_CFD_grid_sensitivity.png)
- [28_structural_elements.png](figures/28_structural_elements.png)
- [29_design_iteration_histories.png](figures/29_design_iteration_histories.png)
- [30_redundancy_and_engine_thrust.png](figures/30_redundancy_and_engine_thrust.png)
- [31_shock_and_sonic_diagnostics.png](figures/31_shock_and_sonic_diagnostics.png)
- [32_pogo_oscillation_history.png](figures/32_pogo_oscillation_history.png)

## Reproducibility

All 14 subsystem CSV inputs are hashed in `mission_summary.json`. `effective_configuration.json` contains the complete generated model. `data/` contains numeric histories, iteration logs, geometry, margins and open gates. `fields/` contains actual finite-volume arrays. `gifs/` and `videos/` are display products. GIF and MP4 field sweeps change operating condition and do not represent physical unsteady flow.

See `../docs/WORK_PLAN.md`, `../docs/METHODS.md`, `../docs/DATA_DICTIONARY.md` and `../docs/SOURCES.md` before changing model fidelity or treating a screen as a design release.
