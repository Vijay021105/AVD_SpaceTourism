# Data dictionary

All input CSVs are required. Unknown or missing CSVs are rejected so a second hidden input table cannot silently become authoritative. Scalar tables have `parameter,value,unit,source` columns. The runtime uses SI values; unit strings are documentation and are not a unit conversion parser. Do not change them without changing numeric units.

| Input table | Meaning |
|---|---|
| `mission.csv` | Payload, diameter, orbit, epoch, launch location, total/spacecraft delta-v allocation, T/W, q/q-alpha/g limits, phasing and approach settings |
| `stages.csv` | Stage order, engine reference/count, assumed structural factors/bounds, fixed dry allowance, residual reserve, mass seeds, sizing Isp and bay/skirt lengths |
| `engines.csv` | One standard M1D record shared by S1/S2; one explicitly identified CE20-based cryogenic proxy; performance, mass, throttle and envelope assumptions |
| `fluids.csv` | Density, nominal liquid temperature, absolute pressure and latent heat |
| `payload.csv` | Component masses summing exactly to 12,500 kg, including spacecraft maneuver propellant |
| `geometry.csv` | Envelope allowances, insulation, ullage, dome ratio, tank gaps, capsule/trunk envelope, nozzle clearance |
| `aerodynamics.csv` | Correlation constants, transition Reynolds assumption, reference incidence, screening Mach bound and Euler convergence thresholds |
| `materials.csv` | Assumed welded-metal and TPS constant properties; not a qualified materials handbook |
| `structures.csv` | Separate safety factors, knockdown, load uncertainty, skin, baffles, joints and mass-estimation assumptions |
| `thermal.csv` | Heating, radiation, TPS limits, insulation, hold time, plume exposure and vent assumptions |
| `dynamics.csv` | Damping, assumed acoustic/buffet spectra, scalar POGO loop, ground wind and separation assumptions |
| `guidance.csv` | Optimization seed, vertical hold, pitch interpolation, coast times, S3 cutoff and q-relief settings |
| `wind.csv` | Strictly increasing geometric altitude and assumed east/north wind; not weather |
| `human_rating.csv` | Channel counts, minimum healthy count, illustrative failure/common-cause values and allowance masses |

The 14 files are hashed before the run. Input hash equality is required for `--reuse-ascent`. A changed source file requires a fresh run.

## Generated subsystem CSVs

| Output | Important fields / interpretation |
|---|---|
| `stages.csv` | Propellant/dry/upper stack, usable/reserve mass, ignition/burnout mass, epsilon budget/effective, ideal delta-v, T/W, full flow and minimum burn time |
| `mass_geometry_coupling.csv` | Each sizing ↔ hardware-floor iteration and remaining estimated mass deficit |
| `structural_fraction_trade.csv` | Conditional factor sensitivity; feasibility column is mandatory when comparing masses |
| `tanks.csv` | `z0_m,z1_m`, inner radius, head depth, barrel, wall, insulation, area, liquid/total volume, capacity and loaded mass |
| `components.csv` | Axial boundaries of engine bays, tanks, skirts, trunk and capsule |
| `outer_profile.csv` / `vehicle_profile.dxf` | Meridian in metres; z measured aft to nose; DXF units set to metres |
| `engine_layout.csv` | Assumed engine-center coordinates, nozzle radius and gimbal swept radius; installation check is separate |
| `ascent.csv` | ECI r/v, attached mass, stage, time, altitude, inertial and air-relative speed, wind, Mach, q, incidence, coefficients, thrust/drag/normal force, pitch/rate, throttle, proper g and velocity integrals |
| `events.csv` | Ignition/cutoff/separation times and pre/post attached masses, with event ECI state |
| `optimized_guidance.csv` | Final command parameters used for the published ascent |
| `ascent_targeting_iterations.csv` | Actual coarse targeting residual evaluations |
| `injection_state.csv` | Achieved ECI state and osculating elements, with time since the common UTC epoch |
| `BAS_ephemeris.csv` | Independently propagated station ECI state and osculating elements |
| `rendezvous.csv` | Crew and BAS ECI r/v, RTN relative components, range, rotating relative speed, ground track and mission phase |
| `rendezvous_burns.csv` | Both impulsive velocity changes; terminal low-thrust cost appears in summary |
| `rendezvous_targeting_iterations.csv` | Transfer shooting evaluations and actual endpoint miss |
| `aerodynamic_database.csv` | Stage-specific engineering Ca/Cd/Cn/Cnalpha/Xcp/roll at reference incidence and nominal reference Reynolds evaluation |
| `aero_polars.csv` | Engineering incidence sweep; Cnalpha column is linear slope per radian |
| `cfd_cases.csv` | Actual state, model applicability, mesh, residual/force criteria, converged/not-converged/excluded status |
| `cfd/*_residual.csv` | Actual Euler iteration, residual, force, mass-flux imbalance and pseudo-time diagnostic |
| `cfd/*_surface.csv` | Zero-incidence computed wall pressure and Cp, including nonconverged diagnostic iterates when so labelled |
| `cfd_grid_sensitivity.csv` | Coefficient changes and convergence status for the grid sequence |
| `windward_leeward_Cp_proxy.csv` | Separate analytical incidence exposure estimate; never labelled CFD |
| `loads_by_case.csv` | Simultaneous axial/shear/bending versus z at selected flight times |
| `loads_envelope.csv` | Componentwise maxima over all ascent samples; maxima need not occur simultaneously |
| `structural_margins.csv` | Metallic yield/ultimate, compression buckling and joint bolt margins; blank means inactive/not calculated, not zero |
| `hardware_mass_breakdown.csv` | Estimated physical item masses, including engines and declared allowances |
| `dry_mass_closure.csv` | Estimated hardware versus assigned dry budget; remaining budget is not proof of detailed completeness |
| `modal_frequencies.csv` / `mode_shapes.csv` | Wet/dry free-free and ground frequencies; residuals and mode classification; normalized shapes |
| `acoustic_buffet_response.csv` | Assumed pressure PSD and one-mode acceleration PSD; dimensional frequency in Hz |
| `pogo_transient.csv` | Test-disturbance response of the same assumed delayed loop, truncated at the linear-range screen |
| `pogo_frequency_response.csv` / `pogo_assessment.json` | Assumed-loop Bode/Nyquist, exact-delay root and stability status |
| `slosh_modes.csv` | Cylindrical first lateral slosh mode at 1 g and assumed damping/baffles |
| `ground_wind_loads.csv` | Wind, line load, cantilever shear, bending and deflection versus height |
| `separation_clearance.csv` | Cold-separation two-body axial clearance and assumed tip-off sweep |
| `redundancy_assessment.csv` | Illustrative function-failure arithmetic; not loss-of-crew probability |
| `thermal_history.csv` | Stagnation correlation, cumulative exposure, radiative equilibrium and assumed base heating |
| `TPS_history.csv` / `TPS_sizing_iterations.csv` | Actual conduction response and thickness search |
| `capsule_heat_distribution.csv` | Assumed spatial exposure envelope and selected ascent TPS thickness |
| `cryo_heat_leak_and_vents.csv` | Conductive leak, pad+flight boiloff, reserve check and ideal-gas vent area |
| `interstage_vents.csv` / `vent_sizing_history.csv` | Quasi-steady bay vent requirement/capacity |
| `design_checks.csv` | PASS or OPEN, value, criterion and evidence scope |
| `book_traceability.csv` | Topic/chapter/implementation mapping |

## Frames, signs and missing values

`z_m` increases from original vehicle aft plane to nose. CFD `x_from_nose_m` increases aft. `Xcp_from_nose_m` uses the CFD nose reference; the loads module converts it to vehicle z.

Rendezvous RTN is based on BAS: R outward, T along orbital motion, N along orbital angular momentum. Along-track approach begins behind BAS. Relative velocity subtracts frame rotation, not only inertial vector difference.

Axial envelope is a positive magnitude; selected axial case force uses the compression-demand convention. Shear and moment in simultaneous cases preserve sign. Coefficient force directions are described in METHODS.md. Blank buckling margins identify no compressive branch; unsupported analyses have explicit OPEN status rather than zeros.

Generated JSON forbids NaN. NPZ files retain computed finite-volume arrays with their true grid, free-stream Mach and conservative state. All array fields have SI dimensions except Mach and the plotted normalized velocities/Cp.
