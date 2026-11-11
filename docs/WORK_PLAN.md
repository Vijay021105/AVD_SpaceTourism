# Work plan and acceptance evidence

This sequence follows the integrated-design approach of B.N. Suresh and K. Sivan, *Integrated Design for Space Transportation System* (2015). The book's actual chapters are used here: orbit 3–4, staging 5, environment 6, mission requirements 7, integration 8, propulsion 9, aerodynamics 10, structures 11, thermal 12, auxiliary systems 13 and guidance/control 14. Chapter 15 re-entry is excluded.

## Phase 1 — mission requirements and configuration control

- [x] Define payload, maximum diameter, orbit, waiting sequence, phasing and ascent-only scope.
- [x] Freeze the 10 standard M1D / 2 standard M1D / 1 cryogenic architecture.
- [x] Create one SI subsystem CSV deck with source and assumption annotations.
- [x] Use a single epoch and a single achieved injection state.
- [x] Separate the 11,400 m/s sizing budget from actual delivered speed and loss histories.
- [ ] Confirm launch site/range, launch window, BAS ephemeris and approved operational orbit.
- [ ] Replace the inherited BAS phase with the agreed station state at the chosen epoch.
- [ ] Confirm payload manifest, adapter, life support, escape hardware and propellant ownership.

Acceptance evidence: `inputs/mission.csv`, `payload.csv`, `engines.csv`, `mission_summary.json` input hashes and `data/injection_state.csv`.

## Phase 2 — astrodynamics and mission geometry

- [x] Propagate Cartesian ECI states under central gravity and J2.
- [x] Derive nominal period, launch-plane geometry, launch azimuth, phase at departure and injection vector.
- [x] Propagate the achieved crew orbit for two nominal periods without resetting it.
- [x] Solve a multi-revolution transfer to the J2-propagated station standoff state.
- [x] Account for both impulses, terminal acceleration demand and spacecraft propellant.
- [x] Plot orbit geometry, ground tracks, relative range, RTN components and approach history.
- [ ] Replace the prescribed approach with relative-navigation, RCS, sensor and docking-contact models.
- [ ] Add navigation/ephemeris dispersions, conjunction constraints, lighting and communications visibility.

Acceptance evidence: injection accuracy, two-body exact solution, J2 node-drift check, endpoint miss, maneuver budget and minimum phasing altitude. Ground tracks are not range clearance. Ballistic impact/descent and debris dispersion are not computed.

## Phase 3 — staging, propulsion and dimensions

- [x] Optimize stage propellant mass at selected structural factors and ignition T/W constraints.
- [x] Apply reserve consistently at burnout and separation.
- [x] Regenerate fuel/oxidizer volumes with density, O/F, ullage and ellipsoidal heads.
- [x] Include insulation in the 4 m exterior envelope.
- [x] Feed estimated shells, insulation, engines, frames, skirts, joints and baffles into dry-mass floors.
- [x] Export exact conceptual tank/component dimensions, meridian DXF and engine-layout coordinates.
- [x] Check static nozzle packing and assumed gimbal swept envelopes.
- [ ] Obtain standard Merlin vacuum-operation, ignition/restart, throttle, mass, envelope and interface data.
- [ ] Replace the CE20 performance proxy with the actual selected cryogenic-stage/engine specification.
- [ ] Resolve any engine-installation gate before freezing the 4 m engine bay.
- [ ] Size feed lines, pumps, pressurant bottles, valves, residuals, purge, ignition and engine isolation.
- [ ] Include supplier-qualified structural mass and component growth margins.

Acceptance evidence: constrained optimizer status, exact tank capacity, dry item budgets, T/W, engine clearance and supplier review. Structural-factor sensitivity is conditional on the assumed factors, not a literature-proven optimum.

## Phase 4 — ascent trajectory and guidance

- [x] Integrate 3-D point-mass ECI motion with J2, atmospheric rotation, wind and stage events.
- [x] Model pressure-dependent thrust, minimum throttle, q relief, q-alpha relief and proper-acceleration limiting.
- [x] Numerically target altitude, radial speed, tangential speed and inclination.
- [x] Retain thrust, mass, pitch, pitch rate, q, coefficients, acceleration, wind and velocity integrals.
- [x] Compare a doubled integration resolution.
- [ ] Add TVC actuator limits, angular dynamics, rate limits, navigation and a closed-loop autopilot.
- [ ] Analyze propulsion failures, engine-out control authority and ascent abort corridors with approved assumptions.
- [ ] Run guidance robustness and Monte Carlo dispersions after nominal physical gates close.

Acceptance evidence: orbit/state error, perigee/apogee, q/q-alpha/g gates, event mass closure, numerical step refinement. A numerical target solution alone does not demonstrate actuator or crew acceptance.

## Phase 5 — aerodynamics and CFD

- [x] Export Ca, Cd, Cn, Cnalpha, Xcp and the symmetric-body rolling coefficient across the engineering Mach range.
- [x] Export incidence polars, friction/wave/base build-ups and trajectory loads.
- [x] Compute actual axisymmetric Euler u, v, p and Mach snapshots on a body-fitted ring grid.
- [x] Export raw fields, geometry, actual residuals and force convergence.
- [x] Run grid sensitivity and explicit continuum/temperature/model-bound checks.
- [x] Keep computed zero-incidence pressure separate from analytical incidence pressure proxies.
- [ ] Obtain laminar/transition and fully turbulent 3-D RANS solutions, with verified wall treatment and mesh convergence.
- [ ] Add incidence, protuberances, engine base/wake, separated-flow, plume interaction and real geometry.
- [ ] Use URANS/LES or experiments for transonic buffet, unsteady pressure spectra and aeroacoustic sources.
- [ ] Add reacting/equilibrium-gas or rarefied methods where perfect-gas continuum assumptions fail.
- [ ] Replace the engineering database with validated force/moment tables, then re-run ascent and structural coupling.

Acceptance evidence: iterative residual + force convergence, mesh/far-field sensitivity, force integration, wind tunnel/benchmark comparison, model applicability and uncertainty. A colorful field image is not this evidence.

## Phase 6 — structures, joints, vibroacoustics and aeroelasticity

- [x] Conserve the distributed attached mass and form flight axial/shear/bending envelopes.
- [x] Compute simultaneous selected load cases with inertial relief and TVC trim.
- [x] Screen tank pressure/membrane, equivalent metallic yield/ultimate, compression buckling and joint bolts.
- [x] Estimate baffles, isogrid/stiffener allowance, thrust frames, skirts and joint-ring masses.
- [x] Compute free-free wet/dry beam modes, mode shapes and fixed-base ground modes.
- [x] Compute linear response to assumed acoustic/buffet spectra, scalar delayed POGO feedback, slosh and ground-wind loads.
- [x] Export two-body separation clearance and assumed tip-off histories.
- [ ] Run shell/solid FEM of tanks, interstages, thrust frames, joints, cutouts and local reinforcements.
- [ ] Size ribs/stringers/isogrid and update mass/stiffness from real section properties.
- [ ] Evaluate weld allowables, fracture, fatigue, proof/burst, thermal stress and imperfection-sensitive buckling.
- [ ] Identify feed/engine/structure transfer functions for POGO and test accumulator mitigation.
- [ ] Couple measured slosh/baffle damping and liquid-acoustic modes.
- [ ] Add two-way aeroelastic and aeroservoelastic models for flutter, divergence, buffet and control stability.
- [ ] Validate ground-wind vortex shedding/lock-in, hold-down loads and acoustic environments.
- [ ] Validate separation with 6-DOF contact, plume, shock and clearance dispersions.

Acceptance evidence: positive margins at combined time-consistent design loads, mass closure, FEM mesh convergence, modal correlation, identified dynamic inputs and qualification tests.

## Phase 7 — thermal and fluid venting

- [x] Compute ascent stagnation heating and a labelled capsule exposure distribution.
- [x] Size a 1-D transient conductive TPS with surface radiation and bondline temperature constraint.
- [x] Export assumed plume/base heat exposure separately from aerodynamic heating.
- [x] Estimate cryogenic insulation heat leak, pad+ascent boiloff and reserve adequacy.
- [x] Size quasi-steady interstage depressurization and ideal-gas cryogen vents.
- [ ] Obtain CFD/experimental convective edge conditions and plume radiation/impingement environments.
- [ ] Model local interstage conduction/radiation and engine-bay shields.
- [ ] Add variable material properties, contact resistance, local hot spots and transient multilayer structures.
- [ ] Model MLI, foam aging, moisture, pressure-control logic, two-phase venting and propellant conditioning.
- [ ] Close TPS and insulation mass against the detailed stage/payload ledger.

Acceptance evidence: heat balance, transient/mesh refinement, allowable temperatures, boiloff/pressure-control capability and thermal qualification. Ascent TPS sizing does not cover re-entry.

## Phase 8 — human-flight assurance and integrated design closure

- [x] Expose factors, mass allowances, liquid reserves and simplified redundant-channel reliability with common cause.
- [x] Maintain pass/open gates and label evidence level in the report.
- [ ] Establish program-specific human-rating requirements and safety acceptance authority.
- [ ] Complete hazards/FMEA/FTA/common-cause analyses, FDIR, engine isolation, abort/escape and safe-haven assessment.
- [ ] Qualify ECLSS, crew accommodation, docking interfaces, software/hardware and ground operations.
- [ ] Complete test plans: engine life/acceptance, structural proof/qualification, vibration/acoustic, thermal-vacuum, modal, separation and integrated GNC.
- [ ] Re-run all multidisciplinary mass/performance loops after higher-fidelity changes.

The research package can be accepted as a reproducible conceptual baseline when its numerical tests pass and each open physical gate has an owner and next evidence item. A vehicle design release requires those physical and program gates to close separately.
