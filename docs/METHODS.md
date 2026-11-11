# Methods and fidelity

## Mission and units

SI units are used internally: m, s, kg, N, Pa, K. Input angles are in degrees and explicitly converted. ECI is an equatorial inertial approximation oriented through a conventional GMST rotation at the input UTC epoch. The surface and launch-site coordinates are spherical with radius 6,378,137 m. This is not a full ITRF/EOP/precession/nutation implementation.

The baseline is a conditional design point, not a statement that the engine cluster, material allowables or human-flight requirements are qualified. In particular, the 4 m envelope, ten lower-stage engines and the chosen human-flight allowances may expose incompatible constraints.

## Staging and dry-mass feedback

The book's structural factor convention is followed:

\[
\epsilon_i={m_{d,i}\over m_{d,i}+m_{p,i}},\qquad
m_{d,i}=\max\left({\epsilon_i\over1-\epsilon_i}m_{p,i}+m_{fixed,i},\;m_{hardware,i}\right).
\]

The selected epsilons are **assumptions** bounded by sensitivity inputs, not measured values from the book. Fixed allowances cover secondary systems and redundancy inside the dry ledger. The estimated hardware floor is regenerated after each propellant optimization; it includes tanks, insulation, baffles, engines, frames, skirts, stiffener allowance and joint rings. Unmodeled detailed hardware is not thereby certified to fit the remaining budget.

For each stage, upper-stage stacks and the spacecraft remain attached during its burn:

\[
\Delta v_i=g_0 I_{sp,i}\ln{m_{upper,i}+m_{d,i}+m_{p,i}\over
m_{upper,i}+m_{d,i}+f_{reserve,i}m_{p,i}}.
\]

SLSQP minimizes liftoff mass in logarithmic propellant variables with launcher delta-v and ignition T/W inequalities. It reports feasibility independently of the optimizer termination message. Minimum T/W values are project choices, not human-rating standards. A high cryogenic-stage T/W requirement can make the mass/performance problem infeasible; lower T/W increases trajectory lofting and gravity-loss requirements. The generated data retain this interaction.

The allocation is 11,100 m/s launcher ideal velocity and 300 m/s spacecraft capability. The trajectory may cut S3 before exhausting usable propellant. Remaining launcher propellant remains on S3 at spacecraft separation and cannot be spent by the spacecraft. No loss term is added to the sizing budget again.

## Propulsion and geometry

At full throttle, engine mass flow follows \(\dot m=F_{vac}/(g_0I_{sp,vac})\). Effective exit area for M1D follows \(A_e=(F_{vac}-F_{SL})/p_{SL}\), and thrust at altitude is \(F=\eta(F_{vac}-p_aA_e)\). A fixed effective area for the cryogenic proxy follows its assumed nozzle diameter. This is a performance interpolation, not a nozzle, turbopump or combustion calculation.

Fuel and oxidizer masses use O/F. Loaded tank volume is \(V=m/[\rho(1-f_u)]\). Two 2:1 ellipsoidal heads and a cylindrical barrel obey:

\[
V=\pi r^2 L_b+{4\pi\over3}r^2h_d,\qquad h_d=0.5r.
\]

When that head pair exceeds the required volume, radius shrinks and barrel length becomes zero. Maximum bare tank radius is derived by subtracting exterior skin, radial clearance, insulation and wall thickness from the 2 m radius. Wall sizing iterates pressure, hydrostatic full-column pressure at an assumed 1.8 g, and thin-cylinder hoop stress. Pressure factors and welded-metal allowables are provisional.

The head area uses an oblate-spheroid area approximation. Engine layout uses a local constrained circle-packing search. It then separately checks pairwise and radial swept gimbal envelopes. That local search does not prove global packing impossibility or installation feasibility. Dimensions must be confirmed with supplier interfaces.

## Ascent equations and controls

The state is integrated in ECI under central plus J2 gravity:

\[
\dot{\boldsymbol r}=\boldsymbol v,\qquad
\dot{\boldsymbol v}=\boldsymbol g_{J2}+{F\hat{\boldsymbol b}+\boldsymbol F_a\over m},\qquad
\dot m=-\eta\dot m_{full}.
\]

Air-relative velocity is \(\boldsymbol v_a=\boldsymbol v-\boldsymbol\omega_E\times\boldsymbol r-\boldsymbol v_{wind}\). Dynamic pressure is \(q=\rho |v_a|^2/2\). Drag is along \(-\hat v_a\); the normal force lies in the plane spanned by the body axis and velocity. The current normal-force model is a small-incidence engineering screen.

RK4 uses consumed-propellant fraction as its independent variable, giving exact prescribed burn masses. Coast segments use time-domain RK4. Optimization uses a coarse fraction grid; the final targeting correction and published results use 300 intervals per full burn, and 600 intervals test integration sensitivity. Exported time intervals are consequently nonuniform.

Guidance has a vertical hold, piecewise pitch versus propellant fraction, stage coast time, and S3 usable-fraction cutoff. A proportional plane-tracking law includes the nominal J2 nodal precession. High-q incidence relief rotates the commanded thrust axis toward the air-relative velocity; throttle relief responds to q. Available throttle is bounded below by engine minimum throttle and above by full thrust. Limits are checked afterward because a saturated command cannot guarantee them. Pitch-rate output is a numerical derivative of commanded thrust elevation; no attitude inertia, actuator, TVC rate or flexible-body controller is simulated.

The atmosphere is ISA through approximately 86 km followed by a declared log-density screening interpolation. Wind is a user-editable assumed east/north profile. It is not launch-day weather. Perfect-gas heating/flow methods are screened at high temperature and rarefaction; drag in the rarefied tail is still only an engineering estimate with small dynamic-pressure contribution.

## Orbit and rendezvous

Numerical DOP853 integration uses Cartesian central/J2 acceleration. The nominal 400 km period is \(T=2\pi\sqrt{a^3/\mu}\); the two-period wait is exactly \(2T\), even though a perturbed osculating orbit has slightly different crossing periods.

At departure the signed BAS phase is measured in the crew RTN frame. A period-based phasing estimate gives transfer time and an initial tangential impulse. The arrival has a separately specified fractional revolution, avoiding a near-return singularity for the slightly non-coplanar transfer. Three-component departure velocity is then solved by least squares so numerical J2 propagation reaches the BAS along-track standoff point. An arrival impulse matches the rotating standoff velocity. This is a fixed-time multi-revolution shooting solution; no impulsive maneuver is assumed free of propellant cost.

The terminal path uses a cubic along-track displacement with zero endpoint relative velocity. RTN-basis motion and J2 gravity differences give an inverse-dynamics acceleration/thrust requirement. The approach endpoint matches by construction. Sensors, collisions, plume contamination, keep-out zones, navigation uncertainty, contact mechanics and flight-software behavior are open work. Spacecraft maneuver propellant follows the rocket equation using spacecraft initial mass and an assumed 320 s Isp.

## Aerodynamic coefficients and computed flow

The engineering database combines laminar and turbulent flat-plate skin-friction estimates, exposed wetted area, assumed base/form terms and a bounded smooth transonic/supersonic wave-drag estimate. Cnalpha is a slender-body-inspired proxy with a transonic enhancement; Cn adds a small crossflow contribution. These are **uncalibrated vehicle-level estimates**. Coefficients use \(A_{ref}=\pi D^2/4\):

\[
C_D=C_A\cos\alpha+C_N\sin\alpha,\qquad C_{N\alpha}=\partial C_N/\partial\alpha\quad[1/\mathrm{rad}].
\]

The normal slope column denotes the linear term; finite-incidence crossflow makes the full local derivative slightly different. Xcp is a forebody approximation referenced aft from the nose; rolling coefficient is zero by imposed symmetry. Real engine/protuberance asymmetry is not resolved.

The inherited conservative Euler kernel was retained after metric/conservation tests and connected to the regenerated geometry. It integrates axisymmetric ring volumes with exact straight-edge face metrics and the radial pressure source. First-order Rusanov flux with positivity step halving and local pseudo-time stepping is used by default. Pseudo-time is not physical unsteady time.

Convergence requires both residual and pressure-force change limits for consecutive checkpoints. Output pressure coefficient excludes base wake and viscous drag. The mesh terminates at the current base plane; no wake recirculation or plume is present. An iteration cap is reported as `NOT_CONVERGED`. High total-temperature, rarefied, unreached-Mach and extrapolated-atmosphere cases are `EXCLUDED`. Grid sensitivity is reported as observed change; no asymptotic-order or grid-independence claim is made.

CFD pressure does not replace the trajectory database in this release. Windward/leeward incidence plots use separate analytical proxies. Full laminar/turbulent fields, transition, separated flow, 3-D force moments, unsteady aeroacoustic sources and hypersonic real-gas/rarefied fields remain unsolved. They require a validated external CFD campaign.

## Loads and structures

Dry masses, engine masses, propellant columns and spacecraft masses are distributed over their attached components with exact mass sums. Each ascent sample generates axial inertia/drag resultants and lateral aerodynamic loading. A single-plane quasi-static TVC trim force balances moment about the CG; translational inertial relief closes total lateral force. Section shear and bending are integrated from the free nose. Selected simultaneous cases are exported in addition to componentwise envelopes.

Tank membrane hoop stress, axial compression/bending and shear form a von Mises screen. Yield and ultimate factors apply to limit stress. Classical axial shell buckling has a separate assumed knockdown; pressure stabilization is not credited. Joint bolts use equal sharing plus a simple moment couple. Envelopes may combine demands from different times and are conservative screening combinations. Local shell modes, fitting/prying, cutouts, weld variability, thermal stress, fracture and fatigue require detailed analysis.

The hardware item mass floor is coupled back to staging. Reinforcement needed to fix negative load margins is **not** automatically added; negative margins therefore keep the design open and require another sizing loop after detailed design.

## Dynamics and human-flight provisions

Free-free Timoshenko beam modes use homogenized sections, consistent mass, rotary inertia and reduced shear integration. The three rigid modes are explicitly projected out. Wet liquid is rigidly attached; a separate first lateral cylindrical slosh estimate varies with fill. Baffle damping is an assumption, not a solved baffle-flow property. Pad modes use a clamped-base version of the matrices.

Acoustic and buffet PSDs are normalized to specified assumed RMS pressure. A one-mode response transfer estimates acceleration PSD; no acoustic-source generation, chamber transmission or SEA validation is performed. A delayed scalar feed/structure loop is checked by exact-delay Nyquist analysis and root residuals. Its assumed poles/gain are not an identified POGO model. A method-of-steps free transient uses a 0.1 mm test disturbance and stops at a 10 mm equivalent-displacement screen instead of extending a linear instability to unphysical amplitudes.

Ground wind uses a power-law profile, an assumed cylinder drag coefficient and gust multiplier to obtain cantilever shear, moment and static deflection. It does not predict lock-in dynamics. Cold separation uses equal/opposite prescribed push forces and an assumed angular tip-off; 6-DOF re-contact is not established.

Redundant-channel reliability uses binomial independent-channel failure plus an illustrative beta common-cause term. Assumed probabilities are placeholders for sensitivity, not measured reliability or crew loss probability. Single-engine upper-stage vulnerabilities, abort capability, engine isolation, software failures, ECLSS and program hazard acceptance remain open.

## Thermal and venting

The ascent stagnation correlation is Sutton–Graves in Earth-air SI form, \(\dot q=k\sqrt{\rho/R_n}V_a^3\). It is turned off above the declared continuum-heating range rather than being applied to the orbital rarefied tail. A cosine exposure proxy distributes the maximum flux over capsule angle. Base exposure is an assumed fraction of engine mixture power and is not a plume solution.

A 31-node one-dimensional finite-volume thermal column includes TPS heat capacity, a 3 mm aluminium backing heat capacity, implicit conduction and linearized surface radiation. Thickness is bisected to meet the maximum bondline temperature at the sampled ascent times. The inner surface is adiabatic. Material properties are constant. The output is an ascent-only conduction TPS estimate; ablation, re-entry, local attachment temperatures and complete TPS mass closure are outside this model.

Cryogenic heat leak follows \(Q=kA(T_a-T_l)/d\) and an energy-limited latent-heat boiloff estimate over pad hold plus ascent. Pad boiloff assumes topping/loading management to preserve the stated liftoff loaded mass; reserve adequacy is a screen and should be replaced by a loading/pressurization simulation. Ideal-gas choked/subcritical orifices size vents with assumed upstream temperature. Bay vent area follows ambient depressurization demand at an allowed differential pressure; the model is quasi-steady, not a transient two-phase or controlled-valve simulation.
