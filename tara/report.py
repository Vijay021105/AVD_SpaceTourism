"""Human-readable evidence report derived from computed values, not templates of success."""
from .io import write_csv

CHAPTERS=[
 ('Mission definition and coupling','2, 7, 8','inputs/mission.csv; tara/pipeline.py','Implemented dependency chain; open gates retained'),
 ('Orbit, J2 and injection','3, 4, 6','tara/orbits.py; tara/rendezvous.py','Numerical ECI J2; spherical surface and GMST'),
 ('Staging and mass allocation','5','tara/sizing.py','Constrained numerical optimization plus hardware dry floors'),
 ('Propulsion','9','inputs/engines.csv; tara/ascent.py','Existing-engine performance proxy; no engine internal design'),
 ('Aerodynamics','10','tara/aerodynamics.py; tara/cfd.py','Engineering database plus inviscid axisymmetric Euler fields'),
 ('Structures and materials','11','tara/structures.py; tara/beam.py','Section cuts, metallic screening, beam modes; detailed shell FEM open'),
 ('Thermal','12','tara/thermal.py','Ascent heating, conduction TPS, heat leak, ideal gas vent sizing'),
 ('Stage auxiliary systems','13','tara/dynamics.py','Separation impulse and tip-off screens'),
 ('Guidance and redundancy','14','tara/ascent.py; inputs/human_rating.csv','Point-mass targeting and illustrative channel reliability; full GNC open'),
 ('Re-entry','15','Not included','Excluded by requested scope')]


def write_report(c,out,s,closure,cases):
    f=s['flight'];r=s['rendezvous'];th=s['thermal']
    text=f'''# TARA_LV computed mission report

This run is a conceptual analysis of a 12,500 kg crew spacecraft launched to the BAS orbital plane, waiting two nominal 400 km periods, then executing a J2-targeted phasing transfer and a prescribed terminal docking approach. Descent and re-entry are outside this run. **Readiness: {s['model_readiness']}.**

## Computed baseline

| Quantity | Result |
|---|---:|
| Liftoff mass | {s['liftoff_mass_kg']/1000:.3f} t |
| Overall conceptual length | {s['height_m']:.3f} m |
| Maximum exterior diameter | {s['diameter_m']:.3f} m |
| Launcher sizing ideal delta-v | {s['sizing']['achieved_dv_mps']:.2f} m/s |
| Spacecraft allocation | {c['mission']['spacecraft_dv_allocation_mps']:.2f} m/s |
| Injection elapsed time | {f['injection_time_s']:.3f} s |
| Injection altitude | {f['final_altitude_m']/1000:.4f} km |
| Injection inertial speed | {f['final_speed_mps']:.4f} m/s |
| Injection radial speed | {f['final_radial_speed_mps']:.4f} m/s |
| Injection inclination | {f['inclination_deg']:.6f} deg |
| Osculating perigee / apogee | {f['perigee_altitude_m']/1000:.2f} / {f['apogee_altitude_m']/1000:.2f} km |
| Peak dynamic pressure | {f['max_q_Pa']/1000:.3f} kPa |
| Maximum proper acceleration | {f['max_proper_g']:.3f} g |
| Doubled-resolution position difference | {f['integration_refinement_position_m']:.3f} m |
| Doubled-resolution velocity difference | {f['integration_refinement_velocity_mps']:.4f} m/s |
| Nominal orbital period | {r['nominal_period_s']/60:.4f} min |
| Two-period coast | {r['two_period_wait_s']/3600:.4f} h |
| BAS phase at departure | {r['actual_phase_at_departure_deg']:.4f} deg |
| Phasing duration | {r['phasing_duration_s']/3600:.4f} h |
| Spacecraft maneuver delta-v | {r['total_spacecraft_dv_mps']:.3f} m/s |
| Spacecraft delta-v reserve | {r['dv_reserve_mps']:.3f} m/s |
| Transfer endpoint miss at 1 km standoff | {r['arrival_position_error_m']:.4f} m |
| Minimum phasing altitude | {r['phasing_min_altitude_m']/1000:.3f} km |
| Prescribed docking elapsed time | {r['docking_time_s']/3600:.3f} h |
| Ascent stagnation heat-flux peak | {th['peak_stagnation_flux_Wm2']/1000:.3f} kW/m2 |
| Ascent TPS thickness | {th['TPS_thickness_m']*1000:.3f} mm |

The 11,400 m/s number is a **sizing budget**, split into launcher ideal capability and spacecraft allocation. It is not orbital speed and is not added to the gravity and drag losses a second time. Actual integrated thrust/mass, gravity projection and drag histories are exported separately. Upper-stage unused propellant is retained on that stage; it is not spacecraft docking propellant.

The BAS initial state is fixed at the common liftoff epoch. The achieved ascent state initializes the crew coast. The transfer is shot to the propagated BAS standoff state under J2. The terminal distance reaches zero because a smooth approach path is prescribed; that endpoint is **not evidence of autonomous docking or contact qualification**.

## Stage ledger

| Stage | Engine | Count | Propellant kg | Dry kg | Ignition T/W | Ideal delta-v m/s |
|---|---|---:|---:|---:|---:|---:|
'''
    for row in s['stages']:
        text+=f'| {row["stage"]} | {row["engine_id"]} | {row["engine_count"]} | {row["propellant_kg"]:.2f} | {row["dry_kg"]:.2f} | {row["ignition_tw"]:.3f} | {row["dv_ideal_mps"]:.2f} |\n'
    text+='\n## Acceptance and open work\n\n| Check | Status | Computed value | Criterion |\n|---|---|---|---|\n'
    for row in s['checks']:text+=f'| {row["check"]} | **{row["status"]}** | {row["value"]} | {row["criterion"]} |\n'
    text+='\n## CFD evidence\n\nThe trajectory uses the engineering database. The Euler snapshots are an independent diagnostic and are **not fed back as a validated drag table**. Slip-wall Euler resolves pressure-wave structure but cannot resolve boundary layers, viscous separation, base recirculation, turbulence, aeroacoustic generation or heat transfer. Axisymmetry at zero incidence cannot identify Cn, Cnalpha or windward/leeward pressure differences. Those outputs use separately labelled engineering estimates.\n\n| Case | Status | Reason |\n|---|---|---|\n'
    for case in cases:text+=f'| {case["case"]} | {case["status"]} | {case.get("reason", "")} |\n'
    if not cases:text+='| Not requested in this run | SKIPPED | Run without --skip-cfd |\n'
    text+='\n## Book-to-implementation map\n\nThe supplied Suresh and Sivan volume is the organizing reference. Its chapter numbering, rather than the unrelated Venus-mission README numbering, is used here. Numerical assumptions are not attributed to the book unless it actually supplies them.\n\n| Topic | Book chapter | Implementation | Scope |\n|---|---|---|---|\n'
    for row in CHAPTERS:text+='| '+' | '.join(row)+' |\n'
    text+='\n## Figures\n\n'
    for name in s['figures']:text+=f'- [{name}](figures/{name})\n'
    text+='\n## Reproducibility\n\nAll 14 subsystem CSV inputs are hashed in `mission_summary.json`. `effective_configuration.json` contains the complete generated model. `data/` contains numeric histories, iteration logs, geometry, margins and open gates. `fields/` contains actual finite-volume arrays. `gifs/` and `videos/` are display products. GIF and MP4 field sweeps change operating condition and do not represent physical unsteady flow.\n\nSee `../docs/WORK_PLAN.md`, `../docs/METHODS.md`, `../docs/DATA_DICTIONARY.md` and `../docs/SOURCES.md` before changing model fidelity or treating a screen as a design release.\n'
    (out/'MISSION_REPORT.md').write_text(text)
    write_csv(out/'data/book_traceability.csv',[dict(topic=a,chapters=b,implementation=d,scope=e) for a,b,d,e in CHAPTERS])
