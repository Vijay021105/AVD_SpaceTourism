# TARA_LV — integrated conceptual crew launch design

TARA_LV rebuilds the supplied BAS_2D_Research work as one Python launch-to-docking pipeline. Its mission is to launch a **12,500 kg crew spacecraft** toward BAS at **400 km / 51.6 degrees**, coast for **two nominal orbital periods**, perform a phasing maneuver, arrive at a 1 km standoff, and follow a prescribed terminal docking approach. The exterior diameter is limited to **4 m**. Descent and re-entry are excluded.

**Start with [`results/MISSION_REPORT.md`](results/MISSION_REPORT.md).** It contains the computed baseline, numerical checks, open design gates, and links to every figure. The results are conceptual engineering studies. Human-rating qualification, viscous separated-flow CFD, detailed structural qualification and autonomous docking are not established by this package.

## Run

Python 3.10 or later is required. NumPy, SciPy, Matplotlib and Pillow provide all mandatory calculations and GIFs. An installed `ffmpeg` executable enables MP4 export.

```powershell
cd TARA_LV
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe full_mission_design.py
```

Linux/macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python full_mission_design.py
```

All paths resolve relative to the driver, so the program can also be launched from another directory. No MATLAB, Windows-specific absolute path, internet access or external ephemeris service is required at runtime. Use a fresh virtual environment to avoid older NumPy builds without `trapezoid`.

| Command | Purpose |
|---|---|
| `python full_mission_design.py` | Complete standard analysis, Euler snapshots/grid sensitivity, figures, GIFs and available MP4s |
| `python full_mission_design.py --profile quick` | Same chain with a smaller Euler mesh and iteration budget; convergence is still checked |
| `python full_mission_design.py --skip-cfd --no-animations` | Mission, sizing, loads, thermal and dynamics without Euler or media generation |
| `python full_mission_design.py --reuse-ascent` | Reuse stored guidance only when all input CSV and ascent-physics hashes match |
| `python full_mission_design.py --output results_trial` | Write a separate experiment |
| `python full_mission_design.py --strict` | Return exit code 2 when any design gate is open |
| `python -m unittest discover -s tests -v` | Independent physics and limiting-case checks |
| `python redraw_figures.py results` | Redraw all PNGs from stored CSV/NPZ data without solving again |
| `python verify_results.py results` | Output integrity, mass staging, common injection state and provenance checks |
| `python regenerate_inputs.py --inputs inputs_copy` | Create a fresh copy of the reviewed nominal CSV deck |

Exit code 0 means the requested calculation completed and wrote its evidence. It does **not** mean all engineering gates passed. Input/integration failures raise an error; `--strict` distinguishes completed calculations with open engineering gates. Standard Euler calculations may take considerably longer than the rest of the pipeline; actual runtime and library versions are recorded in `mission_summary.json`.

`regenerate_inputs.py --overwrite` intentionally resets the CSV deck. The mission driver never invokes it and never silently overwrites edited inputs.

## Configuration authority

| Item | Adopted baseline |
|---|---|
| S1 | 10 standard Merlin 1D, LOX/RP-1 |
| S2 | **2 standard Merlin 1D, not Merlin Vacuum** |
| S3 | 1 `CE32_project_CE20_proxy`, LOX/LH2 |
| Merlin performance | 845 kN sea-level thrust; inherited provisional 282 s sea-level / 311 s vacuum specific impulse |
| Cryogenic performance | CE20-based 22 tonne-force / 442.5 s proxy; `CE32` is a project label |
| Total sizing budget | 11,400 m/s: 11,100 launcher ideal capability + 300 spacecraft allocation |
| Loaded stage masses | Re-optimized; old 340/110/24 tonne values are seeds only |
| Dry masses | Structural-factor estimates with generated hardware-mass lower bounds |
| Tanks | Regenerated from stage propellant mass, O/F, density, ullage, ellipsoidal heads, pressure wall and insulation |
| Epoch | One common `mission.csv` UTC epoch; BAS initial phase is fixed before ascent |
| Launch site | Inherited Kulasekarapattinam vicinity; numerical site only |
| Phasing revolutions | 20 plus a quarter-turn arrival offset; explicit design choice independent of the two-period coast |

The two Merlin stages share the **same engine record**. Vacuum thrust follows the same mass-flow/impulse assumptions; no 348 s MVac performance or 981 kN MVac thrust is imported into S2. Sea-level Merlin ignition in vacuum remains an integration/qualification assumption. Engine dimensions, gimbal swept volume, dry mass and O/F need supplier confirmation. Engine packing failure remains a reported failure rather than being drawn as an accepted installation.

The supplied SpaceX and ISRO engine heritage does not qualify their integration into TARA_LV. The crew exterior is a concept capsule envelope; it is not a reverse-engineered production Dragon drawing.

## Repository layout

| Path | Contents |
|---|---|
| `full_mission_design.py` | Public command-line entry point |
| `inputs/` | The 14 authoritative subsystem CSVs |
| `tara/` | Small, separated physics, coupling, reporting and plotting modules |
| `tests/` | Independent conservation and reference-solution checks |
| `docs/WORK_PLAN.md` | Complete work sequence, acceptance evidence and remaining higher-fidelity tasks |
| `docs/METHODS.md` | Equations, frames, coefficients, numerical methods and model limits |
| `docs/DATA_DICTIONARY.md` | Input/output contracts and units |
| `docs/VALIDATION.md` | Tests performed and open physical/program gates |
| `docs/SOURCES.md` | Book mapping, primary references and provisional assumptions |
| `docs/MIGRATION_AUDIT.md` | Original inconsistencies and their resolution |
| `results/data/` | All generated subsystem CSVs, histories, margins, iteration logs and DXF profile |
| `results/figures/` | Professional PNG figures |
| `results/fields/` | Actual finite-volume solution arrays and grid coordinates in NPZ form |
| `results/gifs/` | GIF animations |
| `results/videos/` | MP4 animations, when ffmpeg is available |
| `results/MISSION_REPORT.md` | Computed report and figure index |
| `results/mission_summary.json` | Machine-readable results, acceptance gates, input hashes and runtime versions |

The original `.git`, merge-conflicted README, old MATLAB drivers, hand-edited tank tables, unrelated mission outputs and previously generated trajectories are not part of the new executable path. The supplied source ZIP remains the historical reference.

## Data and dependency order

1. Read all CSVs, reject unknown/missing tables, check engine count/type, units-by-contract and payload mass sum.
2. Optimize stage propellant allocations with ignition T/W constraints and the ideal-velocity requirement.
3. Generate tanks, structural item estimates and installation dimensions. Feed hardware mass floors back into staging until the item budget closes or the iteration limit is reached.
4. Build a stage-specific engineering aerodynamic database. This exact model is evaluated within ascent integration.
5. Numerically target the 3-D ECI ascent, including J2 gravity, rotating atmosphere, wind, q-relief, incidence relief, throttle floors and acceleration limiting. Report actual violations if the available actuators cannot enforce a limit.
6. Export the achieved injection state. Propagate that state and the independently initialized BAS with J2, coast for two periods, then solve the multi-revolution transfer by shooting.
7. Apply a velocity-match impulse at the 1 km standoff and calculate the acceleration needed to follow the specified smooth final approach.
8. Compute all-sample section-cut envelopes, metallic/buckling/joint margins, mode shapes, assumed acoustic/buffet response, POGO, slosh, ground wind and separation screens.
9. Compute ascent heating and transient TPS thickness, cryogenic heat leak/boiloff, and quasi-steady vent requirements.
10. Run independently labelled Euler snapshots and grid sensitivity, render figures/media, and publish pass/open gates.

Each downstream discipline reads the same in-memory configuration generated from the CSVs. Geometry and dry-mass floors feed back to sizing. The current solver does **not** automatically repair structural margins by adding reinforcement, use Euler coefficients as validated trajectory drag, or couple turbulent unsteady flow to flexible structures. Those are explicit next design loops.

## What the plots mean

The Mach plots of Ca, Cd, Cn, Cnalpha, Xcp and roll cover the engineering screening range. Ca and Cn are body-axis coefficients; Cd is their projection onto the air-relative velocity. Cn is evaluated at the stated incidence, and Cnalpha is per radian. Zero rolling moment is imposed by the smooth symmetric-body assumption. It does not establish roll-control authority.

Euler fields are computed conservative finite-volume iterates with actual residual and force histories. The solver has slip walls, zero incidence, a perfect gas and no base wake domain. It can show pressure waves and shocks within those limitations. Laminar and turbulent friction corrections are flat-plate estimates, not separate laminar/turbulent flow-field solutions. Windward/leeward incidence pressure is an analytical proxy in a separate CSV and figure panel.

Cases beyond the temperature, continuum, altitude-model or reached-Mach bounds are **excluded**, and cases that hit the iteration budget are **not converged**. Their fields can be inspected diagnostically but cannot establish validated loads. No flow-field interpolation is presented as a 3-D RANS/URANS solution.

Mode shapes come from a free-free homogenized Timoshenko beam; ground modes use a clamped base. Slosh, baffles, POGO and pressure excitation use declared reduced-order assumptions. Negative structural margins and unstable assumed loops remain visible. Ascent TPS thickness does not size a re-entry heat shield.

## Human-flight provisions and qualification boundary

The model includes separate yield/ultimate/pressure factors, load uncertainty, buckling knockdown, propellant reserve, hardware mass allowances, multi-channel avionics reliability examples, common-cause sensitivity, acceleration/q limits and a spacecraft maneuver allocation. The reliability arithmetic is illustrative; it is not a probability of crew loss. Allowance masses are included inside stage dry budgets, not added again after sizing.

NASA-STD-5001B and NASA-STD-8719.29 are reference frameworks; a single factor of safety or channel count is not human-rating compliance. The work plan retains hazard analysis, failure detection/isolation/recovery, safe ascent aborts, escape, ECLSS, fracture control, engine-out controllability, range authorization, inspection and qualification testing as open work. Descent analysis is omitted as requested; that omission also means end-to-end crew safety cannot be established by this ascent-only mission model.

## Review and change procedure

Edit CSV inputs, use a new output folder, run the driver and tests, then compare `design_checks.csv`, `stages.csv`, `tanks.csv`, `dry_mass_closure.csv`, `injection_state.csv` and `rendezvous_burns.csv`. Reuse guidance only with matching hashes. Numerical convergence, model verification, physical validation and program qualification are separate forms of evidence; the reports preserve that distinction.
