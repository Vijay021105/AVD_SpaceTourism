# Migration audit

The original archive was inspected before rebuilding. It contains two conflicting workflows, stale generated examples, a root README with unresolved Git merge markers, and independent epoch/geometry assumptions. The unrelated Venus-mission README PDFs are used only as a visual-output reference.

| Original issue | Resolution in TARA_LV |
|---|---|
| MATLAB `mission_config.m` names 2 MVac, uses 348 s, 981 kN and a 6.5 m lower-stage diameter | One shared standard-M1D engine record; fixed 10/2/1 architecture and maximum 4 m |
| Python and MATLAB trajectories are independent | One 3-D Python ECI ascent and one exported achieved injection state |
| Python example states an Earth-intersecting perigee | Ascent re-targeted; actual osculating perigee/apogee and target gates exported |
| Tank coordinates hard-coded to a 75.4033 m meridian | Dimensions regenerated; no fixed overall length |
| 4 m bare tank OD plus added insulation exceeds 4 m | Exterior allowance explicitly subtracts wall/clearance/insulation from tank radius |
| Loaded propellant capacities inconsistent with revised tanks | Capacity is regenerated from the optimized mass, density and ullage |
| Old CE32/CE20/C32 labels mixed | `CE32_project_CE20_proxy` explicitly identifies performance basis |
| Different orbit epochs and manually prescribed crew states | Single UTC epoch; fixed BAS initial state; actual-state rendezvous |
| Three-hour coast or unrelated animation sequence | Two nominal 400 km periods; phasing duration separately computed |
| Frozen structural fractions do not guarantee hardware mass closure | Generated hardware dry floors feed back into stage optimization |
| Axisymmetric pressure CFD used as a broad aero claim | Engineering database, Euler fields and incidence proxies have separate provenance |
| Old figures could survive beside revised inputs | Outputs created under one documented results tree; input hashes and manifest identify this run |
| MATLAB, visualizer catalogs, old docs and `.git` obscure the executable path | New concise package; old archive retained outside this deliverable |

New geometry and all dynamic outputs are calculated; old trajectories, maps or unrelated-mission numerical results are not relabelled as TARA_LV results. A complete redesign may expose infeasibility. Such gates remain open rather than changing user requirements to produce a pass.
