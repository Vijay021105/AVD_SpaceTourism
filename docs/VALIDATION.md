# Verification of the bundled baseline

The full driver completed in the standard profile. The independent physics suite passed 10 tests: exact circular two-body return, J2 node regression, inward gravity, tank volume/envelope, finite-volume ring metrics/static pressure, equal-state numerical flux, stage reserve accounting, nose/heating curvature consistency, zero-heat thermal equilibrium and vent choking/reverse-flow limits.

The 14 nominal CSV input files regenerate byte-for-byte. The final trajectory was re-integrated after removing two unused legacy seed fields and reproduced its stored injection state exactly. Output checks cover stage mass drops, one shared injection vector, timestamps, figure presence, input hashes, finite numeric data and positive CFD density/pressure.

The six admissible standard Euler snapshots meet their residual/force criteria. The grid study is a sensitivity study, not experimental validation or proof of grid independence. Higher-temperature/rarefied or unreached cases are explicitly excluded.

## Open physical and program gates

| Gate | Consequence / next evidence |
|---|---|
| Engine installation | The local packing and assumed swept envelopes do not close within 4 m. Obtain real supplier geometry and resolve layout before geometry freeze. |
| Structural margins | Some shell/buckling/joint screens are negative. Size local reinforcements, joints and stiffeners, then return added mass to staging and ascent. |
| Cryogenic reserve | Pad-plus-ascent LH2 boiloff exceeds the selected 1% allowance under the foam screen. Define topping/loading/pressure-control strategy and verify insulation. |
| Assumed POGO loop | The chosen scalar delayed loop is unstable. The test transient stops at its 10 mm linear-range limit. Identify real engine/feed/structure data and design mitigation. |
| Human rating | Safety/FDIR/abort/escape/ECLSS, qualification and integrated testing are not established. |
| Viscous separated CFD | Full transition/turbulent/separated/base/plume/real-gas/rarefied and incidence solutions are not implemented. |
| Aeroelastic validation | Flexible-flow/control coupling, flutter/divergence and validated buffet response are not established. |

Tests passing means the implemented numerical models satisfy the listed checks. It does not close the gates above. The kinematic ascent/phasing solution is conditional on the idealized vehicle and stated assumptions.
