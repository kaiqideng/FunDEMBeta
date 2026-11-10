# LAMMPS active-tangential-history diagnostic exclusions

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Three original official inputs were run as single-CPU diagnostic probes. They are **not validated cases**, and no GPU or half-step run was promoted for them. Every contact row was replayed with the original official history and `ct=0.5*cn`; endpoint agreement was not used as a substitute.

The normal force prefactor and signed unloading match the pinned source law. Both models use `mu*abs(Fn)` as the cap. The official model caps the total spring-plus-dashpot trial and stores a history that includes damping compensation; native caps elastic history first and limits damping within the remaining Coulomb disk. This difference is accessed by the actual trajectories.

| Official input | Replayed contact rows | Source/native unsaturated rows | First unequal step | Source Ft_x (N) | Native Ft_x (N) |
|---|---:|---:|---:|---:|---:|
| dem04-oblique-grazing | 1480 | 4/0 | 2237 | -0.00453817992281 | -0.0302570775012 |
| dem06-spin-wall | 1047 | 36/0 | 51 | -353.329957106 | 353.329957106 |
| dem06-mindlin-spin-heavy | 1086 | 37/0 | 51 | -337.560121623 | 337.560121623 |

The replay uses the actual saved native trajectory as its prescribed state path. Before the first unequal force, this identifies a concrete accessible constitutive counterexample; it is not claimed to be a separate LAMMPS execution. Subsequent shadow-history rows are additional diagnostics. Source wall torque uses R while native uses R-overlap; pair torque uses its actual shared point. Slip-velocity differences, signs, signed normal forces and every source/native force component are retained in each all-contact CSV.

The spinning cases first touch at extremely shallow positive overlap; the source damping-inclusive history back-calculation can produce a counteracting elastic history and then a force opposite to the native capped force. The grazing case has four source unsaturated rows despite persistent native saturation. Thus gross macroscopic sliding does not prove both algorithms remain on the same Coulomb cap.

Machine-readable proof summary (retained locally; not included) · [First-step rolling counterexamples](cap-degeneracy-diagnostics.md) · Complete coverage ledger (retained locally; not included)

The exclusion concerns faithful reproduction of the exact source-specified constitutive algorithm/configuration. It does not mean FunDEM cannot simulate these physical spinning or grazing conditions. A failed history-equality proof is distinct from a failed official endpoint analytic criterion; this diagnostic does not claim the latter.
