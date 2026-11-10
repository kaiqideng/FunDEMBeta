# LAMMPS official granular verification coverage

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

All **61** input Git blobs at [29dd3fd7642a6fe9e8b4dd99120df4068b221250](https://github.com/lammps/lammps/tree/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests) are preserved and verified. This scope is the complete pinned granular verification directory, rather than every LAMMPS example elsewhere.

Historically completed inputs: **35** (**33 PASS**, **2 FAIL**); remaining without four fully verified retained runs: **26**. The original seventeen input records, including two failures, remain unchanged. New completed execution requires four native runs, independent physics checks and the complete CRC-verified native archive; dem11 additionally needs its geometry preflight.

Machine-readable source/capability/execution ledger (retained locally; not included) · All exact source hashes and parameters (retained locally; not included) · Immutable initial classification (retained locally; not included) · Model and native capability source hashes (retained locally; not included)

| Official input | Analytic model | Official tolerance | Historical comparison status |
|---|---|---:|---|
| dem01-hertz-tsuji-3d-si.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hertz-tsuji/REPORT.md) |
| dem01-hertz-twosphere-3d-si.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hertz-twosphere/REPORT.md) |
| dem01-hooke-2d-lj.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hooke-2d-lj/REPORT.md) |
| dem01-hooke-massvel-3d-si.yaml (retained locally; not included) | collision_restitution | 0.005 | Executed: PASS; [report](../../../../cases/dem01-hooke-massvel/REPORT.md) |
| dem01-hooke-twosphere-3d-si.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hooke-twosphere/REPORT.md) |
| dem01-hooke-twosphere-e02-3d-si.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hooke-twosphere-e02/REPORT.md) |
| dem01-hooke-twosphere-e08-3d-si.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hooke-twosphere-e08/REPORT.md) |
| dem01-hooke-veldamp-3d-si.yaml (retained locally; not included) | collision_restitution | 0.005 | Executed: PASS; [report](../../../../cases/dem01-hooke-veldamp/REPORT.md) |
| dem01-hookehistory-twosphere-3d-si.yaml (retained locally; not included) | collision_restitution | 0.02 | Executed: PASS; [report](../../../../cases/dem01-hookehistory-twosphere/REPORT.md) |
| dem01-legacy-hooke-twosphere-3d-si.yaml (retained locally; not included) | collision_restitution | 0.005 | Executed: PASS; [report](../../../../cases/dem01-legacy-hooke-twosphere/REPORT.md) |
| dem01-legacy-hookehistory-twosphere-3d-si.yaml (retained locally; not included) | collision_restitution | 0.005 | Executed: PASS; [report](../../../../cases/dem01-legacy-hookehistory-twosphere/REPORT.md) |
| dem02-hertz-peak-3d-si.yaml (retained locally; not included) | hertz_peak | 0.01 | Executed: PASS; [report](../../../../cases/dem02-hertz-peak/REPORT.md) |
| dem02-hertz-twosphere-3d-si.yaml (retained locally; not included) | hertz_normal_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem02-hertz-twosphere/REPORT.md) |
| dem02-hertzmaterial-twosphere-3d-si.yaml (retained locally; not included) | hertz_normal_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem02-hertzmaterial-twosphere/REPORT.md) |
| dem02-legacy-hertzhistory-twosphere-3d-si.yaml (retained locally; not included) | hertz_normal_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem02-legacy-hertzhistory-twosphere/REPORT.md) |
| dem03-hertz-peak-wall-3d-si.yaml (retained locally; not included) | hertz_peak | 0.01 | Executed: PASS; [report](../../../../cases/dem03-hertz-peak-wall/REPORT.md) |
| dem03-hertz-tsuji-wall-3d-si.yaml (retained locally; not included) | wall_restitution | 0.005 | Executed: PASS; [report](../../../../cases/dem03-hertz-tsuji-wall/REPORT.md) |
| dem03-hertz-wall-3d-si.yaml (retained locally; not included) | hertz_normal_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem03-hertz-wall/REPORT.md) |
| dem03-hertzmaterial-wall-3d-si.yaml (retained locally; not included) | hertz_normal_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem03-hertzmaterial-wall/REPORT.md) |
| dem03-legacy-hertzhistory-wall-3d-si.yaml (retained locally; not included) | hertz_normal_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem03-legacy-hertzhistory-wall/REPORT.md) |
| dem04-hertz-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem04-hertz-oblique/REPORT.md) |
| dem04-hookehistory-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Executed: PASS; [report](../../../../cases/dem04-hookehistory-oblique/REPORT.md) |
| dem04-legacy-hookehistory-grazing-3d-si.yaml (retained locally; not included) | energy_dissipation | 1e-06 | Unsupported active law |
| dem04-legacy-hookehistory-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-mindlin-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-mindlin-oblique-highvt-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-mindlin-oblique-lowvt-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-mindlin-sticking-3d-si.yaml (retained locally; not included) | energy_dissipation | 1e-06 | Unsupported active law |
| dem04-mindlinforce-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-mindlinrescale-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-mindlinrescaleforce-oblique-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported active law |
| dem04-oblique-grazing-3d-si.yaml (retained locally; not included) | oblique_impact | 0.01 | Unsupported; CPU history replay gives a force counterexample |
| dem05-hookehistory-slip-3d-si.yaml (retained locally; not included) | slip_cessation | 0.01 | Executed: PASS; [report](../../../../cases/dem05-hookehistory-slip/REPORT.md) |
| dem05-legacy-hookehistory-slip-3d-si.yaml (retained locally; not included) | slip_cessation | 0.01 | Unsupported active law |
| dem05-mindlin-slip-3d-si.yaml (retained locally; not included) | slip_cessation | 0.01 | Executed: FAIL; [report](../../../../cases/dem05-mindlin-slip/REPORT.md) |
| dem05-slip-mu02-3d-si.yaml (retained locally; not included) | slip_cessation | 0.01 | Executed: FAIL; [report](../../../../cases/dem05-slip-mu02/REPORT.md) |
| dem05-sliptime-3d-si.yaml (retained locally; not included) | slip_transient | 0.01 | Executed: PASS; [report](../../../../cases/dem05-sliptime/REPORT.md) |
| dem06-hookehistory-spin-3d-si.yaml (retained locally; not included) | spin_impact | 0.07 | Executed: PASS; [report](../../../../cases/dem06-hookehistory-spin/REPORT.md) |
| dem06-mindlin-spin-heavy-3d-si.yaml (retained locally; not included) | spin_impact | 0.02 | Unsupported; CPU history replay gives a force counterexample |
| dem06-spin-identical-3d-si.yaml (retained locally; not included) | spin_no_friction | 1e-08 | Executed: PASS; [report](../../../../cases/dem06-spin-identical/REPORT.md) |
| dem06-spin-wall-3d-si.yaml (retained locally; not included) | spin_impact | 0.02 | Unsupported; CPU history replay gives a force counterexample |
| dem07-hertz-rolling-3d-si.yaml (retained locally; not included) | rolling_decay | 0.01 | Unsupported active law |
| dem07-hooke-rolling-3d-si.yaml (retained locally; not included) | rolling_decay | 0.01 | Executed: PASS; [report](../../../../cases/dem07-hooke-rolling/REPORT.md) |
| dem07-incline-rolling-3d-si.yaml (retained locally; not included) | incline_rolling | 0.02 | Unsupported active law |
| dem07-incline-static-3d-si.yaml (retained locally; not included) | incline_rolling | 0.0001 | Unsupported active law |
| dem07-rolling-mu02-3d-si.yaml (retained locally; not included) | rolling_decay | 0.01 | Unsupported active law |
| dem08-dmt-twosphere-3d-si.yaml (retained locally; not included) | pulloff_dmt | 0.001 | Unsupported active law |
| dem08-dmt-wall-3d-si.yaml (retained locally; not included) | pulloff_dmt | 0.001 | Unsupported active law |
| dem08-jkr-twosphere-3d-si.yaml (retained locally; not included) | pulloff_jkr | 0.005 | Unsupported active law |
| dem09-viscous-linear-3d-si.yaml (retained locally; not included) | terminal_velocity_linear | 0.005 | Executed: PASS; [report](../../../../cases/dem09-viscous-linear/REPORT.md) |
| dem09-viscous-nonlinear-3d-si.yaml (retained locally; not included) | terminal_velocity_schiller_naumann | 0.005 | Unsupported active law |
| dem10-freefall-3d-si.yaml (retained locally; not included) | freefall | 1e-09 | Executed: PASS; [report](../../../../cases/dem10-freefall/REPORT.md) |
| dem10-stack-3d-si.yaml (retained locally; not included) | stack_energy | 0.0005 | Executed: PASS; [report](../../../../cases/dem10-stack/REPORT.md) |
| dem11-regionwall-3d-si.yaml (retained locally; not included) | wall_restitution | 0.01 | Executed: PASS; [report](../../../../cases/dem11-regionwall/REPORT.md) |
| dem12-superellipsoid-3d.yaml (retained locally; not included) | momentum_conservation | 1e-08 | Unsupported active law |
| dem13-twist-marshall-3d-si.yaml (retained locally; not included) | twist_decay_marshall | 0.02 | Unsupported active law |
| dem13-twist-sds-3d-si.yaml (retained locally; not included) | twist_decay | 0.01 | Executed: PASS; [report](../../../../cases/dem13-twist-sds/REPORT.md) |
| dem14-heat-area-3d-si.yaml (retained locally; not included) | heat_equilibration | 0.005 | Unsupported active law |
| dem14-heat-radius-3d-si.yaml (retained locally; not included) | heat_equilibration | 0.005 | Unsupported active law |
| dem15-hertz-oblique-3d-si.yaml (retained locally; not included) | oblique_impact_pair | 0.01 | Unsupported active law |
| dem15-hooke-oblique-3d-si.yaml (retained locally; not included) | oblique_impact_pair | 0.01 | Executed: PASS; [report](../../../../cases/dem15-hooke-oblique/REPORT.md) |

## Mapping scope and limitations

Only the zero-tangential-velocity trajectories in the six extra normal-only mappings permit dormant tangential differences. Their execution does not verify an active tangential-history law. Tsuji dashpots use an exact polynomial coefficient conversion; their material restitution input is an encoding of that coefficient, not a fitted rebound.

The dem11 curved cylinder is retained. Its LAMMPS effective radius is Rp*8/7; E*=(3/4)kn sqrt(8/7) preserves the active Hertz force/dashpot coefficients. Its sampled SDF introduces an explicitly measured geometry error that requires independent mesh preflight.

Unsupported active-law groups remain explicit in the JSON: coupled Hertz/independent kt-or-kr switches; active overlap-dependent tangential dashpot ratios; force/rescale or legacy tangential-history updates; DMT/JKR cohesion; Schiller-Naumann drag; LS-LS area weighting versus a single pair stiffness; Marshall twisting; and heat conduction. No force fitting or solver modification is used to turn these into apparent passes.

Official tolerance meanings are model-specific. spin_no_friction includes absolute transverse speeds, static incline branches use absolute speed bounds, energy_dissipation allows a fractional energy excess, and momentum_conservation uses defined initial momentum scales. Per-input formulas/check step/time are in the JSON. Same-code epsilon1e-10 is not substituted for cross-code acceptance.

## Additional source-level counterexamples and corrections

[Common-cap first-step diagnostics](cap-degeneracy-diagnostics.md) quantify why the two flat Hertz rolling inputs cannot hide their rolling stiffness mismatch behind Coulomb saturation. Active tangential-damping probe candidates remain excluded until an all-step source-history replay proves equality; analytic endpoint agreement alone is insufficient.

Wall modulus mapping correction (retained locally; not included) records that the pinned helper uses E_wall=E/[2(1-nu^2)]. The corrected Tsuji wall case requires E_native=E_source for both materials. Its superseded doubled-E runs are diagnostic evidence and are not counted as source-faithful validation.

Unsupported means that the exact source-specified contact algorithm or configuration combination cannot be reproduced faithfully by the current native interfaces. It does not mean FunDEM cannot simulate that physical class of motion. The three [actual CPU history-replay diagnostics](conditional-history-probe-results.md) refute constitutive equality; their proof failures are not official endpoint-criterion failures and are not added to the two preserved official FAIL cases.

The Hertz/constant-kt/Ct0 oblique input has a narrowly proved exception to the shared-stiffness limitation: total Ft is identical only when every actual source/native contact row stays on the same cap. Its per-mode original-history replays are linked in the JSON. This does not validate unsaturated kt, tangential damping or identical internal history, and source R/native R-overlap wall torque geometry remains disclosed.
