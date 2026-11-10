# LAMMPS official dem01-hooke-twosphere-e02

> Historical evidence from 28–29 October 2026, migrated on 30 October 2026. The outcome and numerical results below are preserved, not a new run of the current core. Historical raw-frame/archive retention statements describe the original campaign, not files distributed here.

**Normal collision: restitution and momentum**  
Date: 29 October 2026. Four main runs: **PASS**. Official input tag: `unvetted`.

## Source and official acceptance

[Pinned official YAML](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hooke-twosphere-e02-3d-si.yaml) | [Analytic checker](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp). Commit 29dd3fd7642a6fe9e8b4dd99120df4068b221250; input blob 2e58af1d814647620c46816fc6429f9e612f7b1a. LAMMPS was not rerun locally. The YAML's same-code epsilon is not a cross-software tolerance.

Model `collision_restitution` is evaluated at original step 2795 and time 8.385e-05s. Half-step runs use twice the native check step at the same physical time. The exact reference, relative/absolute error type and unchanged limit are listed separately for each quantity. Whole histories and source golden states supplement these selected-time assertions.

| Run | Quantity | Reference | Native value | Error type | Error | Limit | Result |
|---|---|---:|---:|---|---:|---:|---|
| CPU dt | Total x momentum | 0 | 0 | absolute | 0 | 1.41371669412e-05 | PASS |
| CPU dt | Restitution | 0.2 | 0.199943294792 | relative | 0.00028352603757 | 0.02 | PASS |
| GPU dt | Total x momentum | 0 | 0 | absolute | 0 | 1.41371669412e-05 | PASS |
| GPU dt | Restitution | 0.2 | 0.199943294792 | relative | 0.000283526037569 | 0.02 | PASS |
| CPU dt/2 | Total x momentum | 0 | 0 | absolute | 0 | 1.41371669412e-05 | PASS |
| CPU dt/2 | Restitution | 0.2 | 0.199682263488 | relative | 0.00158868256169 | 0.02 | PASS |
| GPU dt/2 | Total x momentum | 0 | 0 | absolute | 0 | 1.41371669412e-05 | PASS |
| GPU dt/2 | Restitution | 0.2 | 0.199682263488 | relative | 0.00158868256169 | 0.02 | PASS |

## Inputs and constitutive mapping

- left: R=0.0025 m; position=(-0.00255, 0, 0) m; velocity=(2, 0, 0) m/s; spin=(0, 0, 0) rad/s.
- right: R=0.0025 m; position=(0.00255, 0, 0) m; velocity=(-2, 0, 0) m/s; spin=(0, 0, 0) rad/s.

- left: mass=0.000176714586764 kg; isotropic inertia=4.41786466911e-10 kg m^2.
- right: mass=0.000176714586764 kg; isotropic inertia=4.41786466911e-10 kg m^2.

Linear Hooke Fn_elastic=kn*delta, kn=1000000 N/m. Automatic damping uses the harmonic material restitution and the active normal/tangential stiffness; stored manual dashpots are inactive. Tangential inputs: kt=0 N/m (explicit in linear mode), mu=0. Rolling/twist coefficients: kr=0, mu_r=0, ktwist=0, mu_twist=0.

Gravity=(0, 0, 0)m/s^2. Base/half dt=3e-08/1.5e-08s; base/half numerical steps=2795/5590. Initial poses, velocities and physical source coefficients are preserved. The parameter mapping translates the source model's coefficient definitions; it does not fit measured native results.



Native specification (retained locally; not included) | Parameter mapping (retained locally; not included) | Configuration round-trip (retained locally; not included)

## CPU/GPU and timestep refinement

| Quantity | CPU dt error | CPU dt/2 error | Error type and limit |
|---|---:|---:|---|
| Total x momentum | 0 | 0 | absolute: 1.414e-05 |
| Restitution | 0.00028353 | 0.0015887 | relative: 0.02 |

CPU/GPU whole-history maximum position / velocity / spin / force differences: 0 m / 1.11e-16 m/s / 0 rad/s / 7.11e-15 N.

Maximum individual impulse / position-drift / velocity-kick residuals: 5.96e-19 kg m/s / 2.17e-19 m / 2.12e-16 m/s. Gravity is added separately to saved native contact forces.

Energy histories are supplemental unless selected by the official analytic model. Wall/gravity cases exchange external momentum; total sphere momentum is not assumed conserved. The step-refinement table reports each error directly, without imposing monotonic convergence.

Numerical criteria and measurements are retained below and in the companion CSV tables.

## Complete retained results

All 16,774 native output frames, every integration-step particle/contact CSV, projects, raw initial samples and backend logs are retained. Native output cadence is distinct from integration dt; exact counts/intervals are in the analysis and archive manifests. The four primary frame counts are [2796, 2796, 5591, 5591]. 0additional diagnostic run is retained.

- **CPU dt** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 2,795 integration steps, 2,796 stored frames.
- **GPU dt** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 2,795 integration steps, 2,796 stored frames.
- **CPU dt/2** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 5,590 integration steps, 5,591 stored frames.
- **GPU dt/2** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 5,590 integration steps, 5,591 stored frames.

Analysis JSON (retained locally; not included) | Official metrics CSV (retained locally; not included) | Source checkpoints (retained locally; not included) | Independent physics audit (retained locally; not included) | Run manifest (retained locally; not included) | Full native ZIP (retained locally; not included) | Frame/hash audit (retained locally; not included)

Source checkpoint states are supplemental regression diagnostics. Their differences are reported directly; the YAML same-code epsilon is not imposed as an independent cross-software acceptance limit.

A stable helper executable is statically linked to the current native CPU/CUDA libraries. Historical GUI/library byte identity is not asserted. Existing case reports and records are unchanged. Verified scope follows the active source model and the criteria above; a passing inactive channel does not validate every option of its contact family.

## Repository evidence and native rerun input

- [Dedicated C++ case source](dem01-hooke-twosphere-e02.cpp). The category README documents the runner and model-mapping limits.
- [Compact historical numerical summary](HISTORICAL-RESULTS.md).
- [Source pins, original outcome and retained-artifact hashes](PROVENANCE.md).
- [Compact checkpoint table: official_metrics_comparison.csv](checkpoints/official_metrics_comparison.csv).
- [Compact checkpoint table: source_checkpoint_comparison.csv](checkpoints/source_checkpoint_comparison.csv).

A new runner result must be compared independently with these targets; opening or completing the input is not an automatic PASS. Full native frames, per-step histories, execution logs, vendor inputs, proprietary programs and PDF reports are not included.
