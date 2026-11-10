# LAMMPS official dem03-hertz-peak-wall

> Historical evidence from 28–29 October 2026, migrated on 30 October 2026. The outcome and numerical results below are preserved, not a new run of the current core. Historical raw-frame/archive retention statements describe the original campaign, not files distributed here.

**Hertz compression and force at the source check time**  
Date: 29 October 2026. Four main runs: **PASS**. Official input tag: `unvetted`.

## Source and official acceptance

[Pinned official YAML](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem03-hertz-peak-wall-3d-si.yaml) | [Analytic checker](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp). Commit 29dd3fd7642a6fe9e8b4dd99120df4068b221250; input blob 531d86bb3e2993d869041734a6dce426e3737758. LAMMPS was not rerun locally. The YAML's same-code epsilon is not a cross-software tolerance.

Model `hertz_peak` is evaluated at original step 2138 and time 2.138e-05s. Half-step runs use twice the native check step at the same physical time. The exact reference, relative/absolute error type and unchanged limit are listed separately for each quantity. Whole histories and source golden states supplement these selected-time assertions.

| Run | Quantity | Reference | Native value | Error type | Error | Limit | Result |
|---|---|---:|---:|---|---:|---:|---|
| CPU dt | Overlap at official check time (m) | 2.26806424242e-05 | 2.26806479856e-05 | relative | 2.45203744731e-07 | 0.01 | PASS |
| CPU dt | Force at official check time (N) | 384.052459495 | 384.052492736 | relative | 8.65553132983e-08 | 0.01 | PASS |
| GPU dt | Overlap at official check time (m) | 2.26806424242e-05 | 2.26806479856e-05 | relative | 2.45203744731e-07 | 0.01 | PASS |
| GPU dt | Force at official check time (N) | 384.052459495 | 384.052492736 | relative | 8.65553132983e-08 | 0.01 | PASS |
| CPU dt/2 | Overlap at official check time (m) | 2.26806424242e-05 | 2.26806453362e-05 | relative | 1.28392472616e-07 | 0.01 | PASS |
| CPU dt/2 | Force at official check time (N) | 384.052459495 | 384.052425444 | relative | 8.86616481251e-08 | 0.01 | PASS |
| GPU dt/2 | Overlap at official check time (m) | 2.26806424242e-05 | 2.26806453362e-05 | relative | 1.28392472616e-07 | 0.01 | PASS |
| GPU dt/2 | Force at official check time (N) | 384.052459495 | 384.052425444 | relative | 8.86616481251e-08 | 0.01 | PASS |

## Inputs and constitutive mapping

- sphere: R=0.0025 m; position=(0, 0, 0.00255) m; velocity=(0, 0, -3.9) m/s; spin=(0, 0, 0) rad/s.

- sphere: mass=0.000458148928649 kg; isotropic inertia=1.14537232162e-09 kg m^2.

Hertz Fn_elastic=(4/3)*E_eff*sqrt(R_eff)*delta^1.5; native E_eff=5.3333325e+10 Pa. Stored linear kn is inactive. Automatic damping uses the harmonic material restitution and the active normal/tangential stiffness; stored manual dashpots are inactive. Tangential inputs: kt=0 N/m (explicit in linear mode), mu=0. Rolling/twist coefficients: kr=0, mu_r=0, ktwist=0, mu_twist=0.

Gravity=(0, 0, 0)m/s^2. Base/half dt=1e-08/5e-09s; base/half numerical steps=3038/6076. Initial poses, velocities and physical source coefficients are preserved. The parameter mapping translates the source model's coefficient definitions; it does not fit measured native results.

- Direct Hertz kn -> E*=3kn/4, both native materials E=2(1-nu^2)E*. Wall reduced radius remains R.

Native specification (retained locally; not included) | Parameter mapping (retained locally; not included) | Configuration round-trip (retained locally; not included)

## CPU/GPU and timestep refinement

| Quantity | CPU dt error | CPU dt/2 error | Error type and limit |
|---|---:|---:|---|
| Overlap at official check time (m) | 2.452e-07 | 1.2839e-07 | relative: 0.01 |
| Force at official check time (N) | 8.6555e-08 | 8.8662e-08 | relative: 0.01 |

CPU/GPU whole-history maximum position / velocity / spin / force differences: 0 m / 8.88e-16 m/s / 0 rad/s / 0 N.

Maximum individual impulse / position-drift / velocity-kick residuals: 1.6e-17 kg m/s / 2.17e-19 m / 4.38e-16 m/s. Gravity is added separately to saved native contact forces.

Energy histories are supplemental unless selected by the official analytic model. Wall/gravity cases exchange external momentum; total sphere momentum is not assumed conserved. The step-refinement table reports each error directly, without imposing monotonic convergence.

Numerical criteria and measurements are retained below and in the companion CSV tables.

## Complete retained results

All 18,232 native output frames, every integration-step particle/contact CSV, projects, raw initial samples and backend logs are retained. Native output cadence is distinct from integration dt; exact counts/intervals are in the analysis and archive manifests. The four primary frame counts are [3039, 3039, 6077, 6077]. 0additional diagnostic run is retained.

- **CPU dt** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 3,038 integration steps, 3,039 stored frames.
- **GPU dt** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 3,038 integration steps, 3,039 stored frames.
- **CPU dt/2** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 6,076 integration steps, 6,077 stored frames.
- **GPU dt/2** (main CPU/GPU/time-step comparison): loadable project (retained locally; not included), complete native result (retained locally; not included); 6,076 integration steps, 6,077 stored frames.

Analysis JSON (retained locally; not included) | Official metrics CSV (retained locally; not included) | Source checkpoints (retained locally; not included) | Independent physics audit (retained locally; not included) | Run manifest (retained locally; not included) | Full native ZIP (retained locally; not included) | Frame/hash audit (retained locally; not included)

Source checkpoint states are supplemental regression diagnostics. Their differences are reported directly; the YAML same-code epsilon is not imposed as an independent cross-software acceptance limit.

A stable helper executable is statically linked to the current native CPU/CUDA libraries. Historical GUI/library byte identity is not asserted. Existing case reports and records are unchanged. Verified scope follows the active source model and the criteria above; a passing inactive channel does not validate every option of its contact family.

## Repository evidence and native rerun input

- [Dedicated C++ case source](dem03-hertz-peak-wall.cpp). The category README documents the runner and model-mapping limits.
- [Compact historical numerical summary](HISTORICAL-RESULTS.md).
- [Source pins, original outcome and retained-artifact hashes](PROVENANCE.md).
- [Compact checkpoint table: official_metrics_comparison.csv](checkpoints/official_metrics_comparison.csv).
- [Compact checkpoint table: source_checkpoint_comparison.csv](checkpoints/source_checkpoint_comparison.csv).

A new runner result must be compared independently with these targets; opening or completing the input is not an automatic PASS. Full native frames, per-step histories, execution logs, vendor inputs, proprietary programs and PDF reports are not included.
