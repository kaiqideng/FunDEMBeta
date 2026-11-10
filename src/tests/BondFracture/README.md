# Bond fracture regression notes

## Current implementation and regression scope

This guide describes the current regression drivers. Recorded results under
the historical headings apply only to their own source snapshots. Current
mechanics are documented in the
[Bond user guide](../../../docs/UserGuide.md) and
[interaction overview](../../interaction/README.md).

The current driver uses `dn=2*UN/(A*openingStress)`,
`ds=2*US/(A*shearStress)` with positive-denominator guards, ordinary
`d=hypot(dn,ds)`, and `df=2*Gc/hypot(openingStress,shearStress)` frozen at
strength onset. It also clips the evaluated torsional energy to zero when
nonpositive. Retired square-root/stiffness and radial-traction laws are not
alternative implementations or current reference formulas.

The native regression sources use the current helper signatures, including
three section-stress outputs (axial, bending, shear) and no retired
displacement/stiffness inputs to `updateBondDamage()`. Their current numerical
references must follow the stress-based mapping above. A rebuilt executable
or updated reference is not itself a CPU/GPU, refinement, energy-work or
physical-validation pass; retain the actual command logs and verdicts for
the tested source. Original historical values and PASS/FAIL outcomes below
are deliberately unchanged.

## Constitutive scope

These prescribed-pose checks exercise the native CPU and CUDA kernels, not a
dynamic integrator or a calibrated rock specimen. The model is an EVM
extension inspired by Suarez et al. (2025), not an exact reproduction of that
paper's geometric interface-gap law.

The existing B3/B4 bond-point torque supplies the circular-section strength
calculation. All bending energy contributes to opening, even under axial
compression; torsion shares shear toughness. Partial damage preserves axial
compression, while complete Bond removal also releases that branch. Native
contact is not an energy-conserving handover of a removed spring. Torque
cancellation may prevent initiation despite nonzero endpoint couples.

The focused tests also cover raw-stiffness lifecycle, rigid rotation,
unloading, frozen BK mixity, all four energy contributions, compression
preservation, and zero-energy/brittle boundaries. Full definitions and
preconditions belong to the linked user guide rather than a second law here.

## Build and run

Keep the terminal in the directory containing `FunDEMBeta/` and generated
files outside source. These programs are enabled by `BUILD_TESTING`,
independently of the external-reference validation option:

~~~bash
cmake -S FunDEMBeta -B build -DBUILD_TESTING=ON -DFUNDEM_BUILD_VALIDATIONS=OFF -DFUNDEM_ENABLE_CUDA=ON
cmake --build build --target bondFractureTests fundem_bond_fracture_test --parallel
ctest --test-dir build -R '^FunDEM.BondFracture$' --output-on-failure
ctest --test-dir build -R '^FunDEM.BondFracture.Smoke.*CPU$' --output-on-failure
python3 FunDEMBeta/src/tests/BondFracture/runVerification.py \
  --build /absolute/build --output /absolute/new/result-directory
~~~

The five single-configuration executables are in
`build/src/tests/BondFracture/`; the focused unit-test executable remains in
`build/src/tests/`. The `bondFractureTests` target builds the five programs;
it does not implicitly execute a campaign.

CTest additionally registers five `FunDEM.BondFracture.Smoke.<case>.CPU`
checks. A CUDA-enabled build registers the corresponding `.GPU` checks;
those require an actual usable device and are not reported as CPU fallbacks.
GPU smoke checks serialize access to the shared CUDA device. Each invocation
uses a new timestamp/random-suffix directory in the build's
`src/tests/BondFracture/ctest/`, so rerunning CTest never overwrites evidence.
The shared runner checks loading completion, non-healing finite damage,
immediate failed-load clearing, and records actual particle work separately
from radial traction work. Only pure axial opening requires `work=A*Gc`;
other modes retain that deviation as a diagnostic. These are
lightweight single-increment constitutive smoke checks, not a convergence or
general fracture-energy claim. A native failure returns nonzero. The recorded
work checks do not establish increment convergence: the separate strict
campaign below additionally decides history/refinement acceptance.

The script evaluates five programs on actual CPU/GPU at base/half/quarter
loading increments: 30 native runs. It retains every state, command log,
English case/aggregate report, comparison, onset check, pre-run source
snapshot and hash manifest. CUDA execution requires a CUDA-capable
build/device. `--backends cpu` runs a 15-run CPU-only campaign and reports
CPU/GPU parity as `NOT RUN`, never PASS. `--binary-dir` supplies the exact
executable directory for multi-configuration builds when needed. Python uses
only its standard library and is not required for ordinary CTest smoke runs.

`--base-refinement N` changes only sampling resolution. Each case uses
`2000*N`, `4000*N` and `8000*N` intervals; its poses, stiffnesses, strengths,
toughnesses and acceptance tolerances remain unchanged. For example:

~~~bash
python3 FunDEMBeta/src/tests/BondFracture/runVerification.py \
  --build /absolute/cpu-build --backends cpu --base-refinement 8 \
  --output /absolute/new/fine-cpu-results
~~~

The strict campaign compares h with h/2, h/2 with h/4, and h with h/4 at
matching prescribed poses, with the same original numerical tolerances.
Energy or refinement FAIL makes the campaign exit nonzero; a missing/broken
native run also fails with its evidence retained. Output must be a fresh
directory: old PASS and FAIL campaigns remain unchanged. No tolerance is
relaxed and no model or fixture parameter is tuned to make these checks pass.
Git HEAD, tracked-tree status and exact source SHA-256 values identify the
tested implementation; no hard-coded checkpoint is used. Changing the
audited sources while the campaign runs makes its final verdict FAIL.

| Program/source | Prescribed motion |
| --- | --- |
| bond_fracture_tension / tension.cpp | Axial opening 0 to 0.012 m |
| bond_fracture_shear / shear.cpp | Unit-radius center arc 0 to 0.06 rad |
| bond_fracture_bending / bending.cpp | Opposite Z rotations to +/-0.05 rad, kb=200 |
| bond_fracture_torsion / torsion.cpp | Slave X twist 0 to 0.05 rad |
| bond_fracture_bending_zero_torque / bendingZeroTorque.cpp | Original kb=100 bend with native torque cancellation |

Options are --cpu, --gpu, --refinement N, --half-step, --quarter-step, --smoke and
--output NEW_DIRECTORY. `--smoke` labels the single-increment scope without
changing fixture parameters or force calculations.
Default parameter-time increments are 0.0005, 0.00025 and 0.000125 s
(2000/4000/8000 intervals).
These are constitutive-kernel evaluations at prescribed poses, not dynamic
integration.

## Work and verification scope

True mechanical work integrates particle forces and endpoint torques against
prescribed motion. It must balance independently integrated damage
dissipation plus retained elastic energy. Radial scalar traction work
integrates `(1-D)*2*(UN+US)/(A*d)` against `d` and is recorded separately.
Only pure axial opening requires mechanical work `A*Gc` and scalar area `Gc`.
The current section-stress-based failure length does not impose that budget
for general loading modes, even at fixed modal direction. Intact zero-torque
bending instead checks external work against stored elastic energy.

The Python onset audit independently reconstructs section stresses and
component separations from the first native crossing. It checks frozen
onset values and irreversible damage history. Updating an oracle does not
make a failed refinement or energy comparison pass.

Original limits remain 2e-9 normalized CPU/GPU error, 0.005 loading-increment
error, 1e-10 absolute noise floor and 2e-7 diagnostic-area/applicable-work
error. No failure is hidden by changing tolerances.
Three increments provide numerical consistency/refinement evidence, not an
assumed convergence order, physical validation, or a calibrated rock model.

## Fixture parameters and historical results

L=A=1 in SI units; kn=1000 N/m; ks=400 N/m; kt=100 N m/rad.
Tensile strength=cohesion=1 Pa; friction=0; GIc=GIIc=0.005 J/m^2; eta=1.75.
Kb=100 except the separate nondegenerate bending-fracture case (kb=200).

For the original kb=100 bending fixture, B3+B4/2=kb-ks*L^2/4=0. Its native
TBT cancels despite nonzero endpoint couples and stored energy, so it stays
elastic under the retained strength convention. These are verification
values, not calibrated rock properties.

Extreme-weight checks include subnormal values, overflow-prone sums and
high-stiffness low-energy onset. CPU/GPU execution and actual work remain
separate observations; read the executed test output before claiming a pass.

Earlier onset-energy, stress-inverse, and four-component separation campaigns
remain historical records. Some earlier campaigns used subsequently rejected
force/energy-function changes; others tested restored original expressions
under the preceding damage law. Preserve each campaign's actual PASS/FAIL
status, source snapshot, hashes and scope. None is relabelled as an execution
of the current stress/energy-based separation law. New verification requires a fresh
output directory and its own source provenance.

### Historical Euclidean-coordinate checks on 2026-11-03

- All 19 selected lightweight CPU CTests passed. The CUDA focused test
  passed, including seven native CPU/GPU loading paths and the added
  constant-ratio tension/shear path. All ten CPU/GPU constitutive smoke
  checks also passed. These focused/smoke results are not a full strict
  CPU/GPU refinement campaign.
- Native work-balance test oracles use 8,000 and 16,000 loading intervals,
  without changing physical parameters or the `2e-7` work tolerance.
  For the nonproportional mixed path, trapezoidal `integral(Ua*dD)` error
  was about `1.039996e-6 J` at 4,000 intervals, `2.59999e-7 J` at 8,000,
  and `6.49998e-8 J` at 16,000. The coarser estimates remain failed at
  the unchanged tolerance; the fine oracle is not evidence that all
  loading-increment comparisons have converged.
- That CPU-only h/h2/h4 campaign completed 15 native runs and
  70,015 states. Native energy-balance, onset/history and applicable
  constant-ratio work checks passed. Loading-increment refinement failed,
  so the overall campaign verdict remains **FAIL**. Full-campaign
  CPU/GPU parity was **NOT RUN** for this Euclidean-coordinate campaign;
  it must not be inferred from the focused or smoke passes.
- These are constitutive implementation checks, not physical specimen
  validation. That Euclidean and preceding weighted-coordinate
  evidence remain separate, each with its own source scope and verdict.

### Historical weighted-coordinate checks on 2026-11-03

The following retained results tested the preceding stiffness-weighted
coordinate, not the later Euclidean separation and frozen radial-traction
definition or the live stress/energy-based law. Their original source snapshots, metrics and verdicts remain
unchanged. They must not be quoted as executed PASS results for the current law.

- All 19 selected lightweight CPU CTests passed; the long Fiber integration
  test was not run. The CUDA focused test, including native CPU/GPU loading,
  unloading and compression-retention comparisons, passed.
- The default CPU/GPU campaign completed 30 native runs and 140,030 states.
  Pure-mode mechanical fracture work, conjugate scalar work, onset/history
  audits and CPU/GPU parity passed. The largest damage difference between
  CPU/GPU was about `5.62e-11`; the largest pure-mode work-budget error was
  about `1.31e-8 J`. Loading-increment refinement failed, so its overall
  verdict remains **FAIL**.
- A separate CPU-only campaign with base refinement 8 completed 15 runs and
  560,015 states. Mechanical work and onset/history audits passed. The only
  remaining failed comparisons were the small shear-energy component in
  bending: maximum absolute difference about `2.25e-10 J`, exceeding the
  unchanged absolute/normalized refinement criteria. Its overall verdict
  also remains **FAIL**; CPU/GPU parity was **NOT RUN** in that campaign.
- Both campaigns retained unchanged source hashes during execution and
  preserve their own source snapshots. Neither establishes specimen-scale
  physical validation or permits a blanket all-checks-passed claim.

## Reference

Suarez et al. (2025), *Multiscale modeling of rock fracture in comminution —
A comparative study of FEM accuracy and DEM scalability*, Minerals Engineering
232, 109488. https://doi.org/10.1016/j.mineng.2025.109488.
