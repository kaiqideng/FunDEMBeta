# `src/tests/`

> **Purpose:** Verify stable infrastructure, key physical invariants, and host/device consistency.

## Core Contents

| Test | File | Coverage |
| --- | --- | --- |
| `FunDEM.Physics` | `PhysicsTest.cpp` | Mass/inertia, contacts, bonds, materials, geometry, execution modes, SPH wall-reaction conservation, and important degenerate states |
| `FunDEM.MaterialPair` | `MaterialPairTest.cpp` | Four-property bulk materials, two group switches, explicit zero overrides, LS unset sentinels, stable symmetric pair indices, append-only growth, and move semantics |
| `FunDEM.ContactPairResponse` | `ContactPairResponseTest.cpp` | Hertz–Mindlin defaults and energy, final pair overrides, group-controlled dashpots, LS nodal damping weights, and missing LS-pair configuration |
| `FunDEM.RandomShape` | `RandomShapeTest.cpp` | Seeded radial deformation of sphere/superellipsoid bases, finite mesh distances, outward surfaces, input guards, and nonaccumulating centroid correction |
| `FunDEM.HostContainerBinding` | `HostContainerBindingTest.cpp` | Stable ownership/index binding through append, vector relocation, copy/move and assignment, plus concurrent const reads |
| `FunDEM.LSGeometry` | `LSGeometryTest.cpp` | Centroid frame, grid/query agreement, cached volume/inertia, fixed-wall defaults, geometry insertion without a second correction, and regeneration invalidation |
| `FunDEM.LSPrimitive` | `LSPrimitiveTest.cpp` | Finite primitive topology, outward watertight surfaces and edge-spacing targets, exact distance cases, cavity preservation, integrated volume, local-frame bounds, rebuilds, and fixed policy |
| `FunDEM.PointGrid` | `PointGridTest.cpp` | Exact self/cross radius pairs against brute force, inclusive thresholds, stable original indices, sparse/extreme-coordinate cases, thread-count determinism, and all CPU contact-family regressions |
| `FunDEM.LSContactParallel` | `LSContactParallelTest.cpp` | Serial/parallel LS contact agreement, weighted pair masses, cache lifecycle, and stable node/history behavior |
| `FunDEM.ContactMaterial.GPU` | `ContactMaterialTest.cu` | 26 CPU/CUDA contact fixtures covering sphere–sphere, sphere–LS, LS–LS, infinite mass, both group switches, pure dampers, history, energy, and force/torque scattering on a nondefault stream |
| `FunDEM.Fiber` | `FiberTest.cpp` | Eleven-sphere bonded cantilever with one fixed end, zero gravity, a ramped 100 kN endpoint load, 0.001 global damping, and a `1e-9 m/s` maximum particle-speed limit |
| `FunDEM.Math` | `MathTest.cu` | Vectors, matrices, quaternions, and host/device mathematical consistency |
| `FunDEM.HostAoSDeviceSoA` | `HostAoSDeviceSoATest.cu` | Device layout, field access, size/copy semantics, and suffix append preserving an evolving Device prefix across growth |
| `FunDEM.SPHState.CPU` | `SPHCPUStateTest.cpp` | CPU SPH neighborhoods, leading multirate state, current-time output, continuation, and stability limits |
| `FunDEM.SPHCoupling` | `SPHCouplingTest.cpp` | Held-load DKD/KDK momentum budgets, pressure/viscous reaction, rotating solids, dense output, nested/throwing observations, public host copies, and split continuation on CPU/GPU/Hybrid |
| `FunDEM.SPHJet` | `SPHJetTest.cpp` | True-inlet buffer generation and growing ownership, aperture masks, source discretization, one-way particle release, completion, multiple jets, time-step changes, snapshots and continuation |
| `FunDEM.SPHRobustness` | `SPHRobustnessTest.cpp` | Invalid computed-density diagnostics, search-skin expiry, safe source discretization, and integer-aligned acoustic-step limits |
| `FunDEM.SPHState` | `SPHStateTest.cu` | Hybrid CUDA SPH kernels, leading multirate state, output snapshots, host/device consistency, and CPU/Hybrid numerical agreement |
| `FunDEM.SPHRobustness.GPU` | `SPHRobustnessTest.cu` | Device displacement guards, CPU/GPU invalid-density diagnostics, and approaching-pair recovery during a reused advection interval in GPU/Hybrid modes |

## Key Rules

- Pure C++ unit tests, including `FunDEM.MaterialPair` and `FunDEM.ContactPairResponse`, and the CPU integration test `FunDEM.Fiber` are always built; CUDA-specific tests are added when CUDA is enabled.
- All test targets enable the project's standard compiler warnings.
- CTest labels use `unit;cpu`, `unit;cuda`, or `integration;cpu`.
- Long physical validation belongs in an explicitly labeled integration test or the repository-root [`example/`](../../example/README.md), not in the unit-test suite.

When the ON-DEM validation check targets are built, CMake additionally registers
`FunDEM.ONDEM.PeriodicShear`, `FunDEM.ONDEM.ShearIO`, and
`FunDEM.ONDEM.ExternalContactForce` with the `validation;cpu` labels. Their
presence is conditional; the standard unit suite does not require downloaded
benchmark coordinates. See the [validation guide](../../validation/README.md)
for benchmark/reference inputs and reporting scope.

## Build and Run

Open a terminal in the parent `workplace/` directory containing `FunDEMBeta/` and keep it there for all commands:

```bash
cmake -S FunDEMBeta -B build -DBUILD_TESTING=ON -DFUNDEM_ENABLE_CUDA=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The CMake source directory remains the repository root. Test sources live in `FunDEMBeta/src/tests/`, and single-configuration builds place their executables in `build/src/tests/`. For example, run `build/src/tests/fundem_physics_test` directly when diagnosing `FunDEM.Physics`.

Run only one backend label:

```bash
ctest --test-dir build -L cpu --output-on-failure
ctest --test-dir build -L cuda --output-on-failure
```

## Extension Guide

- Use `fundem_add_cpp_test()` for pure C++ tests.
- Use `fundem_add_cuda_test()` for tests that require kernels or device fields.
- Prefer formula invariants, index bounds, valid sizes after copies, and CPU/GPU consistency over implementation-specific behavior.

## Test Matrix

| Configuration | Always-available coverage | Additional CUDA coverage |
| --- | --- | --- |
| CPU-only | Physics and contact formulas, material pairs, point-grid search, ownership, LS geometry/primitives, CPU SPH state/coupling/jets/robustness, fiber integration | None |
| CUDA-enabled | All C++ tests; guarded C++ fixtures also exercise available GPU/Hybrid paths | Math device parity, container round trips, contact-pair response parity, Hybrid SPH state, device search/density robustness |

`FunDEM.Physics` includes small system-level checks such as backend tags, ownership rejection, geometry-driven contact-family selection, explicit pair response, wall acceleration, CPU parallel assembly, and fluid-wall reaction conservation in every available execution mode. `FunDEM.MaterialPair` verifies that disabling either group retains its values, explicit zero differs from a default, and LS stiffness configuration does not depend on a switch. `FunDEM.ContactPairResponse` independently checks the analytic Hertz–Mindlin formulas, normal-energy derivative, and conservation of pair damping under LS surface-node weighting. `FunDEM.ContactMaterial.GPU` compares the shared constitutive response and particle loads against CPU evaluation, including the two group switches and zero-spring damping. `FunDEM.SPHState.CPU` verifies that output cadence and split solves do not change CPU SPH state. The CUDA SPH test compares the same state against Hybrid execution.

## Labels and Intended Runtime

- `unit;cpu`: fast host tests suitable for every edit.
- `unit;cuda`: fast device tests requiring a CUDA-capable runtime.
- `integration;cpu`: longer physical evolution such as the fiber convergence test.
- `validation;cpu`: conditional external-reference benchmark checks.

Labels classify test targets, not every branch in a test: C++ tests such as
`FunDEM.SPHCoupling` and `FunDEM.SPHJet` also execute guarded GPU/Hybrid cases
when compiled with CUDA. A registered test or coverage description is not a
claim that it has passed on the current machine; retain the actual CTest output
when reporting verification.

Use `ctest --test-dir build -L unit` during normal development. Run integration tests before publishing changes to time integration, force assembly, damping, bonds, or convergence criteria. Examples and tutorials are not registered as tests because they are user workflows and may produce substantial output.

## Writing a Focused Test

1. Test an observable invariant rather than a private implementation detail.
2. Use the smallest particle count that activates the required path.
3. Check both zero/degenerate and normal inputs where applicable.
4. Avoid output unless output lifecycle is the subject of the test.
5. For CPU/GPU parity, construct identical solvers and compare host-visible state after synchronization.
6. Remove temporary output directories even when a test creates multiple frames.

Useful invariants include action-reaction balance, momentum conservation, finite state reductions, stable host/device sizes, quaternion normalization, infinite-mass immobility, preserved histories, and identical results across output schedules.

## Failure Diagnosis

Run the failing executable directly when CTest output is insufficient. For CUDA failures, confirm the configured architecture and current device before changing numerical tolerances. A CPU/GPU mismatch should first be reduced to one interaction or a few SPH particles; do not hide a systematic difference by widening tolerance globally.
