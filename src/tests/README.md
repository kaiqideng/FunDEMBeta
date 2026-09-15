# `src/tests/`

> **Purpose:** Verify stable infrastructure, key physical invariants, and host/device consistency.

## Core Contents

| Test | File | Coverage |
| --- | --- | --- |
| `FunDEM.Physics` | `PhysicsTest.cpp` | Mass/inertia, contacts, bonds, materials, geometry, execution modes, SPH wall-reaction conservation, and important degenerate states |
| `FunDEM.Fiber` | `FiberTest.cpp` | Eleven-sphere bonded cantilever with one fixed end, zero gravity, a ramped 100 kN endpoint load, 0.001 global damping, and a `1e-9 m/s` maximum particle-speed limit |
| `FunDEM.Math` | `MathTest.cu` | Vectors, matrices, quaternions, and host/device mathematical consistency |
| `FunDEM.HostAoSDeviceSoA` | `HostAoSDeviceSoATest.cu` | Automatic device layout, field access, size semantics, and bidirectional copies |
| `FunDEM.SPHState.CPU` | `SPHCPUStateTest.cpp` | CPU SPH neighborhoods, deferred time-step state, output snapshots, continuation, and stability limits |
| `FunDEM.SPHJet` | `SPHJetTest.cpp` | Jet value validation, virtual-pipe geometry, source discretization, generated particles, inlet enforcement, and completion state |
| `FunDEM.SPHState` | `SPHStateTest.cu` | Hybrid CUDA SPH kernels, deferred time-step state, output snapshots, host/device consistency, and CPU/Hybrid numerical agreement |

## Key Rules

- `FunDEM.Physics`, `FunDEM.SPHState.CPU`, and the CPU integration test `FunDEM.Fiber` are always built; CUDA-specific tests are added when CUDA is enabled.
- All test targets enable the project's standard compiler warnings.
- CTest labels use `unit;cpu`, `unit;cuda`, or `integration;cpu`.
- Long physical validation belongs in an explicitly labeled integration test or the repository-root [`example/`](../../example/README.md), not in the unit-test suite.

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
| CPU-only | Physics invariants, CPU SPH state, fiber integration | None |
| CUDA-enabled | All CPU tests | Math device parity, container round trips, Hybrid SPH state |

`FunDEM.Physics` includes small system-level checks such as backend tags, ownership rejection, sphere-LS stiffness selection, wall acceleration, CPU parallel assembly, and fluid-wall reaction conservation in every available execution mode. `FunDEM.SPHState.CPU` verifies that output cadence and split solves do not change CPU SPH state. The CUDA SPH test compares the same state against Hybrid execution.

## Labels and Intended Runtime

- `unit;cpu`: fast host tests suitable for every edit.
- `unit;cuda`: fast device tests requiring a CUDA-capable runtime.
- `integration;cpu`: longer physical evolution such as the fiber convergence test.

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
