# `src/execution/cpu/`

> **Purpose:** Execute generic integration and interaction-force assembly over host containers for stages selected by the solver.

## Core Contents

| File | Responsibility |
| --- | --- |
| `particleFunctions.h` | Clear particle loads, integrate rigid motion, build temporary interaction adjacency, and assemble contact/bond forces and torques |
| `parallelPolicy.h` | Define OpenMP thresholds of 1024 particles, 256 interactions, and 4096 candidate master-surface-node queries; LS query blocks contain at most 1024 nodes |

The helpers are header-only templates in `fundem::cpu`. OpenMP parallelizes suitable loops when enabled by the consuming target. For larger workloads, contact and bond assembly first evaluates individual interactions and then gathers loads through temporary particle-to-interaction adjacency, so workers do not concurrently accumulate into the same particle. Small interaction sets use direct serial accumulation.

CPU LS–LS detection in `interaction/contactSearch` uses the total candidate master-surface-node work, rather than pair count alone, to select parallel execution. Ordered pair/node ranges are flattened into blocks of at most 1024 nodes, allowing even a few high-resolution pairs to share the work at the 4096-query threshold. Parallel stages use one OpenMP level, without nested pair/node regions. Contact construction reuses block buffers; final merging preserves pair and surface-node order. See the [interaction guide](../../interaction/README.md#cpu-search-flow) for area normalization and history rules.

## Layer Boundary

Shared numerical formulas at the root of `execution/` operate on values and are usable by CPU and CUDA. The helpers here iterate over host containers, invoke particle and interaction methods, and may allocate temporary assembly storage. `execution/cuda/` provides the corresponding device kernels and launchers.

The CPU search caches pair rotations and relative positions for one search, then samples the slave SDF in local coordinates through the shared `sampleLevelSetContact` helper. This runtime CPU cache does not add device-layout state; the shared world-space contact-query wrapper used by CUDA is preserved.

Persistent contact histories, SPH neighbor lists, displacement references, and virtual-boundary coupling accumulators belong to `interaction/`. Its CPU SPH engine applies shared formulas to its own cached neighborhoods; it includes `parallelPolicy.h` for the common loop threshold without depending on generic rigid-particle assembly helpers.

Solvers select execution modes and own initialization, step order, simulation clocks, output, and snapshots. These CPU helpers perform only the requested operation and retain no solver lifecycle state. Applications normally use `FunDEM::FunDEM`, which supplies the particle and interaction implementations needed by these templates.

## Extension Guide

- Keep new shared equations in the root numerical headers.
- Keep persistent search or coupling state in `interaction/`.
- Preserve stable particle indices and race-free force/torque assembly.
- Keep stage ordering and snapshot decisions in the solver.
