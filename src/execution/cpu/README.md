# `src/execution/cpu/`

> **Purpose:** Execute generic integration and interaction-force assembly over host containers for stages selected by the solver.

## Core Contents

| File | Responsibility |
| --- | --- |
| `particleFunctions.h` | Clear particle loads, integrate rigid motion, build temporary interaction adjacency, and assemble contact/bond forces and torques |
| `parallelPolicy.h` | Define OpenMP thresholds of 1024 particles, 256 interactions, and 4096 candidate master-surface-node queries; LS query blocks contain at most 1024 nodes |
| `pointGrid.*` | Share sorted occupied-cell indexing between contact broad phase and exact fixed-radius point-pair searches |

The container assembly helpers are header-only templates in `fundem::cpu`; `pointGrid.cpp` is compiled into `FunDEM::execution`. OpenMP parallelizes suitable loops when enabled. For larger workloads, contact and bond assembly first evaluates individual interactions and then gathers loads through temporary particle-to-interaction adjacency, so workers do not concurrently accumulate into the same particle. Small interaction sets use direct serial accumulation.

CPU LS–LS detection in `interaction/contactSearch` uses the total candidate master-surface-node work, rather than pair count alone, to select parallel execution. Ordered pair/node ranges are flattened into blocks of at most 1024 nodes, allowing even a few high-resolution pairs to share the work at the 4096-query threshold. Parallel stages use one OpenMP level, without nested pair/node regions. Contact construction reuses block buffers; final merging preserves pair and surface-node order. See the [interaction guide](../../interaction/README.md#cpu-search-flow) for area normalization and history rules.

## Fixed-Radius Point Queries

`pointGrid::build(points, cellSize)` indexes finite host positions without owning
or reordering the input array. `forEachCandidate(position, radius, visitor)` is a
read-only conservative query: the caller still performs its exact geometry test.
It visits cells in z/y/x order and original indices in increasing order within
each cell, preserving the contact broad phase's traversal. Large sparse ranges
scan occupied cells rather than allocating or iterating a dense empty domain.

`findPointPairs(points, radius)` returns self pairs with `firstIndex_ <
secondIndex_`. Its three-argument overload accepts separate first/second arrays
and returns local indices into those arrays. Both use the inclusive exact test
`std::hypot(dx, dy, dz) <= radius` and sort by first/second index independently of
OpenMP scheduling. They have no particle radius, mass, material, or contact/SDF
semantics: separated particles may form a distance-based Bond if their centers
satisfy the radius criterion. Contact search retains its separate bounding-sphere
and narrow-phase tests.

Query radii must be finite and nonnegative; zero finds distinct coincident points.
Cell size must be finite and positive. Invalid inputs are rejected before pair
workers run. Queries union their original cell neighborhood with an outward-rounded
world-space box, mapped through the same quantization as the index. This protects
nearby points from large-origin subtraction and cell-boundary rounding. The grid
checks signed 64-bit coordinate conversion and conservatively falls back to all
candidates if a coordinate or query extent is unrepresentable; the final exact
test still decides pair inclusion. Inputs cannot exceed the
`int` local-index representation. Concurrent const queries share one built grid
without nested parallel regions; rebuilding or clearing it requires exclusive use.

```cpp
#include "execution/cpu/pointGrid.h"

const std::vector<fundem::math::Vec3> centers{{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}};
const auto pairs = fundem::cpu::findPointPairs(centers, 0.75);
// pairs[0] refers to original centers[0] and centers[1]; no contact is implied.
```

The point grid is nonperiodic and CPU-only. It neither creates periodic images
nor applies `math::periodicDomain`; an external boundary policy must provide any
additional image queries. It is separate from the CUDA particle `spatialGrid`
and from the SPH engine's neighborhood cache. Sharing this broad-phase index does
not replace the existing sphere/LS narrow-phase or contact history.

`FunDEM.PointGrid` compares self/cross results with exact brute force, including
threshold endpoints, coincident points, sparse scenes and extreme coordinates.
It also checks all three CPU contact families against direct search and verifies
stable history and deterministic ordering across OpenMP thread counts.

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
