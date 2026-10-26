# `src/execution/cuda/`

> **Purpose:** Provide GPU kernels and host launch interfaces that operate only on device views.

## Core Contents

| File | Responsibility |
| --- | --- |
| `atomicKernel.cuh` | Atomically accumulate `Real` values and vector components |
| `integrationKernel.*` | Integrate linear/angular velocity and position/orientation |
| `spatialGridKernel.*` | Build hashes, sorted indices, and cell begin/end ranges |
| `contactDetectionKernel.*` | Detect neighbors and contacts for three particle combinations and write prefix results |
| `contactForceKernel.*` | Evaluate contact force/torque and atomically accumulate particle loads |
| `bondForceKernel.*` | Refresh bond copies, evaluate damage, and accumulate force/torque |
| `sphInteractionKernel.*` | Query SPH background grids and evaluate density, fluid–fluid, and fluid–wall interactions |
| `sphJetKernel.*` | Impose inlet velocity/reference density on active layers of a short finite-duration jet buffer |

## Key Rules

- `.cu/.cuh` filenames and `__global__` function names end in `Kernel`.
- Public launch functions do not use the `Kernel` suffix.
- Every launch interface accepts `cudaStream_t` explicitly and performs no device-wide synchronization; upper layers synchronize after a complete copy chain when host results are required.
- Use atomic accumulation only when multiple threads write the same particle force or torque. Private contact, bond, and SPH state is written directly to its own slot.
- Declare public lightweight device views in `.cuh` files; views private to one implementation may remain in a `.cu` file.
- Containers and solvers own device allocation, host/device copies, and solve-stage transitions.

Contact detection exposes three final launch flows: sphere–sphere, sphere–LS, and LS–LS. LS–LS includes a separate area-weighted effective-mass stage. CUDA SPH uses a background grid and rebuilds it only once per SPH advection step.

## Extension Guide

- Let launchers handle zero-length work, and retain valid-count bounds checks inside kernels.
- Keep sorted-slot indices separate from original particle indices; particle reads use `particleIndex_[sortedIndex]`.
- Never retain host pointers in kernels or assume identical Host AoS and Device SoA memory layouts.
- Do not hide backend checks inside launchers; calling a launcher means the upper layer has already selected the GPU path.
- After changing a host `DeviceLayout`, update the corresponding device views, launchers, and round-trip tests.

## Kernel and Launcher Structure

```text
solver-selected GPU stage
    -> launchFunction(..., stream)
        -> zero-length guard and launch dimensions
        -> XxxKernel<<<grid, block, 0, stream>>>()
        -> cudaGetLastError validation
```

Launchers enqueue work and return without synchronizing. A following kernel on the same stream observes prior writes automatically. Synchronization is required only before the CPU consumes results or when a synchronous public copy promises completion.

## Spatial-Grid Data Flow

1. Calculate one hash and original particle index per particle.
2. Sort hash/index pairs by hash while leaving the particle container untouched.
3. Clear cell ranges.
4. Scan sorted hashes to write `sortedParticleIndexBegin_` and `sortedParticleIndexEnd_`.
5. During search, map `sortedIndex` back through `particleIndex_[sortedIndex]` before reading particle fields.

The multilevel rigid grid reserves a final level for infinite-mass particles. SPH uses uniform fluid and virtual-particle grids; their interaction kernels query cells directly rather than materializing a device neighbor list.

## Main Launch Pipelines

### DEM contact

```text
build grid -> count candidates -> prefix sum -> grow contact storage
           -> write contacts -> LS effective-mass correction when required
           -> contact force -> atomic particle accumulation
```

### SPH advection/acoustic cycle

```text
each acoustic step:
    active-jet velocity
    -> [advection boundary only: build grids -> free surface -> density reinitialization -> viscous prior force]
    -> density half step -> position half step -> pressure force
    -> velocity step -> active-jet velocity -> position half step -> density half step
```

Before each density or pressure stage, the solver checks actual fluid/boundary displacement against the reference positions owned by `interaction/SPHNeighborhood`. It rebuilds the grids if those references are invalid or the maximum displacement exceeds half the search skin. A rebuild within an advection interval preserves density, viscous prior force, and the current advection phase.

Wall reaction uses atomic additions because many fluid particles may contribute to one virtual boundary point. Each fluid thread writes only its own density, pressure, prior force, and total force.

Finite-duration jets initially create only a short upstream buffer. Later source layers append actual particles at common acoustic starts; device capacity growth preserves the already integrated prefix and uploads only new tail fields on the solver stream. The particle grid is resized and neighbor references invalidated before the new layer is used. Sources retain stable append-only membership, so layers from different sources can be interleaved in the global container.

`launchApplySPHJetVelocity()` receives the stable `[begin, begin + count)` range of one currently eligible inlet layer and cannot modify unrelated particles. The kernel combines the source's timed activity with its short-pipe geometry, writes the `constrained` marker, and imposes velocity, reference density, zero equation-of-state pressure, and zero density rate while active. Density reinitialization and continuity integration skip marked particles. Constraint activity is evaluated at the acoustic start, after velocity integration at the midpoint, and at the endpoint; it is not frozen to one start-time decision for the whole interval. Released layers form a permanent prefix and are not visited by subsequent inlet launches. Once every source ends, the completion fast path suppresses later launches. The four-argument `addSPHJet` overload remains a static downstream column and never launches this kernel.

Inlet openings affect only virtual-wall quadrature. The mask launcher uses the finite upstream buffer cylinder plus half-spacing clearance, setting volume to zero inside and restoring nominal volume outside. Changed samples have force, viscous prior force, and cached coupling torque cleared; the solver removes the old owner reaction before this launch and refreshes advection preparation when the mask changes. No signed-distance grid, triangle surface, or downstream jet path is altered. GPU/Hybrid update moving sample positions before evaluating the mask on the same solver stream; output uses a temporary host mask without changing device quadrature.

## CUDA Error and Stream Rules

- Every function that reads or writes device data accepts `cudaStream_t`.
- Use `host_device_detail::checkCuda()` for runtime API errors and validate every kernel launch.
- Never call `cudaDeviceSynchronize()` inside a launcher.
- Use asynchronous field copies when only force, torque, or another subset is needed.
- Do not use a pointer from a host object as though it addressed device memory.
- Rebuild a device view after any operation that can reallocate its owning container.
