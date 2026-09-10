# `execution/cuda/`

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
| `sphJetKernel.*` | Restore prescribed velocity inside one finite-duration jet's particle range and upstream pipe |

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

Wall reaction uses atomic additions because many fluid particles may contribute to one virtual boundary point. Each fluid thread writes only its own density, pressure, prior force, and total force.

Finite-duration jet particles are allocated once on the host before device initialization. `launchApplySPHJetVelocity()` receives one source's stable `[begin, begin + count)` range and cannot modify unrelated particles; its kernel also checks that each ranged particle is still inside the upstream pipe and updates the particle's internal `constrained` marker. Density reinitialization and continuity integration kernels skip marked particles. Both applications use the solver stream and the same substep-start activity decision. The marker is cleared when a particle leaves the pipe or when the source ends. After the final source ends, the `jetsCompleted` fast path suppresses all later jet launches. The legacy four-argument `addSPHJet` overload is only a static downstream column and never launches this kernel.

## CUDA Error and Stream Rules

- Every function that reads or writes device data accepts `cudaStream_t`.
- Use `host_device_detail::checkCuda()` for runtime API errors and validate every kernel launch.
- Never call `cudaDeviceSynchronize()` inside a launcher.
- Use asynchronous field copies when only force, torque, or another subset is needed.
- Do not use a pointer from a host object as though it addressed device memory.
- Rebuild a device view after any operation that can reallocate its owning container.
