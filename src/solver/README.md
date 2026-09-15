# `src/solver/`

> **Purpose:** Coordinate physical objects, compute backends, time integration, data synchronization, and output.

## Core Contents

```text
solver
└── LSDEM
    ├── SphereDEM
    └── SPHDEM
```

- `solver`: boundaries, gravity, time, device, output schedule, and common lifecycle.
- `LSDEM`: materials, level-set geometry, LS particles, LS–LS contacts, and LS–LS bonds.
- `SphereDEM`: adds spheres, sphere–sphere interactions, and sphere–LS interactions.
- `SPHDEM`: adds CPU/GPU SPH, virtual wall particles, and SPH–LS reaction forces; ordinary spheres are disabled.
- `solvers.h`: unified user include entry point.

Reusable CPU loops live in [`execution/cpu/`](../execution/cpu/README.md). SPH neighborhoods, interaction evaluation, and wall-coupling accumulators live in [`interaction/`](../interaction/README.md), while `particle/SPHJet.h` defines the inlet source values. Solvers own these components' lifetimes and decide when each stage runs.

## Key Rules

### Backend Selection

The solver owns one `executionMode` selected before initialization:

- `LSDEM`: `CPU` or `GPU`.
- `SphereDEM`: `CPU`, `GPU`, or `Hybrid`; `Hybrid` runs spheres on the GPU and LSParticles on the CPU.
- `SPHDEM`: `CPU`, `GPU`, or `Hybrid`; `Hybrid` runs SPH on the GPU and LSParticles on the CPU.

Each solver supplies complete CPU, GPU, and Hybrid initialization and stepping templates for the modes it supports. The selected template owns the full control flow, so physical stages and CUDA launchers do not repeat backend checks. Internal overloads use `cpu::mode`, `gpu::mode`, and `hybrid::mode` instead of backend suffixes. Every helper that processes device data receives its CUDA stream explicitly. SphereDEM Hybrid transfers the required LSParticle state for GPU sphere–LS interaction. SPHDEM Hybrid transfers averaged virtual-boundary kinematics to the GPU and returns one resultant force and torque per CPU-side LSParticle only when SPH advances.

### Lifecycle

Typical call order:

1. Add materials, geometries, particles, and bonds.
2. Configure boundaries, gravity, time step, and output.
3. Select the solver execution mode.
4. Call `solve(numberOfSteps)`.

The first `solve()` initializes automatically. At completion, device state and interactions are synchronized to the host and `initialized()` returns to false. Objects may then be added, and the next solve rebuilds storage for the new container sizes. Configuration that would invalidate device state cannot change while initialized.

Solver-owned containers are exposed to users only as const references. Existing elements must never be swapped, reordered, replaced, or erased; add valid objects through the solver's `add...` interfaces so every insertion index remains stable.

Each DEM step uses a velocity-Verlet structure. Force assembly order is clear state → same-type contacts/bonds → cross-type contacts/bonds → user external fields. Protected virtual external-field interfaces receive a host container on CPU and the corresponding `device_type` on GPU.

### SPH Time Integration

`SPHDEM::setSPHProperties(spacing, smoothingLength, referenceDensity, dynamicViscosity)` configures uniform discretization and fluid properties shared by all SPH particles. Individual `SPHParticle` objects contain only dynamic state.

`addSPHJet(inletCenter, velocity, radius, length)` retains the legacy behavior of adding one static downstream particle column without an ongoing inlet constraint. `addSPHJet(SPHJet{outletCenter, direction, radius, speed, duration})` instead allocates all particles once in an upstream virtual pipe of length `speed * duration`; solving never injects particles or resizes that source. Each constraint owns one stable particle range, and only particles from that range that remain inside its pipe receive the prescribed velocity.

CPU, GPU, and Hybrid execution apply an active jet at the acoustic-substep start before optional advection preparation and again immediately after velocity integration. Particles that remain inside the virtual pipe carry an internal `constrained` marker, retain their current density and pressure during density reinitialization, and skip continuity-equation density integration. The marker clears as soon as a particle leaves the pipe or the source ends. A constraint active at the start remains active for the complete substep, so its end time is rounded upward by at most one acoustic substep. After the last jet ends, the `jetsCompleted` fast path skips all remaining range and pipe checks; `SPHJetsCompleted()` reports that completion state.

The SPH acoustic time step is an integer multiple of the DEM time step without exceeding the acoustic limit; its minimum is one DEM step. If one DEM step already exceeds the acoustic or advection/viscous limit, initialization rejects the configuration. At the start of each advection step, `SPHDEM` builds the fluid and virtual-wall neighborhood on the selected SPH backend, updates the free-surface state, reinitializes density, and stores the viscous force as a prior force. Neighbor queries include a relative-displacement buffer derived from the maximum fluid and virtual-wall velocity and acceleration. Acoustic substeps reuse the search structures while `interaction/SPHNeighborhood` confirms that actual displacement remains within half that buffer. Before each density or pressure stage, exhausted or invalid references trigger a search rebuild without reinitializing density, refreshing the viscous prior force, or restarting the advection step. Virtual-wall kinematics accumulated over multiple DEM steps are averaged by step count. Particle-count or SPH-configuration changes invalidate deferred state and rebuild it from the current host state.

### Output

After `setOutputDirectory()` and `setOutputStepInterval()`, each output frame reports the calculated step count, frame index, physical time, and device memory usage. VTU output is binary by default; use `setVTUOutputFormat(vtuFormat::ascii)` for ASCII.

- Sphere defaults: position, radius, and velocity.
- Movable LS-particle defaults: surface-node position and velocity. Fixed LS particles are written separately as `fixedLSParticle`; distance to center is optional.
- SPH defaults: position, velocity, mass, and density.
- Contact/bond defaults: point and normal.
- `energy.dat`: solid-particle kinetic, gravitational, contact, and bond energy. SPH-fluid energy is intentionally excluded.

Select additional fields through `setSphereVTUFields()`, `setLSParticleVTUFields()`, `setSPHParticleVTUFields()`, `setContactVTUFields()`, and `setBondVTUFields()`.

## Extension Guide

- Keep backend decisions at solver stages instead of repeating them inside launchers.
- Extend the physical hierarchy for new particle systems instead of creating separate CPU and GPU solver files.
- Implement user external fields through the protected virtual interfaces: host containers for CPU and `device_type` for GPU.
- Preserve force assembly order: clear → internal interactions → cross-type interactions → user external fields.

## Solver Ownership

| Solver layer | Added persistent resources |
| --- | --- |
| `solver` | Domain, gravity, base time step, execution mode, device/stream state, counters, and output configuration |
| `LSDEM` | Materials, geometry, LS particles, LS interactions, CPU search, and LS spatial grid |
| `SphereDEM` | Sphere particles, sphere interactions, sphere-LS interactions, CPU search, and sphere grid |
| `SPHDEM` | SPH particles, virtual boundary particles, SPH grids/CPU neighborhoods, and virtual-particle coupling |

The hierarchy shares construction and output interfaces through inheritance but keeps particle-specific work in the derived layer. `SphereDEM` and `SPHDEM` are alternative extensions of `LSDEM`; a single solver never contains both ordinary spheres and SPH fluid.

## Execution-Mode Data Ownership

| Solver/mode | Host-calculated state | Device-calculated state | Cross-boundary data |
| --- | --- | --- | --- |
| `LSDEM::CPU` | LS motion, contacts, bonds | None | None |
| `LSDEM::GPU` | User setup/output | LS motion, contacts, bonds | Final LS/contact/bond state |
| `SphereDEM::CPU` | All sphere and LS work | None | None |
| `SphereDEM::GPU` | User setup/output | All sphere and LS work | Final state and interactions |
| `SphereDEM::Hybrid` | LS motion and LS-LS work | Sphere motion, sphere-sphere, sphere-LS | Required LS state and resulting LS loads |
| `SPHDEM::CPU` | LS and SPH work | None | None |
| `SPHDEM::GPU` | User setup/output | LS and SPH work | Final solid/fluid state |
| `SPHDEM::Hybrid` | LS motion and LS-LS work | SPH and fluid-wall work | Averaged virtual-boundary kinematics and owner loads |

Backend tags on containers describe where calculation occurs; they do not by themselves allocate memory. Initialization performs the copies and allocations required by the selected complete solver template.

## Lifecycle State Machine

```text
construction/reconfiguration
    -> initialize on demand
    -> initialized stepping phase
    -> solve completion snapshot and host synchronization
    -> construction/reconfiguration allowed again
```

`solve(n)` owns a complete stage and returns with current host-visible state. `step()` leaves the solver initialized so repeated manual stepping avoids reinitialization. While initialized, operations that invalidate container sizes, mode selection, time integration, or device state are rejected.

Repeated `solve()` calls preserve time, total step count, frame count, valid contact history, and existing objects. Appending objects or changing prescribed LS motion invalidates deferred integration state so the next initialization rebuilds working storage from the current host objects.

## Force Assembly Order

Every DEM stage follows the same order:

```text
clear force and torque
    -> same-type contacts and bonds
    -> cross-type contacts/bonds or SPH wall reaction
    -> user-defined external force and torque
    -> second velocity half-step
```

This ordering guarantees that user hooks add to, rather than silently replace, internal physics. GPU hooks receive a device view and stream; CPU hooks receive the mutable host container selected by the execution template.

## SPH Stage Coordination

`SPHDEM` uses `interaction/sphInteraction.*` for CPU neighbor lists, WCSPH stages, integration, and reductions. That component's header remains private implementation and is not installed. `interaction/SPHNeighborhood.*` tracks search validity for both CPU and CUDA, and `interaction/virtualParticleCoupling.*` retains host boundary kinematics and resultant owner loads.

The solver owns uniform SPH configuration, time-step selection, advection/acoustic stage ordering, and snapshot restoration. It passes current values to the interaction components for each stage and decides when a neighborhood must be rebuilt. The components do not advance the solver clock or choose an output schedule.

The CPU and CUDA sequences are deliberately parallel: free-surface classification, density reinitialization, viscous prior force, two density half-steps, pressure force, velocity integration, and wall reaction use the same formulas from `execution/sphFunctions.h`.

## SPH Deferred-Time Handling

SPH acoustic time can represent several DEM steps. Between SPH updates, the solver accumulates virtual-boundary kinematics and applies the most recent owner force/torque to LS particles. Output or solve completion may occur in the middle of such an interval. A temporary snapshot flushes SPH to the exact visible time for output while retaining the unflushed continuation state used by the next solve stage. Changing configuration discards that deferred state and restarts from the visible host state.

## Extension Checklist

1. Decide which solver layer owns the new particle type or coupling.
2. Add one full stage for every supported execution mode rather than scattering backend checks.
3. Keep device helpers stream-explicit.
4. Define exactly which data crosses the boundary in Hybrid mode.
5. Preserve stable insertion indices and solve continuation rules.
6. Add output defaults, memory accounting, and at least one invariant test.
7. Update this README and the User Guide when the public workflow changes.
