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
- `LSDEM`: unified materials, material-pair coefficients, level-set geometry, LS particles, LS–LS contacts, and LS–LS bonds.
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

1. Add materials, configure material pairs, then add geometries, particles, and bonds.
2. Configure boundaries, gravity, time step, and output.
3. Select the solver execution mode.
4. Call `solve(numberOfSteps)`.

The first `solve()` initializes automatically. At completion, device state and interactions are synchronized to the host and `initialized()` returns to false. An unchanged next solve resumes the retained continuation; appending objects or changing model configuration invalidates it and initializes the required storage from the synchronized host state. Configuration that would invalidate device state cannot change while initialized.

Solver-owned containers are exposed to users only as const references. Existing elements must never be swapped, reordered, replaced, or erased; add valid objects through the solver's `add...` interfaces so every insertion index remains stable.

`addMaterial(material{E, nu, restitution, density})` appends all new unordered
material pairs automatically, including the diagonal. Configure a standalone
`materialPair{i,j}` with stiffness/friction/damping setters and pass it to
`setMaterialPair()`. The pair table is read through `materialPairs()` and uploaded
alongside bulk materials on the solver stream. Appending materials preserves
existing pair indices and overrides. Pair changes are rejected while initialized.

Initialization checks configured nonnegative per-area stiffness for LS material combinations
that could interact, including same-material contacts. Two fixed-only LS bodies
do not impose irrelevant contact-law requirements. Sphere-sphere and sphere-LS
may instead use their Hertz-Mindlin defaults. These checks are not repeated as
full model-validation passes inside every simulation step.

Each pair has only two switches: `setLinearStiffnessEnabled()` selects all four
ordinary sphere-contact stiffnesses, and `setDampingEnabled()` selects all four
explicit dashpots. Both default to false. Disabling a group preserves its stored
values. LS per-area stiffness has no switch and must be assigned explicitly;
zero is configured, while the initial `-1` sentinel is unset.

Each DEM step uses a velocity-Verlet structure. Force assembly order is clear state → same-type contacts/bonds → cross-type contacts/bonds → user external fields. Protected virtual external-field interfaces receive a host container on CPU and the corresponding `device_type` on GPU.

### SPH Time Integration

`SPHDEM::setSPHProperties(spacing, smoothingLength, referenceDensity, dynamicViscosity)` configures uniform discretization and fluid properties shared by all SPH particles. Individual `SPHParticle` objects contain only dynamic state.

The solver groups its state by responsibility, rather than keeping inlet flags among integration parameters:

| Internal type | Responsibility |
| --- | --- |
| `SPHProperties` | Uniform physical inputs and derived constants; constructs current stage parameters without caching a second physical configuration. |
| `SPHTimeInterval` | Stability limit, quantized duration and integer interval; reused for acoustic and advection steps. |
| `SPHStepState` | Stability intervals, reduced fluid statistics, current acoustic interval length/completed DEM-step counts, and advection/search phase. |
| `SPHJetState` | Finite source schedules, short inlet-buffer indices, append-only emitted-particle ownership, and completion fast path. |
| `SPHIntervalStart_` | Only acoustic-start position, velocity, and density for dense observations; retained on the device in GPU/Hybrid mode. |
| `outputSnapshot` | Leading fluid/boundary host buffers temporarily replaced by a synchronized observation; no device or integration-state copy is restored. |

CPU and GPU collect statistics independently but share acoustic-limit checking, time quantization, search-buffer calculation and advection-phase advancement. CPU interaction inputs are constructed in one place. Cylindrical source generation shares the same regular cross-section lattice; a static downstream column and a time-dependent inlet remain distinct operations.

The current acoustic interval retains its own integer length because the stability calculation may already have selected a different length for the next interval. Its elapsed sample count, remaining substeps, and physical duration are derived from the current interval's two counters and the DEM step. A second fluid clock and initialized-population caches are unnecessary: new fluid intervals start only when the solid clock reaches the preceding acoustic endpoint, and changed models use normal initialization. CPU and Hybrid share host-solid advancement and force assembly; device SPH stages consume the common interaction-parameter record.

Source state likewise stores only essential topology and ownership. End time is derived from start time and duration; the public handle is the first owned index; the created-layer count is derived from membership and cross-section size. One released-prefix index replaces a separate per-layer metadata list. Initial seeding and later emission share the same layer initializer and world-position formula. Acoustic-start free-surface flags need no separate copy because classification is unchanged during a DKD interval.

`addSPHJet(inletCenter, velocity, radius, length)` creates one static downstream particle column without an ongoing inlet constraint. `addSPHJet(SPHJet{outletCenter, direction, radius, speed, duration})` creates an actual finite-duration inlet. With spacing `s`, smoothing length `h`, speed `U`, and duration `T`, its total axial layer count is `floor(U*T/s)`, but only `B = min(totalLayers, ceil(2h/s)+1)` upstream layers exist initially. Later layers are appended into the short inlet buffer at intervals `s/U`, until the finite total has been created. Already emitted particles are never recycled, reordered, or overwritten.

CPU, GPU, and Hybrid constrain only the active short-buffer particles: prescribed velocity, reference density, and equation-of-state pressure are imposed before optional advection preparation and again after velocity integration. These particles skip their own density integration while still interacting with neighbors. Crossing the outlet releases a particle permanently; returning into the pipe does not reapply its inlet constraint. Once the source duration ends no new layers or constraints are applied, and the completion fast path avoids later source work.

Births occur only at common acoustic starts, before saving the acoustic-start state. The acoustic schedule is shortened to the first DEM grid point at or after the next birth, so a birth is at most one DEM step late. Its insertion position compensates that delay to retain the regular lattice phase. No birth occurs during an observation. Each acoustic interval therefore has equal start/end particle counts, and dense output cannot expose a future layer early. GPU growth preserves existing device state and initializes only the appended tail on the solver stream; it does not upload stale host particles over the evolving fluid.

`addSPHJet()` returns the first seed particle index. Use `SPHJetParticleIndices(firstIndex)` to obtain the source's current append-only membership. Multiple sources and ordinary blocks can be interleaved, so later source particles are not necessarily one contiguous interval. `SPHJetsCompleted()` reports completion at the public visible time, independently of the leading fluid endpoint. Output and paused continuation preserve source schedules and membership instead of reinjecting particles.

An inlet opening is applied only to virtual SPH wall samples. Samples inside the finite upstream buffer cylinder, expanded by half a particle spacing axially and radially, receive zero quadrature volume. Their force, viscous prior force, and coupling torque are cleared when the mask changes. The nominal volume is restored when a moving sample leaves the opening. The opening remains after emission ends, never extends along the downstream jet path, and does not alter LS signed-distance grids, surface triangles, or DEM contacts. Moving-wall mask changes restart SPH advection preparation so no old viscous wall force survives the changed quadrature; the preceding held hydraulic load is removed before masked load fields are cleared. Observation applies this rule only to its reversible host view. With an SPH periodic domain, the virtual opening is also tested in the corresponding translated cell images; its source schedule and stored fluid coordinates remain unchanged.

The SPH acoustic time step is an integer multiple of the DEM time step without exceeding the acoustic limit; its minimum is one DEM step. If one DEM step already exceeds the acoustic or advection/viscous limit, initialization rejects the configuration. At the start of each advection step, `SPHDEM` builds the fluid and virtual-wall neighborhood on the selected SPH backend, updates the free-surface state, reinitializes density, and stores the viscous force as a prior force. Neighbor queries include a relative-displacement buffer derived from the maximum fluid and virtual-wall velocity and acceleration. Acoustic substeps reuse the search structures while `interaction/SPHNeighborhood` confirms that actual displacement remains within half that buffer. Before each density or pressure stage, exhausted or invalid references trigger a search rebuild without reinitializing density, refreshing the viscous prior force, or restarting the advection step.

SPH retains its drift-kick-drift (DKD) position/velocity update and two density half-steps. DEM retains kick-drift-kick (KDK). At a common acoustic boundary, the solver refreshes virtual-wall geometry from the current solid pose, uses completed-interval average wall velocity/acceleration (instantaneous kinematics at startup), and advances the fluid over `H = N * DEM_dt` first. The resulting pressure and viscous reaction force and world torque replace only the old hydraulic part of each solid's cached load, before the first solid half kick. DEM then takes `N` ordinary KDK substeps with this same hydraulic load; contacts, bonds, gravity, and external fields retain their normal DEM schedule. No endpoint velocity or angular-velocity impulse correction is applied. This is explicit partitioned coupling, not a claim of exact energy conservation, unconditional added-mass stability, or global second-order accuracy.

### Optional SPH Periodic Neighbors

`SPHDEM::setSPHPeriodicDomain(const math::periodicDomain&)` enables selected
axes for fluid/fluid and fluid/virtual-wall interactions on CPU, GPU and Hybrid.
The default domain disables every axis, preserving ordinary behavior. For example,
after configuring the usual SPH properties:

```cpp
fundem::math::periodicDomain periodic;
periodic.minimum_ = {0.0, 0.0, 0.0};
periodic.length_ = {1.0, 1.0, 1.0};
periodic.periodicAxes_[0] = true;
periodic.periodicAxes_[1] = true;
solver.setSPHPeriodicDomain(periodic);
```

Each enabled period L must have finite valid bounds and satisfy `4h < L`,
where h is the smoothing length: physical support `2h` is strictly less than
half a period. Both this setter and `setSPHProperties()` enforce that condition.
Changing the domain invalidates deferred continuation/search state and is
rejected while `initialized()` is true. Setting an identical value is a no-op,
including while initialized. Apply actual changes at the ordinary synchronized
lifecycle boundary, not concurrently with a running stage.

The API changes SPH neighbor images and translated virtual inlet openings only. Fluid integration, inlet constraints and
dense-output endpoint buffers keep continuous real coordinates, without ghost
particles or mass. An application that wants wrapped display/file coordinates
must wrap a detached observation/output copy; never feed those wrapped values
back into the retained integration or inlet buffers. This SPH setting does not
automatically make LS/DEM contacts or Bonds periodic; a boundary adapter must
provide the corresponding solid policy. All device grid/interaction helpers
use the caller's stream; the public setter configures metadata, not device work.

### Output

After `setOutputDirectory()` and `setOutputStepInterval()`, each output frame reports the calculated step count, frame index, physical time, and device memory usage. VTU output is binary by default; use `setVTUOutputFormat(vtuFormat::ascii)` for ASCII.

Core solvers enable VTU files by default. `setVTUOutputEnabled(false)` suppresses all particle and interaction VTU files without changing synchronization, frame indices, output cadence, selected optional fields, or `energy.dat`. Re-enable it to emit the normal default fields plus the retained optional selections. A graphical frontend may use `observeCurrentState()` for its own result format without enabling core VTU output; that frontend's file policy is not a core-library default.

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

## SPH Leading-Time Observations and Continuation

While solids catch up, the committed internal SPH state is at the end of the acoustic interval and may be ahead of `time()`. Use `observeCurrentState(callback)` to read a common-time host view rather than reading integration buffers directly. For an interior fraction `s` of an acoustic interval of length `H`, the DKD quadratic is `x(s) = (1-s)*x0 + s*x1 - 0.5*H*s*(1-s)*(v1-v0)`. It uses the stored endpoint positions and velocities and recovers the constant-acceleration trajectory. Velocity and density use linear interpolation, and pressure is recomputed from density. Interacting-fluid values inside the interval are a dense-output numerical approximation, not a separately reevaluated SPH solution. At common acoustic endpoints the stored SPH solution itself is used.

Observation updates visible virtual-wall geometry from the actual current solids, evaluates inlet flags at the visible time, and never reruns SPH forces, advances contact/bond history, or alters device state. Nested or throwing observers retain the same continuation. Output cadence therefore does not split the integration interval. `solve()` completion leaves the synchronized host view available and retains the leading endpoint privately; a subsequent unchanged solve restores that endpoint and finishes the remaining solid substeps. Particle-count or configuration changes discard the private continuation and reinitialize from the current visible host state. The same contracts apply to CPU, GPU, and Hybrid execution.

Pausing and resuming the same solver object is different from creating a new solver from a saved frame. In-memory continuation retains the private leading endpoint, interval counters, held loads, and search phase. A frame containing only visible particle/contact/bond values does not contain all of that private state; loading it starts a new synchronized interval and is not promised to reproduce an uninterrupted run bit for bit. Observation buffers reuse their allocations on restoration rather than repeatedly allocating full particle copies.

## Extension Checklist

1. Decide which solver layer owns the new particle type or coupling.
2. Add one full stage for every supported execution mode rather than scattering backend checks.
3. Keep device helpers stream-explicit.
4. Define exactly which data crosses the boundary in Hybrid mode.
5. Preserve stable insertion indices and solve continuation rules.
6. Add output defaults, memory accounting, and at least one invariant test.
7. Update this README and the User Guide when the public workflow changes.
