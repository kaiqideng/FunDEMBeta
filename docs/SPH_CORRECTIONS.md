# SPH State, Coupling, and Search Corrections

Updated October 24, 2026. The search and robustness notes originated in the September 2026 review; the coupling schedule below supersedes its endpoint impulse correction.

This note distinguishes committed numerical integration from current-time observations. FunDEM keeps SPH drift-kick-drift (DKD) and DEM kick-drift-kick (KDK); it does not claim to implement every feature of SPHinXsys.

## 1. One current-time observation scope

Affected code: `solver::observeCurrentState`, `solver::writeOutput`, and the SPH output-snapshot hooks in `src/solver/SPHDEM.cpp`.

An SPH acoustic interval spans `N` DEM steps. The fluid now advances first to the acoustic endpoint, then the solids catch up using `N` ordinary substeps. During those substeps the internal fluid state can therefore be ahead of the public DEM clock. Reading integration buffers directly is not a synchronized observation.

`observeCurrentState(callback)` provides one common-time host view for scientific output and viewport data. Nested output reuses that view. On exit, including callback exceptions, leading host buffers are restored without changing the numerical continuation. Reentrant `step`, `solve`, and `initialize` are rejected during an observation. The callback receives a const solver; mutating it through a separately captured alias is outside the observer contract.

For an observation inside an acoustic interval, the solver retains its acoustic-start state and uses:

- the DKD quadratic from stored endpoint positions and velocities: `x(s) = (1-s)*x0 + s*x1 - 0.5*H*s*(1-s)*(v1-v0)` for interval fraction `s`;
- linear velocity and density interpolation, with pressure recomputed from the equation of state;
- the unchanged acoustic-interval free-surface classification and interval-average density rate; only position, velocity, and density need acoustic-start storage because free-surface classification is updated before that storage is captured;
- jet flags evaluated at visible time, and virtual-boundary geometry rebuilt from the actual current solid pose;
- current rigid-point velocity for visible virtual particles, with acceleration averaged over the completed DEM samples in the current interval (instantaneous kinematics when no samples exist), rather than exposing raw device accumulator sums.

The position interpolant recovers constant-acceleration DKD motion and preserves both stored positions. For interacting fluid, these interior values are a dense-output numerical approximation, not a new SPH evaluation at the observation time. At common acoustic endpoints the committed solution itself is used. Observations do not evaluate another fluid force, modify contact/bond histories, rebuild searches, change advection phase, or write device state. Their frequency therefore does not alter the actual trajectory.

On `solve()` completion the synchronized host view remains available, while the leading endpoint is retained privately for continuation. Resuming an unchanged solver restores the leading buffers and finishes the remaining DEM substeps. Configuration changes or appended particles discard that continuation and restart from the current visible state; they do not resume an old acoustic interval with incompatible inputs. A saved frame inside an interval is a common-time state, not a serialized copy of the complete hidden multirate continuation.

Consequently, pausing and resuming the same solver object preserves its numerical continuation; constructing a new solver from visible frame values does not restore the private future endpoint, held reaction history, and search phase. The latter begins a new interval from the frame time and is not a bitwise continuation guarantee. Output itself never replaces the private state by the interpolated view.

## 2. A common held-load convention

Affected code: `src/interaction/virtualParticleCoupling`, virtual-particle load fields, and CUDA wall-load collection/scattering.

All backends hold the newly evaluated world force **and world torque** throughout the following solid interval. GPU caches the world torque when the fluid load is collected, rather than recomputing it from a different rotated lever arm at each solid substep. CPU, GPU, and Hybrid therefore use the same held-load convention.

Only the old hydraulic component is removed from the cached total force/torque when an interval begins. The new reaction is added before the first solid half kick. Replacing this component does not reevaluate contact/bond history, repeat gravity, or invoke external-field callbacks twice. Both pressure and the cached viscous prior force participate in the reaction delivered to the solids. No new-minus-held load increment or velocity-correction cache is retained.

The GPU acoustic-start snapshot is copied device-to-device. Ordinary GPU intervals do not add a complete fluid download/upload. Hybrid still transfers boundary kinematics and the wall reaction required by its CPU-owned solids; output reads the data needed for its host view on the solver stream.

## 3. Complete, consistently timed wall kinematics

Affected code: SPHDEM `advanceCPU`, `advanceGPU`, `advanceHybrid`, and virtual-boundary kinematics sampling.

At an acoustic boundary, prescribed-motion preparation runs before fluid evaluation. Virtual geometry is then refreshed from the current solid position and orientation. Boundary velocity and acceleration use averages accumulated over the preceding completed solid interval; startup uses instantaneous kinematics because no completed interval exists. This is an explicit, lagged boundary treatment, not a jointly implicit or exact midpoint fluid-solid solve.

For finite-mass solids, new kinematic samples are taken at each solid substep's left endpoint, after the new hydraulic interval load has been installed and before the first velocity kick or force clearing. The sampled total includes contact, fluid, and external contributions. The external-force callback remains last in ordinary force assembly. Prescribed infinite-mass motion uses differences of sampled world velocities rather than force divided by mass; constant prescribed velocity must not create fictitious acceleration.

## 4. Match reaction impulse through stage ordering, not correction

Affected code: SPHDEM `beginSPHInterval`, all three `advance...` methods, and hydraulic-load collection/scattering.

Let `H = N * delta_t` and let `F` be the fluid reaction collected from that leading SPH interval. Both solid velocity half kicks of every substep use that same hydraulic force:

```text
fluid: DKD over H, collecting pressure and viscous wall reaction F
solid: N * [K(delta_t/2, F), D(delta_t), K(delta_t/2, F)]
I_solid = N * delta_t * F = H * F.
```

The opposite linear momentum exchange is already represented in the fluid kick. The previous late-fluid arrangement instead mixed the old hydraulic force into the first solid half kicks and repaired that mismatch afterward. Its endpoint impulse time, load-difference members, and correction methods/kernels have been removed. Normal velocity integration still uses physical force; no additional velocity or angular-velocity correction is applied.

This removes the identified linear reaction-impulse mismatch at common acoustic endpoints. It does not establish exact energy conservation, exact global angular-momentum conservation for arbitrary rotating bodies, global second-order accuracy, or unconditional stability for large fluid-solid added-mass ratios. Those require separate spatial/model analysis and convergence studies. Observation inside an interval is subject to the dense-output approximation described above.

## 5. Retain computed-invalid-state information

Affected code: host density-rate assignment, CUDA density-rate kernel, and existing state-statistics reductions.

The public particle setter continues to reject invalid user input. Internal computed density rates are written to the same field used by the device path instead of silently falling back to a previous finite value. Once a computed rate is nonfinite it remains latched until initialization, so a later finite split stage or inlet constraint cannot conceal the earlier failure.

The existing end-of-acoustic-step statistics also inspect density rate, pressure, and position. No full `isValid()` pass was added inside neighbor interactions. Emergency density limits remain a last-resort numerical bound, not a substitute for reporting invalid calculations.

## 6. Reuse neighbors only while their displacement allowance is valid

Affected code: `src/interaction/SPHNeighborhood` and the `SPHDEM::ensureSPHNeighborhood` overloads.

Reference positions are stored separately from physical particle state at each search build. Before a density or pressure stage consumes neighbors, the maximum actual fluid/boundary displacement is compared with half the search skin. Since the relative motion of any pair is bounded by twice that maximum, a stage cannot silently use a list whose original search allowance has been exhausted.

Normally the advection-step rebuild is sufficient. If rapidly changing motion exhausts the skin sooner, only the search structures are rebuilt. The physical kernel still has support `2h`; there is no extra density reset or viscous-force refresh. CPU uses a parallel reduction and GPU reduces on the supplied stream. Reference-storage capacity is included in reported device memory.

The guard adds a linear reduction before the relevant stages. It is a correctness cost, not a claimed speedup; GPU latency and workload scaling should be profiled separately from the tests below.

## 7. Safe intervals and short-buffer source injection

Acoustic and advection interval ratios are capped to representable positive int values before conversion. Capping makes the chosen interval shorter, not longer than its stability limit. Oversized cylindrical jet lattices are rejected before unsafe casts or products; ordinary lattice dimensions and particle counts are unchanged. These checks occur at conversion/source setup boundaries, not in particle-pair calculations.

The DEM-versus-stability-limit comparison uses relative tolerance, not an absolute tolerance in seconds. This prevents extremely small acoustic limits from being treated as equal to much larger DEM steps. Circular source counts use exact lattice row counts with rounding correction, avoiding quadratic candidate enumeration before rejecting a source that cannot fit its indexed storage.

A finite-duration `SPHJet` now seeds only a short buffer, at most `ceil(2h/spacing)+1` layers, instead of preloading its entire nominal stroke. The complete finite layer count remains `floor(speed * duration / spacing)`. New layers are appended at their scheduled inlet times, quantized to a common DEM node by clipping the acoustic plan. Their positions compensate the at-most-one-DEM-step insertion delay. Insertion happens before the acoustic-start snapshot, never during output; start and end snapshots therefore contain the same active particle set within an interval.

Particle IDs remain stable and append-only. Device capacity growth preserves the old device prefix and initializes only the new tail; a stale host copy must never replace the evolving fluid. Source ownership is recorded explicitly because several jets and ordinary blocks can be interleaved. `SPHJetParticleIndices(firstIndex)` exposes the particles already created by one source. Active short-buffer particles receive prescribed velocity and reference density. Outlet release is permanent, source completion stops later replenishment, and a paused/resumed solver preserves its birth counters and emitted ownership.

Source end time, public handle, and created-layer count are derived from the source value and append-only ownership, rather than copied into separate mutable fields. A single released-layer prefix marks the remaining inlet suffix; per-layer global starts follow directly from ownership. Initial layers and later births use the same initializer and axial-position expression.

Nozzle clearance changes virtual-wall quadrature only. A finite cylinder spanning the upstream buffer is expanded by `spacing/2` in the radial and axial directions; virtual samples inside it receive zero volume and no SPH wall load. The mask is independent of emission completion and does not cut a downstream target wall. It does not modify an LS SDF, triangle mesh, or DEM contact boundary and is not a CAD-level hole. Moving samples regain their nominal lattice volume outside the mask. A changed mask clears its old sample loads and refreshes advection preparation, including viscous prior force. The owner-level old hydraulic load is removed before clearing those sample loads. Current-time observation changes only its temporary host mask, never device state or the actual neighbor/advection phase.

Finite inlets mask virtual-boundary quadrature rather than changing LS geometry: a sample inside any short upstream buffer, with half-spacing radial/axial clearance, has zero effective volume. Slots and owner indices remain stable, and samples outside the finite cut retain their usual volume. The mask persists after emission ends. Moving-boundary mask changes invalidate the neighborhood and restart advection preparation so previously cached viscous loads cannot survive a newly opened or restored sample. Observation updates only the displayed host mask; it never changes integration phases. The LS grid, surface geometry, and DEM contact detection are not cut by this SPH-only policy.

## 8. Workbench history and resumed fluid properties

These changes belong to FunDEMSoftware rather than the solver library:

- Recorded display frames and full restart checkpoints are stored in a process-owned temporary disk spool. Frame time/step metadata stays resident. The initial implementation accompanying this 1.0.1 review used a byte-bounded playback cache with a 32 MiB default. Oversized individual frames could be loaded without being permanently cached.
- Geometry resources are stored separately and reused. Both frame files must be written successfully before a frame is published. Records have format, count, and checksum checks; read or disk-full errors are reported rather than silently dropping history.
- Reset/project replacement releases the old spool after any in-progress reader finishes. Normal application shutdown removes it. This is temporary session storage, not a substitute for exporting a project or scientific files; abnormal process termination can leave temporary files behind.
- A frame checkpoint records its original uniform SPH properties. Those properties are read-only when resuming its saved SPH model, and project compilation rejects mismatches. Old JSON files without this metadata establish the baseline from their stored solver settings on import; earlier manual tampering cannot be reconstructed.

### Subsequent workbench update: FunDEMSoftware 0.3.27

The current workbench uses a configurable multi-frame LRU playback cache with a 256 MiB default. **Settings > Display Storage** exposes separate **Viewport** and **Playback Cache** budgets; the latter applies immediately, persists as an application preference across Reset/project replacement, and accepts **Off** to disable retention. Cached frames reuse immutable LS geometry, whose retained allocation is counted once across frames. Reducing the limit evicts cache entries, while oversized frames remain readable without retention. Complete render/checkpoint disk records and the temporary-spool cleanup policy are unchanged.

This follow-up changes workbench memory retention, not the numerical corrections recorded for core version 1.0.1 above. The default is an engineering tradeoff rather than a measured optimum: more retained frames can avoid repeated disk reads, but uncached reads and large frame copies may still delay playback. It is not a total-process RAM or GPU-memory limit and does not introduce lossy scientific history.

## Focused verification

The regression suites cover:

- `SPHCPUStateTest` / `SPHStateTest`: analytic current-time positions, nested and throwing observers, output/acoustic misalignment, and continuation across advection intervals;
- `SPHCouplingTest`: held torque under rotation, changing-load reaction impulse with a known external-force budget, first-half-kick installation, complete wall kinematics, and non-mutating observations;
- `SPHRobustnessTest`: invalid-rate latching, approaching particles crossing cells, boundary displacement, and extreme interval/jet conversion;
- `SPHJetTest`: bounded initial buffers, scheduled growth to exact finite counts, two-source/block ownership, preserved device prefixes, no early births during observation, permanent outlet release, non-grid births after a DEM-step change, finite virtual openings with unchanged LS grids/stable sample slots, static/moving-wall mask updates, and uninterrupted versus observed/segmented continuation;
- software engine/JSON checks: paused and recorded positions, checkpoint time consistency, and frozen resumed SPH properties;
- `frameHistoryTest`: serialization, bounded cache, failure atomicity, corrupt/truncated records, concurrent access, and spool lifetime.

The September review verified its earlier CPU/CUDA corrections and software history behavior. Those historical results are not evidence that a later coupling change has been tested: the leading-fluid schedule must pass the focused suites for the actual CPU, GPU, and Hybrid builds being delivered. Native Windows and macOS package checks are additionally required before publishing their binaries. Tiny regression tests are not substitutes for a full hydrostatic, impact-force, or long-duration FSI convergence study.

## Deliberately unchanged model choices

Density reinitialization still uses `max(rho_sum, rho0)`, matching the SPHinXsys FreeSurface scalar update at [commit 5dec24e](https://github.com/Xiangyu-Hu/SPHinXsys/blob/5dec24e4731e5ad8e45cc7f7191cda19871f410f/src/shared/particle_dynamics/fluid_dynamics/density_summation.hpp#L29-L32). Older paper variants and the NearFreeStream branch are different choices, not automatic replacements.

No unverified gradient-correction or particle-shifting model was enabled. The methodological references remain the [dual-criteria paper](https://arxiv.org/abs/1905.12302) and [multi-resolution FSI paper](https://arxiv.org/abs/1911.13255). Constrained jet particles still skip their own density evolution while contributing to neighbors' interactions.
