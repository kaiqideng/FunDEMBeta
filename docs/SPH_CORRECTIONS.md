# SPH State, Coupling, and Search Corrections

Version 1.0.1 — September 11, 2026.

This note describes the implementation corrections identified by the September 2026 SPH review. It distinguishes committed numerical integration from temporary observations and does not claim that FunDEM implements every feature of SPHinXsys.

## 1. One current-time observation scope

Affected code: `solver::observeCurrentState`, `solver::writeOutput`, and the SPH output-snapshot hooks in `src/solver/SPHDEM.cpp`.

SPH acoustic updates can span several DEM steps. Between updates the stored fluid state represents an earlier time than the DEM clock. Previously, writing VTU temporarily advanced the fluid and then restored that deferred state; a frontend reading the particle container afterward obtained old positions with a new time label.

`observeCurrentState(callback)` now provides a single scope in which all observers see the same current-time state. The software captures scientific output, viewport data, and the restart checkpoint inside that scope. Nested output reuses the active scope. On exit, including callback exceptions, deferred state is restored so output frequency does not change the acoustic integration schedule. Reentrant `step`, `solve`, and `initialize` are rejected during an observation. The callback receives a const solver; mutating it through a separately captured alias is outside the observer contract.

An observation can temporarily rebuild search structures. Restoring the saved particles therefore invalidates the search references as well. The next interaction rebuilds neighbors if needed without reinitializing density, resetting prior force, or changing the advection phase.

## 2. A common held-load convention

Affected code: `src/interaction/virtualParticleCoupling`, virtual-particle load fields, and CUDA wall-load collection/scattering.

Between acoustic updates, all backends hold the last evaluated world force **and world torque**. Previously CPU/Hybrid held torque while GPU recomputed it from the current rotated lever arm. That difference was unrelated to floating-point reduction order.

The GPU now caches the world torque when the fluid load is evaluated and scatters that cached value until the next fluid update. New-minus-held load increments are retained for impulse matching. The added fields are included automatically in virtual-particle storage and in temporary SPH snapshots; host owner accumulators are included in coupling snapshots.

## 3. Complete, consistently timed wall kinematics

Affected code: SPHDEM `advanceCPU`, `advanceGPU`, `advanceHybrid`, and virtual-boundary kinematics sampling.

For finite-mass solids, kinematics are sampled at the left endpoint of each DEM step, before the first velocity kick and before force clearing. At that point the prior endpoint force contains contact, fluid, and external contributions. The external-force callback remains the last part of force assembly.

This fixes the previous use of partially assembled contact-only force to reconstruct wall acceleration. The averaged boundary acceleration used by wall-pressure reconstruction now has an explicit time convention. Prescribed infinite-mass boundary motion uses its sampled velocity history rather than a force divided by an infinite mass; constant velocity must not create fictitious acceleration.

## 4. Match transient reaction impulse at acoustic boundaries

Affected code: `src/execution/sphCouplingFunctions.h`, coupling load increments, and the final rigid-body velocity stage.

Let an acoustic interval contain N DEM steps of duration delta_t. With an old held fluid reaction F0 and a newly evaluated reaction F1, the split solid integration receives

```text
I_applied = delta_t * ((N - 1/2) * F0 + (1/2) * F1).
```

The corresponding fluid update uses the new reaction over the acoustic interval. After the final solid half kick, the implementation therefore adds

```text
Delta_I = (acoustic_dt - DEM_dt / 2) * (F1 - F0)
Delta_v = inverse_mass * Delta_I.
```

The analogous centroidal torque increment is mapped through the current world inverse inertia tensor to an angular-velocity correction. Gravity, external loads, and the gyroscopic term are not integrated a second time; endpoint position and orientation are not changed by this correction.

Only actual acoustic updates schedule this correction. Observation-only flushes do not modify solid velocity. They are current-time projections, not additional committed coupling steps. This removes the identified linear reaction-impulse imbalance; it is not a claim of exact energy conservation or exact global angular-momentum conservation for arbitrary rotating, moving bodies. Those remain subject to the chosen rigid-body integration and spatial discretization.

## 5. Retain computed-invalid-state information

Affected code: host density-rate assignment, CUDA density-rate kernel, and existing state-statistics reductions.

The public particle setter continues to reject invalid user input. Internal computed density rates are written to the same field used by the device path instead of silently falling back to a previous finite value. Once a computed rate is nonfinite it remains latched until initialization, so a later finite split stage or inlet constraint cannot conceal the earlier failure.

The existing end-of-acoustic-step statistics also inspect density rate, pressure, and position. No full `isValid()` pass was added inside neighbor interactions. Emergency density limits remain a last-resort numerical bound, not a substitute for reporting invalid calculations.

## 6. Reuse neighbors only while their displacement allowance is valid

Affected code: `src/interaction/SPHNeighborhood` and the `SPHDEM::ensureSPHNeighborhood` overloads.

Reference positions are stored separately from physical particle state at each search build. Before a density or pressure stage consumes neighbors, the maximum actual fluid/boundary displacement is compared with half the search skin. Since the relative motion of any pair is bounded by twice that maximum, a stage cannot silently use a list whose original search allowance has been exhausted.

Normally the advection-step rebuild is sufficient. If rapidly changing motion exhausts the skin sooner, only the search structures are rebuilt. The physical kernel still has support `2h`; there is no extra density reset or viscous-force refresh. CPU uses a parallel reduction and GPU reduces on the supplied stream. Reference-storage capacity is included in reported device memory.

The guard adds a linear reduction before the relevant stages. It is a correctness cost, not a claimed speedup; GPU latency and workload scaling should be profiled separately from the tests below.

## 7. Safe interval and source discretization conversion

Acoustic and advection interval ratios are capped to representable positive int values before conversion. Capping makes the chosen interval shorter, not longer than its stability limit. Oversized cylindrical jet lattices are rejected before unsafe casts or products; ordinary lattice dimensions and particle counts are unchanged. These checks occur at conversion/source setup boundaries, not in particle-pair calculations.

The DEM-versus-stability-limit comparison uses relative tolerance, not an absolute tolerance in seconds. This prevents extremely small acoustic limits from being treated as equal to much larger DEM steps. Circular source counts use exact lattice row counts with rounding correction, avoiding quadratic candidate enumeration before rejecting a source that cannot fit its indexed storage.

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
- `SPHCouplingTest`: held torque under rotation, changing-load reaction impulse with a known external-force budget, complete wall kinematics, and non-mutating observations;
- `SPHRobustnessTest`: invalid-rate latching, approaching particles crossing cells, boundary displacement, and extreme interval/jet conversion;
- software engine/JSON checks: paused and recorded positions, checkpoint time consistency, and frozen resumed SPH properties;
- `frameHistoryTest`: serialization, bounded cache, failure atomicity, corrupt/truncated records, concurrent access, and spool lifetime.

Local verification passed all six focused CPU/CUDA tests with CUDA 13.2 on an RTX 5060, all three CPU-only software engine/history/policy checks, and the standalone history tests under ASan/UBSan. Native Windows and macOS package checks are additionally required before publishing their binaries. Tiny regression tests are not substitutes for a full hydrostatic, impact-force, or long-duration FSI convergence study.

## Deliberately unchanged model choices

Density reinitialization still uses `max(rho_sum, rho0)`, matching the SPHinXsys FreeSurface scalar update at [commit 5dec24e](https://github.com/Xiangyu-Hu/SPHinXsys/blob/5dec24e4731e5ad8e45cc7f7191cda19871f410f/src/shared/particle_dynamics/fluid_dynamics/density_summation.hpp#L29-L32). Older paper variants and the NearFreeStream branch are different choices, not automatic replacements.

No unverified gradient-correction or particle-shifting model was enabled. The methodological references remain the [dual-criteria paper](https://arxiv.org/abs/1905.12302) and [multi-resolution FSI paper](https://arxiv.org/abs/1911.13255). Constrained jet particles still skip their own density evolution while contributing to neighbors' interactions.
