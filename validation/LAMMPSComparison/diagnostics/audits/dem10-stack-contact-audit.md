# dem10-stack: source and contact mapping audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

The reference is the LAMMPS granular verification input dem10-stack-3d-si.yaml (retained locally; not included), pinned at [29dd3fd7642a6fe9e8b4dd99120df4068b221250](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-stack-3d-si.yaml). Its exact Git blob is `405583bd3e27294ab46506f6ba242485a756de33`. This source carries the **granular unvetted** tag. It is an official repository test input, with that qualification retained.

The source manifest (retained locally; not included) records thirteen live fetched official files, their Git blobs and SHA256 hashes, and fifteen read-only snapshots of current native contact and integration sources. snapshot_dem10_stack_sources.py (retained locally; not included) validates the exact input bytes and snapshots. The already used newer LAMMPS pin is reused. The historical 24-input catalog remains a separate fixed collection.

## Parameters and coordinate transformation

| Quantity | Official value | Native mapping |
|---|---:|---|
| Sphere radius | 0.0005 m, both | Same |
| Lower / upper density | 20000 / 10000 kg/m³ | Same |
| Lower / upper initial center | y = 0.00045 / 0.00135 m | z = same values |
| Floor / ceiling | y = 0 / 0.0018 m | z = 0 / 0.0018 m |
| Gravity | (0, −9.81, 0) m/s² | (0, 0, −9.81) m/s² |
| Normal stiffness | 1000 N/m | Linear kn = 1000 N/m for all three active pairs |
| Restitution | 1 | Automatic normal damping = 0 |
| Tangential law | linear_nohistory, stiffness 0, friction 0 | kt = 0, tangential damping = 0, friction = 0 |
| Original step / duration | 2 μs / 10000 steps = 20 ms | Same |
| Half step / duration | — | 1 μs / 20000 steps = 20 ms |

The proper rotation is `(x,y,z) -> (x,-z,y)`, with determinant +1. Its inverse maps native vectors to `(x,z,-y)`. Thus native center z, velocity vz and force Fz correspond to official y, vy and Fy. Spheres are invariant under this rotation. There is no reflection or changed physical orientation.

The two finite spheres have masses 1.047197551196598e−5 and 5.235987755982989e−6 kg. Their moments of inertia are 1.047197551196598e−12 and 5.235987755982990e−13 kg·m². The floor and ceiling are immovable plane surfaces. Their finite mesh extent is 0.01 m and grid spacing is 0.0005 m; the centered vertical particle motion stays on the flat central region.

## Official acceptance and saved-state diagnostics

The analytic checker (retained locally; not included) selects `stack_energy` at zero-based segment 4, after five 2000-step segments. The shared driver (retained locally; not included) calls the analytic check after each segment and the checker evaluates it only at the selected final segment.

With native z equal to official y, define positive overlaps

- δf = max(0, R − z1);
- δp = max(0, 2R − (z2 − z1));
- δc = max(0, z2 + R − H).

The mechanical energy is

`E = 0.5 m1 v1² + 0.5 m2 v2² + g(m1 z1 + m2 z2) + 0.5 kn(δf² + δp² + δc²)`.

The reference uses the live masses/radii, the prescribed initial centers and zero initial velocities. The original initial overlaps are **50 μm / 100 μm / 50 μm**. Initial spring energy is 7.5e−6 J, gravitational energy is 1.155713397439345e−7 J and total **E0 = 7.615571339743934e−6 J**.

The actual official analytic condition is **|E(20 ms) − E0| / |E0| ≤ 0.0005**, equivalent to 0.05%. This final-time test is preserved for all four native runs. Whole-trace energy variation, force reconstruction and exact normal-mode trajectories are additional diagnostics.

The YAML's `epsilon: 1e-10` belongs to same-code saved position/velocity/torque/spin regression checks. It is not used as the cross-code acceptance tolerance. Native states at all five saved segment times are also exported for comparison after applying the inverse coordinate rotation.

## Force law and initialization

The pinned Hooke implementation (retained locally; not included) returns `k * delta`. Restitution 1 gives zero damping in the restitution damping implementation (retained locally; not included). Native contact material pairs enable the linear normal law and automatic restitution damping, with zero manual damping. No tangential, rolling or twisting coefficients are active.

The floor, pair and ceiling normal forces are `kn δf`, `kn δp` and `kn δc`. Their initial magnitudes are 0.05 / 0.10 / 0.05 N. Initial net contact forces on lower and upper spheres are −0.05 and +0.05 N along native z. Gravity is added separately by the velocity integrator and is excluded from the contact force accumulator.

The fixed walls exchange momentum with the spheres and do no work. Sphere momentum alone is therefore not conserved; floor and ceiling impulses plus gravity complete its balance. Mechanical energy includes both spheres' gravitational potential and each of the three spring energies once.

**The native initial CSV callback precedes solver initialization.** It retains zero force/contact channels at t = 0 although the supplied geometry is already compressed. The actual SphereDEM initialization snapshot (retained locally; not included) calls `assembleForceAndTorque(0.0)` before its first half kick on both CPU and GPU. Raw CSV data remain unchanged. Energy at t = 0 is reconstructed directly from the supplied geometry, and the effective initial contact force used for integration checks is computed from that geometry. The independent audit verifies the first displacement and full Verlet recurrence against those physical initial forces; it must not silently treat the raw preinitialization zero as the first-step force.

There are no fitted forces, initial-condition adjustments, mass changes, manual force injection, pre-runs or source modifications. The raw preinitialization sample and the physical bootstrap are disclosed separately.

## Independent checks and retained evidence

audit_dem10_stack_physics.py (retained locally; not included) independently reconstructs every recorded contact from geometry, checks contact action/reaction and each sphere's force, computes total energy, and verifies gravity-inclusive impulses, drift and half-kick updates. It also checks zero transverse motion, spin and torque. CPU/GPU comparisons and the two time steps are evaluated separately.

While all three contacts remain compressed, the system has the exact constant matrix `K = kn [[2,-1],[-1,2]]` and load `b = kn [-R,R+H] − M g`. A mass-normalized symmetric eigenproblem yields a continuous normal-mode reference without fitted parameters. The same eigenmodes with `cos(theta) = 1 − omega² dt²/2` give an exact discrete velocity-Verlet reference, including its correct physical bootstrap. This supplemental reference is valid only if every recorded overlap stays positive, which the audit checks.

Each backend and time step keeps its original project, native logs, every-step trajectory/contact CSVs, native `run.json` and frames. A complete native ZIP is retained on Windows. The compact WSL mirror contains configurations, logs, CSVs, source snapshots, reports, plots and audits; large native frame trees and their ZIP stay at the Windows delivery root. Old reports, source snapshots, results and their published statuses remain unchanged.

