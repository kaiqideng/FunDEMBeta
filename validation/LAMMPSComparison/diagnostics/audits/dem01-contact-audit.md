# DEM01 contact-model mapping audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. Scope: independent, read-only source review. No solver source, material parameter, or project parameter was changed by this review. This note does not claim a simulation pass; measured results belong in the case report.

## Reference and initial conditions

The reference is [LAMMPS DEM01, Hertz two-sphere collision](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem01-hertz-twosphere-3d-si.yaml), pinned to commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`. The local YAML has Git blob SHA `9029229be2a7123a3a5c9aef4d04e2df5a0bc9e0`.

| Input | Official value |
|---|---:|
| Sphere diameter | 0.005 m |
| Density | 2700 kg/m3 |
| Young's modulus | 7.0e10 Pa |
| Poisson's ratio | 0.33 |
| Restitution input | 0.5 |
| Sliding friction | 0.4 |
| Initial centers | (-0.002550, 0, 0), (+0.002550, 0, 0) m |
| Initial velocities | (+2, 0, 0), (-2, 0, 0) m/s |
| Initial angular velocities | (0, 0, +22), (0, 0, -22) rad/s |
| Time step | 1.0e-8 s |
| Run segments | 2600, 1500, 1000 steps |

Derived, rather than fitted, properties are mass per sphere `1.7671458676442588e-4 kg`, reduced mass `8.835729338221294e-5 kg`, reduced radius `0.00125 m`, effective Young's modulus `3.927729772191673e10 Pa`, and effective shear modulus `7.87897888433659e9 Pa`.

## Normal elasticity and automatic damping

Let `a = sqrt(R_eff * overlap)` and `K_sec = (4/3) E_eff a`. Both implementations use Hertz elastic force `F_el = K_sec * overlap` and the same reduced mass and radius. In FunDEM, the sphere pair must have `linearStiffnessEnabled = false` and `dampingEnabled = false`. The latter selects restitution-derived contact damping; it does not disable the physical dashpot. Any separately configured global velocity damping must be zero.

For the identical materials, FunDEM's harmonic-mean restitution mixing returns the same `e = 0.5`. Its dimensionless factor is

`beta = -ln(e) / sqrt(pi^2 + ln(e)^2) = 0.2154537619662468`.

The automatic normal dashpot is

`c_n = 2 sqrt(5/6) beta sqrt(m_eff * (3/2) K_sec)`.

This is algebraically identical to the pinned LAMMPS `hertz/material ... damping coeff_restitution` implementation. In particular, the `3/2` factor multiplies the secant stiffness to obtain the Hertz tangent stiffness. The dashpot varies with overlap; `0.21545` is not a constant dashpot in N s/m.

FunDEM source: `src/execution/contactMaterialFunctions.h`, lines 109-140; `contactFunctions.h`, lines 89-103; and `contactDetection.h`, lines 19-29, relative to the historical `FunDEMBeta/` source tree.

Pinned primary sources: [LAMMPS damping implementation, lines 147-183](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_damping.cpp#L147-L183), Git blob `fa5addccc5d3ba9f6aa1462c7753b89c4cd8e9b8`; [Hertz material implementation, lines 189-237](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_normal.cpp#L189-L237), Git blob `5ecbb5cb562e93b9dba827f6d6a4f18c6b86f495`.

## Unloading-force convention

DEM01 does not enable `limit_damping`. The pinned LAMMPS model defaults that switch to zero and therefore retains the signed sum of the elastic and viscous normal forces while overlap is positive. FunDEM's current `calculateLinearContactForce` also retains that signed force; its normal-force truncation line is commented out. These conventions agree for this reference.

FunDEM reports an absolute normal-force magnitude separately from the signed force vector. Use the force vector projected onto the contact normal when checking an attractive unloading tail. An unsigned magnitude cannot show the force sign. This tail follows the chosen dashpot convention and is not an added cohesive material law.

Sources: FunDEM `contactFunctions.h`, lines 244-267; [LAMMPS force assembly, lines 480-495](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L480-L495), Git blob `bab3889de43c3ac559977058b13b3a5fa1d92524`.

## Why sliding and tangential damping are inactive

The contact point lies midway between the identical spheres. On the left sphere, `omega = +22 ez` and its contact lever arm points in `+ex`. On the right sphere, `omega = -22 ez` and its lever arm points in `-ex`. Their rotational contact-point velocities are equal, both in `+ey`. Their translational velocities have only normal components, so the relative tangential velocity is zero. Equal and opposite spins therefore remain unchanged when rolling and torsional friction are zero, as in the reference.

The LAMMPS input specifies Mindlin tangential stiffness, `c_t = 0.5 c_n`, and sliding friction 0.4. FunDEM uses the same default Mindlin stiffness but its automatic `c_t / c_n = sqrt(K_t / S_n)`, approximately `0.8957651531` here. This is a real difference in the general tangential model. It has no effect on this symmetric normal collision because tangential velocity and tangential history remain zero. Preserve the official opposite spins and friction 0.4, and explicitly set rolling and torsional friction to zero. Passing DEM01 validates normal restitution and automatic normal damping, not the tangential damping law or sliding friction.

Sources: FunDEM `contactFunctions.h`, lines 68-86 and 162-194; [LAMMPS contact kinematics, lines 449-478](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L449-L478); [official granular documentation](https://docs.lammps.org/pair_granular.html).

## Time integration and interpretation of the reference

FunDEM CPU and GPU both perform half velocity kick, full position drift, force evaluation, and final half velocity kick. Normal viscous forces use the velocity available at force evaluation, normally the intermediate half-step velocity. The pinned [LAMMPS sphere integrator](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/fix_nve_sphere.cpp#L98-L137) uses the same ordering; its Git blob is `21520d7a695b3256b090194d66d99b20a5a2af31`.

One protocol detail can matter for very strict numerical comparisons: the official [test driver, lines 118-133](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_dem_common.cpp#L118-L133) executes each segment with `run N post no`, leaving the default `pre yes`. LAMMPS therefore performs run setup and force evaluation at each segment boundary. An uninterrupted FunDEM run retains the previously assembled force across those boundaries. With velocity-dependent damping, this can introduce a small finite-time-step discrepancy. This is a protocol consideration, not a demonstrated cause of any measured discrepancy; the case report must identify whether its runner mirrors those segment setups.

The YAML gives checkpoints at cumulative times 26, 41, and 51 microseconds. At 26 microseconds the spheres are still in contact and have not completed restitution. At 41 and 51 microseconds, their recorded x velocities are `-0.99988365289612324` and `+0.99988365289612324 m/s`, giving the official stored restitution `0.4999418264480616`. The input target remains `0.5`. These two values must be distinguished in the report.

The official analytic check uses `e_measured = -(v_left_x - v_right_x)/(2 vin)` with relative tolerance `0.02` and a separate total x-momentum bound. Its YAML `epsilon = 1e-10` applies to regression against the stored trajectory, rather than to the analytic restitution target. Both checks can be reported, with separate pass/fail status. Time-step refinement should be compared with the analytic target as well as the finite-step stored values. Rotational kinetic energy is initially nonzero and remains constant, so restitution-energy checks should use translational relative kinetic energy rather than total kinetic energy.

Driver Git blob: `db954a77c073fad55f9078df9257810dd507f26c`. Analytic-check source: [test_analytic_models.cpp](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp), Git blob `42bcbbb906377cf32f146ffbb608073b8f0ec793`.

## Recovering the dashpot from native saved frames

`SimulationSession::step` calls the native step before it captures an interval frame. `captureInteractions` copies native contact geometry and forces without recalculating them. Consequently, in a completed-step DEM01 frame at time `t_n`, position, overlap, and the contact force belong to the force assembly after the position drift; particle velocities have already received the final half kick. Recomputing the dashpot force using those completed-step velocities directly introduces a time-staggering error. CPU and GPU observations merely synchronize and copy this state; the DEM solver does not perform a new contact-force evaluation for a snapshot.

Use each contact's explicit master/slave references. Its unit normal `n` points from slave to master, and its force vector `F_c` acts on the master; the opposite force acts on the slave. These statements follow directly from `detectSphereContact`, `contact.h` field contracts, and `contactForceKernel.cu`, lines 131-138. Define `f_n = dot(F_c, n)`, retaining its sign, and let `u_n = dot(v_master_n - v_slave_n, n)` use the saved completed-step velocities.

For this isolated, force-free two-sphere environment, undoing the final half kick gives the exact force-evaluation normal speed, up to floating-point error:

`u_eval_n = u_n - (dt / (2 m_eff)) f_n`.

More generally, undo the half kick on each particle using its saved total force and gravity, then subtract the reconstructed velocities and project onto the normal. Equal gravity cancels in the relative velocity. This generalized step is needed if additional forces are present; substituting only the contact force would then be insufficient.

The independently inferred damping coefficient is

`c_n_inferred = (F_el_n - f_n) / u_eval_n`,

where `F_el_n = (4/3) E_eff sqrt(R_eff * overlap_n) overlap_n`. Compare this inferred value against `2 sqrt(5/6) beta sqrt(m_eff * 2 E_eff sqrt(R_eff * overlap_n))`. This reconstruction uses the recorded force and completed state, rather than calling FunDEM's damping function. Keep only positive overlaps and speeds sufficiently far from zero; the coefficient is unidentifiable at zero normal speed. Report the speed/overlap thresholds and the number of retained contact samples.

Two independent timing checks are available when every step is exported:

1. Reconstruct each particle's force-evaluation velocity as `(x_n - x_(n-1)) / dt`. In DEM01 the positions are unwrapped and there are no prescribed-velocity modules, so this is exactly the drift velocity. Differencing nearby positions amplifies text-rounding noise, so the exported data must retain double precision.
2. Verify `v_n - v_(n-1) = (dt/2) (a_(n-1) + a_n)`. Equivalently, reconstruct the evaluation velocity as `v_(n-1) + (dt/2) a_(n-1)` and compare it with `v_n - (dt/2) a_n`. At a real reinitialization boundary, use the force that was actually assembled for that boundary rather than a stale pre-boundary force.

A well-conditioned companion check is the signed force residual `f_n - (F_el_n - c_n_model u_eval_n)`. Normalize it by a documented force scale, for example the maximum of the elastic-force magnitude, damping-force magnitude, and a small absolute floor. This remains useful around maximum compression, where division by a nearly zero speed makes the inferred-coefficient error unstable. Do not use `normalForceMagnitude` as `f_n`, and do not average two completed-step velocities to guess the half-step velocity.

Snapshot sources: `archived application-adapter provenance`, lines 283-310, 351-433, and 708-742; `archived application-adapter provenance`, lines 37-60. Native ordering: `FunDEMBeta/src/solver/SphereDEM.cpp`, lines 194-217; `src/solver/solver.cpp`, lines 244-271 and 383-394. Native geometry: `src/execution/contactDetection.h`, lines 193-224. Component-prefixed paths are relative to the historical Platform source root; the other paths continue within the component named immediately before them.
