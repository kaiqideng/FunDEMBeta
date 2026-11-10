# DEM04 Hooke-history oblique-impact contact audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. Scope: independent source review and parameter mapping. No solver source or simulation input was modified by this audit. Measured validation results are reported separately in the case report.

## Pinned official test and acceptance criterion

Reference: [LAMMPS dem04-hookehistory-oblique-3d-si.yaml](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem04-hookehistory-oblique-3d-si.yaml), commit f9720ce2d176d615495d1d6c31349b1f71be1c73, Git blob 6ab54dcfe7fcef2d4dcc16eb430a29ba8ebe8f3a. The test records LAMMPS 4 Jul 2026 and a generation date of 21 Jul 2026.

One solid sphere approaches a stationary horizontal plane at z = 0. The setup uses SI units and three spatial dimensions.

| Quantity | Official input |
|---|---:|
| Diameter | 0.005 m |
| Density | 2700 kg/m3 |
| Initial center | (0, 0, 0.002550) m |
| Initial velocity | (2.5069, 0, -2.9876) m/s |
| Initial angular velocity | zero |
| Normal Hooke stiffness kn | 2.0e7 N/m |
| Tangential history stiffness kt | 5.7e6 N/m |
| Restitution input e | 0.98 |
| Sliding friction coefficient mu | 0.092 |
| Tangential damping multiplier | 0.0 |
| Gravity, rolling resistance, torsion, cohesion | none |
| Original time step | 1.0e-8 s |
| Run segments | 1724, 1133, 900 steps |
| Total duration | 37.57 microseconds |

The official analytic check is oblique_impact at analytic_segment = 1, the second cumulative checkpoint: 2857 steps, 28.57 microseconds. The first checkpoint is at 17.24 microseconds while the sphere is still in contact; the final checkpoint is at 37.57 microseconds. A half-time-step run should keep these physical times and double the corresponding step counts.

The [official analytic source](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp#L94-L111) checks all three of the following values with individual relative tolerance 0.01:

- vz_out = e * vz_in = 2.927848 m/s.
- vx_out = vx_in - mu * (1 + e) * vz_in = 1.962678784 m/s.
- omega_y_out = (5 / 2) * mu * (1 + e) * vz_in / R = 544.221216 rad/s.

All three must pass. The YAML epsilon = 1e-10 is a separate same-code stored-trajectory regression tolerance; it is not the cross-code acceptance criterion.

The reference stored rebound is vx = 1.9626463236052916 m/s, vz = 2.9278279037320356 m/s, omega_y = 544.25367639471085 rad/s. Report differences from those finite-step values separately from the analytic test.

## Exact normal and tangential mapping

Derived sphere mass m = (4 pi / 3) rho R^3 = 0.00017671458676442588 kg. Solid-sphere inertia I = (2 / 5) m R^2 = 4.4178646691106474e-10 kg m2. LAMMPS wall effective mass is the entire particle mass m, because the wall has no dynamics. FunDEM must retain the same sphere mass and inertia with a fixed wall.

LAMMPS Hooke elasticity is Fn_el = kn * delta. For damping coeff_restitution, its nondimensional prefactor for Hooke is

beta = -ln(e) / sqrt(pi^2 + ln(e)^2),
cn = 2 * beta * sqrt(m * kn).

For this case beta = 0.006430588502366038 and cn = 0.7645958085329603 N s/m. There is no sqrt(5/6) Hertz factor in the Hooke formula.

The official tangential damping is an xt multiplier on the normal damping prefactor. Here xt = 0, so ct = 0 N s/m exactly. FunDEM's automatic option derives a nonzero tangential dashpot whenever kt and beta are nonzero; therefore the correct mapping for this particular official test uses:

- linearStiffnessEnabled = true, normalStiffness = 2.0e7, slidingStiffness = 5.7e6.
- dampingEnabled = true, normalDampingCoefficient = 0.7645958085329603, slidingDampingCoefficient = 0.
- Sliding friction = 0.092; rolling/torsional friction and dashpots = 0.
- Global velocity damping = 0; gravity = 0.

Here dampingEnabled = true selects the explicit physical dashpot coefficients. It does not mean that the normal coefficient has been fitted: cn is calculated directly from the official e, m, and kn. This case validates the mapped explicit Hooke damping and tangential friction; it does not independently exercise FunDEM's automatic-damping branch. The previous DEM01 case addresses automatic damping.

[Official damping code, lines 147-183](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_damping.cpp#L147-L183) and [Hooke code, lines 145-169](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_normal.cpp#L145-L169) establish this mapping.

## Normal-force sign and tangential history

LAMMPS defaults limit_damping = 0 and the official input does not enable it. Therefore the signed force Fn = kn delta - cn vn is retained for positive geometric overlap, including a small attractive tail during separation. Its Coulomb limit uses mu * abs(Fn).

FunDEM's calculateLinearContactForce has normal-force truncation commented out. It returns the same signed normal force while separately reporting normalForceMagnitude = abs(Fn). The reported magnitude must not be substituted for signed normal force in impulse integration.

LAMMPS linear_history increments the spring history using the tangential contact velocity, calculates the elastic plus viscous force, caps the total force to the Coulomb disk, and writes back the history corresponding to the capped elastic part. FunDEM projects and increments the spring history, caps its elastic force, and then limits the added dashpot force to the same disk. These are different general algorithms when ct is nonzero. In this benchmark ct = 0 and the plane normal is fixed, so both reduce to the same radial return:

xi_trial = xi_old + dt * u_t,
Ft = clip(-kt * xi_trial, -mu abs(Fn), +mu abs(Fn)),
xi_new = -Ft / kt when capped.

LAMMPS's rotation-rescaling of history and FunDEM's simple tangential projection also coincide here because the wall normal remains exactly +z. This benchmark does not test rotating contact-normal transport or nonzero tangential dashpots.

Sources: [linear_history implementation, lines 121-188](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_tangential.cpp#L121-L188), [signed force and friction scale](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L484-L510), FunDEMBeta src/execution/contactFunctions.h lines 105-267 and contactMaterialFunctions.h lines 94-140.

## Contact lever-arm difference

The implementations retain a real geometric convention difference:

- LAMMPS wall contact: u_t = vx - R omega_y and tau_y = -R Fx.
- FunDEM sphere/level-set plane: contact point is x_center - (R - delta) n, so u_t = vx - (R - delta) omega_y and tau_y = -(R - delta) Fx.

FunDEM retains nominal R for mass, inertia, and constitutive properties. Only the contact kinematics and the friction torque arm use R - delta. Consequently good agreement within the official 1% tolerance must not be presented as exact algorithmic equivalence.

The source does not document why that convention was originally selected. Its reuse of a shared contact point is an implementation fact; it is not evidence that the convention is physically more accurate than the LAMMPS choice.

Sources: [LAMMPS contact kinematics, lines 464-477](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L464-L477) and [wall torque, lines 498-510](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L498-L510); FunDEMBeta src/execution/contactDetection.h lines 181-185, contactFunctions.h lines 70-85, and execution/cpu/particleFunctions.h lines 112-116. CPU and CUDA both add (point - center) cross force.

## Additional impulse and spin checks

These are supplementary diagnostics, distinct from the official three-output tolerance.

With zero gravity and no initial spin, ideal gross sliding predicts:

- Signed normal impulse Jn = m * (1 + e) * vz_in = 0.0010453459488464496 N s.
- Signed tangential impulse Jt = -mu * Jn = -9.617182729387336e-5 N s.
- Nominal-radius angular impulse Ly = -R * Jt = 2.404295682346834e-7 N m s.

Native finite-step identities should hold to floating-point precision:

- Jn_integrated = m * (vz_final - vz_initial).
- Jt_integrated = m * (vx_final - vx_initial).
- Integral tau_y dt = I * (omega_y_final - omega_y_initial).
- For each positive-overlap sample, abs(Fx) <= mu * abs(Fz).
- During capped sliding, abs(Fx) / abs(Fz) is 0.092.
- tau_y = -(R - delta) Fx for the native sphere, with y velocity, x/z spin and out-of-plane force/torque remaining zero.

Use a signed trapezoidal force/torque sum when every numerical step is exported. Initial loading can contain an elastic-history buildup before reaching the Coulomb cap, and the untruncated normal force has an attractive separation tail; hence abs(Jt)/Jn need not equal mu to machine precision. The official analytic formula idealizes gross sliding.

The measured torque-weighted arm, integral tau_y dt / (-Jt), can be compared with R and with the impulse-weighted R - delta. This quantifies the lever-arm contribution without modifying the solver.

## Force-evaluation timing of native saved samples

Native SphereDEM uses velocity-Verlet in both CPU and CUDA paths: first half velocity/angular-velocity kick, position/orientation drift, contact detection and history/force assembly, then the final half kick. SimulationSession archives the completed step; captureInteractions copies existing contact fields and does not recalculate forces.

For this zero-gravity, one-contact case, recover the velocities used in the saved force assembly with:

v_eval = v_saved - dt * F_saved / (2m),
omega_eval = omega_saved - dt * torque_saved / (2I).

For the native plane, reconstruct:

vn_eval = vz_eval,
ut_eval = vx_eval - (R - delta) * omega_y_eval,
Fn_model = kn delta - cn vn_eval.

The completed-step velocity-increment identities are:

v_n - v_(n-1) = dt * (F_n + F_(n-1)) / (2m),
omega_n - omega_(n-1) = dt * (tau_n + tau_(n-1)) / (2I).

These provide a timing/serialization cross-check. Do not evaluate damping using completed-step velocities directly and do not use contact.torque as the entire particle friction torque: contact.torque contains the extra rolling/torsional resistance moment, which is zero in this case. Particle torque includes the moment of the contact force.

Sources: FunDEMBeta src/solver/SphereDEM.cpp lines 194-217, src/interaction/contact.cpp lines 49-63, and src/execution/cuda/contactForceKernel.cu lines 130-138; archived application-adapter provenance lines 283-310 and 708-742.

## Source retention and limitations

The files in reference/lammps/dem04-source are complete snapshots fetched from the pinned official GitHub commit. source-manifest.json records their official and independently computed Git blob hashes and SHA-256 hashes; all eight Git hashes match.

Nine relevant FunDEM source files have been copied into reference/fundem/dem04-source without editing the originals. Their byte counts and SHA-256 hashes are also recorded in the manifest. These snapshots document the source inspected for this audit. The run manifest must independently document the actual executable and native libraries: source inspection alone does not prove byte identity with an installed GUI executable.

LAMMPS uses periodic x/y boundaries in this input. The sphere moves less than 0.1 mm in x, much less than its 0.4 m periodic domain length, with no lateral contacts or forces. Therefore the periodic image convention is inactive over this trajectory; it should still be disclosed if the FunDEM project uses a nonperiodic exterior.

The scope is this fixed-plane, zero-spin initial condition, linear Hooke normal model, zero tangential dashpot, and gross-sliding oblique rebound. It does not establish the behavior of general Mindlin variants, nonplanar contacts, or static-to-dynamic friction transitions.


The official test driver calls a separate run command for every segment (test_dem_common.cpp, lines 118-119). The fixed-wall setup method calls post_force again while disabling history advancement (fix_wall_gran.cpp, setup and post_force). In a velocity-dependent force model, that re-evaluation can make a tiny difference from a single continuous event. This is a source-based explanation of the very small translational reference discrepancy, not an independently rerun LAMMPS result.


## Independent completed-run cross-check

All four completed CSV runs were read independently by tools/audit_dem04_physics.py. The machine-readable output is tools/dem04-physics-audit.json. No main analysis script or report was changed.

| Run | vx at 28.57 us (m/s) | vz (m/s) | omega_y (rad/s) | Spin relative error | Official three checks |
|---|---:|---:|---:|---:|---|
| CPU dt | 1.962646252085700 | 2.927827884655243 | 542.750639025781 | 0.270216767% | PASS |
| GPU dt | 1.962646252085699 | 2.927827884655242 | 542.750639025781 | 0.270216767% | PASS |
| CPU dt/2 | 1.962651011469064 | 2.927892057493618 | 542.745838072904 | 0.271098936% | PASS |
| GPU dt/2 | 1.962651011469064 | 2.927892057493617 | 542.745838072904 | 0.271098936% | PASS |

CPU dt signed impulses are Jn = 0.0010453423941716134 N s, Jt = -9.617757615766688e-5 N s, and Ly = 2.3977988722892204e-7 N m s. The impulse-weighted arm is 0.0024930955510444916 m. CPU dt/2 gives Jn = 0.0010453537344482257 N s, Jt = -9.617673510519946e-5 N s, Ly = 2.3977776623291203e-7 N m s, and an arm of 0.002493095299717128 m. Each impulse agrees with the corresponding particle momentum change to floating-point precision.

Across the four runs, the maximum absolute reconstruction residual is 2.85e-14 N for the signed normal force, 2.78e-17 N for the tangential spring/radial return, and zero for tau_y + (R - delta) Fx. The completed-step velocity-increment residual is below 4.42e-16 m/s, and the angular-velocity-increment residual is below 1.05e-13 rad/s. These values establish consistency of force/torque serialization, half-step timing, and the measured native state.

Original-dt contact has 934 samples: 929 reach the Coulomb cap and five have spring buildup before saturation. Half-dt contact has 1867 samples: 1857 capped and ten uncapped. The maximum abs(Fx)/abs(Fz) is 0.09200000000000001 for every run, while the average over all contact samples is slightly lower because of that buildup. Therefore the average is not an estimate of a different material friction coefficient. Slip at force assembly stays positive, with minimum above 0.6058 m/s; there is no sticking or slip reversal. Four original-dt and seven half-dt samples retain negative signed normal force in the separation tail.

The half step changes the final spin by about 0.00480095 rad/s. It does not remove the approximately 0.27% nominal-radius analytic discrepancy, as expected from the retained R - delta lever arm. All official tolerances still pass. These are endpoint acceptance results for this case, not a claim of exact equivalence to the full LAMMPS stored trajectory.
