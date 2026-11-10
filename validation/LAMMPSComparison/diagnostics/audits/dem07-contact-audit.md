# dem07 Hooke rolling resistance: independent model audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Audit date: 2026-10-28. This document reviews equations, parameter mapping, and output semantics. It does not claim a simulation pass; run outcomes belong in the independently generated analysis and case report.

## Official source and acceptance

Reference: [dem07-hooke-rolling-3d-si.yaml](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem07-hooke-rolling-3d-si.yaml), commit f9720ce2d176d615495d1d6c31349b1f71be1c73. The YAML records LAMMPS 4 Jul 2026 and a reference-generation time of 21 Jul 2026.

A solid sphere starts at rest in translation and rotates about +y on a fixed horizontal wall. Gravity is 1 m/s^2, rather than 9.81 m/s^2. Sliding friction is disabled. No torsional resistance is defined. This isolates resistance to rotation about a contact-tangent axis; it is not a pure-rolling translation test.

The official run segments are 3000, 4000, and 3000 steps at dt = 1e-6 s. Cumulative checks occur at 3, 7, and 10 ms. The official analytic_model is rolling_decay, analytic_segment = 1 (the second segment, at 7 ms), and analytic_tol = 0.01. The analytic check is:

    omega_y(t) = omega0 - (5 mu_r g / (2 R)) t
    abs(omega_actual(0.007) - 0.65) / 0.65 <= 0.01

See [test_analytic_models.cpp, rolling_decay](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp#L112-L124). The 10 ms value of 0.5 rad/s is an additional diagnostic, not the configured official analytic checkpoint. The YAML epsilon = 1e-10 is for LAMMPS same-code regression states and should not be relabeled as a general cross-code tolerance.

## Direct parameter mapping

| Quantity | LAMMPS | FunDEM setting |
|---|---:|---|
| Diameter | 0.005 m | R = 0.0025 m |
| Density | 2700 kg/m^3 | Sphere material density = 2700 |
| Initial centre | (0, 0, 0.002499982) m | Unchanged |
| Initial velocity | (0, 0, 0) m/s | Unchanged |
| Initial angular velocity | (0, 1, 0) rad/s | Unchanged |
| Normal stiffness | 10000 N/m | normalStiffness = 10000 |
| Sliding stiffness | 0 N/m | slidingStiffness = 0 |
| Rolling SDS stiffness | 10000 N/m | rollingStiffness = 10000 |
| Torsional stiffness | absent | torsionalStiffness = 0 |
| Sliding coefficient | 0 | slidingFrictionCoefficient = 0 |
| Rolling coefficient | 0.05 | rollingFrictionCoefficient = 0.05 |
| Torsional coefficient | absent | torsionalFrictionCoefficient = 0 |
| Normal restitution | 1 | Both materials restitution = 1 |
| Rolling dashpot | 0 | Automatic value = 0 |
| Gravity | (0, 0, -1) m/s^2 | Unchanged |
| Original dt | 1e-6 s | 10000 steps, 1001 archived frames; 10001 trajectory CSV rows |
| Half dt | 5e-7 s | 20000 steps, 1001 archived frames; 20001 trajectory CSV rows |

FunDEM must set linearStiffnessEnabled = true. It selects all four configured sphere-contact stiffnesses together, so kr is not overwritten by Hertz/Mindlin defaults. Set dampingEnabled = false to use restitution-derived automatic dashpots. This field means no explicit dashpot override. With e = 1, beta = -ln(e)/sqrt(pi^2 + ln(e)^2) = 0. The linear formula c = 2 beta sqrt(m_eff k) then gives cn = cs = cr = ctors = 0 exactly. This is a zero-damping configuration, not a nonzero automatic-damping calibration test.

Material elastic E and nu do not control this linear contact law. The wall has zero inverse mass and no prescribed motion. The periodic x/y LAMMPS boundaries have no active effect for a single stationary sphere with no transverse force; this boundary-equivalence scope should be stated if a finite FunDEM domain is used.

## SDS rolling law and sign

LAMMPS uses a length-scaled rolling history, a pseudo force, and a moment:

    v_r = R_eff (omega_i - omega_j) x n
    xi_r <- tangential projection of xi_r + dt v_r
    F_r,trial = -k_r xi_r - gamma_r v_r
    |F_r| <= mu_r |F_n|
    M_r = R_eff n x F_r

For the present gamma_r = 0 and flat, fixed normal, FunDEM implements the same equations. FunDEM first limits the spring force, then admits damping inside the Coulomb disk; LAMMPS limits the combined spring-dashpot force and adjusts history. Those choices are equivalent here because the rolling dashpot is zero. They must not be assumed equivalent for arbitrary nonzero damping or rotating contact normals.

At the wall, n = +z and omega = +omega_y y. Therefore v_r is +x, F_r is -x, and M_r is -y, opposing the initial rotation. Both laws use R_eff = R for sphere-plane contact. Thus, during capped positive rotation:

    M_y = -mu_r R F_n
    k_r R^2 = 0.0625 N m/rad

FunDEM's sphere-level-set contact point lies at the wall and gives the ordinary force moment an R - delta arm. However, its separate rolling resistance law receives effectiveRadius = masterSphere.radius(), and CPU/GPU assemble the rolling torque directly. Since the sliding force is zero and the normal force is collinear with the centre-to-contact vector, that ordinary force moment is zero. The dem05 R - delta sliding-torque issue does not enter this dem07 test.

LAMMPS sources: [SDS history and cap](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_rolling.cpp#L77-L143); [sphere-wall radius](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L403-L410); [rolling velocity and torque](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L518-L532). Full immutable snapshots and Git blob IDs are retained under reference/lammps/dem07-source/.

## Analytical values and normal-load oscillation

For a solid sphere, m = (4/3) pi rho R^3 and I = (2/5) m R^2:

| Quantity | Value |
|---|---:|
| Mass m | 1.767145867644259e-4 kg |
| Isotropic moment of inertia I | 4.4178646691106484e-10 kg m^2 |
| Mean normal load mg | 1.767145867644259e-4 N |
| Mean capped rolling moment mu_r R mg | 2.208932334555324e-8 N m |
| Mean angular deceleration | 50 rad/s^2 |
| Static equilibrium overlap mg / kn | 1.767145867644259e-8 m |
| Initial overlap R - z0 | approximately 1.8e-8 m |
| Normal natural frequency sqrt(kn/m) | 7522.527780636749 rad/s |
| Normal oscillation period | 8.352491995247562e-4 s |

The supplied initial overlap is slightly larger than its equilibrium value. Since the normal damping is zero, the normal load oscillates around mg with about 1.85916% amplitude. This oscillation is present in the supplied official torque states and is not a rolling-friction error. A per-frame cap check must use the real F_n, not mg. Using mg to require every torque sample within 1% would incorrectly reject the intended case.

The approximate mean-load law gives omega_y = 1.0, 0.85, 0.65, and 0.5 rad/s at 0, 3, 7, and 10 ms. The official reference states are:

| Time (ms) | z (m) | vz (m/s) | omega_y (rad/s) | M_y (N m) |
|---:|---:|---:|---:|---:|
| 0 | 0.002499982 | 0 | 1 | 0 (empty initial history) |
| 3 | 0.0024999826039845000 | -1.3471417605498602e-6 | 0.85009282187891866 | -2.1745019375082365e-8 |
| 7 | 0.0024999825690977923 | 1.6832868579600465e-6 | 0.64994130044799026 | -2.1788627759672385e-8 |
| 10 | 0.0024999820048911287 | -4.2486660149128432e-7 | 0.50004670812095975 | -2.2493886089164878e-8 |

The 0 ms row follows the input and zero-history initialization; the other three rows are copied numerically from the official YAML. The physical spring displacement needed to reach the mean cap is mu_r mg / kr = 8.835729338221296e-10 m, reached after approximately 3.5343e-7 s at the initial R omega. This is shorter than either test time step, so the first computed force sample after initialization should already be capped.

## Native frame and integration semantics

Both FunDEM CPU and GPU initialize by assembling force/torque with historyTimeStep = 0. Empty rolling history and zero dashpot yield M_0 = 0 despite nonzero normal load. Each advance performs a half velocity/angular-velocity update, a full position/orientation update, force/history assembly with dt, and the second half velocity/angular-velocity update. A frame captured after step() therefore stores full-step omega_n together with endpoint torque M_n.

For the sphere's isotropic inertia, the independent force-to-motion check is:

    omega_n - omega_(n-1) = dt (M_(n-1) + M_n) / (2 I)

or cumulatively:

    omega_n = omega_0 + dt/I [M_0/2 + sum(M_1 ... M_(n-1)) + M_n/2].

The first step must use half the first nonzero torque, not a full impulse. Check endpoint torque against endpoint F_n. The rolling history increment uses the angular velocity after the first half kick, so checking its uncapped update against the final full-step angular velocity would be a timing error.

Additional diagnostics: x/y positions and velocities stay zero; omega_x and omega_z stay zero; the fixed wall remains fixed; F_n = kn delta; sliding force is zero; rolling history x is positive; M_y is negative while omega_y > 0; M_y/(R F_n) = -0.05 for capped frames. CPU and GPU should agree at matching time steps. Compare the original and half-step curves at common physical times rather than matching step numbers.

## Reviewed FunDEM implementation

Read-only source files under //wsl$/Ubuntuhistorical-archive::

- execution/contactFunctions.h:89-102, 162-193, 271-305: beta, dashpot, capped history, rolling velocity and moment.
- execution/contactMaterialFunctions.h:94-141: the shared linear switch and automatic damping.
- interaction/contactSearch.cpp:249-289: sphere-level-set point, mass, and physical effective radius.
- interaction/contact.cpp:73-87: independent rolling torque and relative spin.
- execution/cpu/particleFunctions.h:109-115: ordinary force moment plus independent rolling torque.
- execution/cuda/contactForceKernel.cu:112-129: shared CPU/GPU rolling-law routine.
- execution/cuda/contactDetectionKernel.cu:477: sphere-level-set effective radius equals master radius.
- solver/SphereDEM.cpp:173-217: zero-history initialization and velocity-Verlet CPU/GPU advances.

Source copies and SHA256 hashes are retained in tools/dem07-contact-source/. They identify the files reviewed at audit time; binary build identity remains the responsibility of the native run manifest. No solver source, physical parameter, or integration method was changed by this audit.

## Independent recomputation of completed CSVs

Read-only recomputation of all four native runs is saved in dem07-independent-check.json, with the script in dem07-independent-check.py. At 7 ms, original-step CPU/GPU omega_y = 0.6499413004459774 rad/s (0.0090307006% from the official mean-load ideal), and half-step CPU/GPU omega_y = 0.6499285591756649 rad/s (0.0109908961%). Both are within 1%. The original step agrees with the official YAML state to 2.0128e-12 rad/s.

The half-step value is slightly farther from the mean-load ideal, so this result must not be described as monotonic convergence toward 0.65. The approximate analytic law omits the real normal-load oscillation and the initial zero-history force timing. CPU/GPU agreement and the force-to-motion identities give stronger direct checks here.

Across all four runs, the maximum absolute residual in M_y/(R F_n) + 0.05 is 6.94e-18; in M_y + R kr xi_rx it is 6.62e-24 N m; Fn - kn delta is exactly zero in the exported doubles. The largest per-step trapezoid angular-impulse residual is 1.103e-16 rad/s, and the cumulative residual is at most 9.55e-15 rad/s. Transverse translation, other spin components, and sliding force remain zero.

The authored time-zero trajectory snapshot has zero stored forces and no exported contact row. Native initialization occurs before the first advance. Therefore perform F_n = kn delta and rolling-cap checks on the active contact rows (steps 1 onward), not the uninitialized authored frame. The initial rolling torque is zero in either representation, and the trapezoid angular-impulse check including time zero passes.
