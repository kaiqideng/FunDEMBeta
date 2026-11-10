# DEM13 SDS twisting-resistance contact audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. Scope: independent official-source and native-model review. No solver source, old benchmark input, or historical source pin was modified.

## Correct official scenario and provenance

The earlier DEM06 spin-impact scenarios spin about an axis tangent to the contact plane and exercise sliding friction; they contain no twisting sub-model. The actual official twisting-resistance case is [dem13-twist-sds-3d-si.yaml](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem13-twist-sds-3d-si.yaml), commit 29dd3fd7642a6fe9e8b4dd99120df4068b221250, Git blob a6b317c73f34232186ceb3207a9422c02686f657. This is a separate, newer source pin; the previous 24-input historical catalog remains unchanged. The official YAML records LAMMPS 4 Jul 2026, generation date 28 Jul 2026, and tags granular unvetted. Its origin and acceptance rule are official, but the unvetted tag must be disclosed.

Eleven live official source files were saved with verified Git-blob hashes and SHA256, together with eleven native source snapshots and SHA256. Their source manifests and snapshots are retained in the original campaign archive, not distributed in this repository.

## Official input and acceptance rule

| Quantity | Official value |
|---|---:|
| Sphere diameter / density | 0.005 m / 2700 kg/m3 |
| Initial center | (0, 0, 0.0024999998) m |
| Initial velocity | zero |
| Initial angular velocity | (0, 0, 10) rad/s |
| Fixed floor | z = 0 |
| Gravity magnitude | 1 m/s2, down |
| Hooke normal stiffness | 1e6 N/m |
| Restitution input | 1.0 |
| Sliding history stiffness | 2.8e5 N/m |
| Sliding damping multiplier / friction | 0 / 0.5 |
| Twisting SDS stiffness ktw | 1 N m/rad |
| Twisting SDS dashpot gtw | 1e-6 N m s/rad |
| Twisting SDS friction parameter mut | 1e-4 m |
| Original time step | 1e-6 s |
| Run segments | 30000 + 20000 steps |
| Total duration | 0.05 s |

The official check in test_analytic_models.cpp:515-528 is twist_decay, analytic_segment = 0: after the first 30000 steps, at 0.03 s. It checks only omega_z against

omega_z(t) = omega0 - [5 mut g/(2 R^2)] t = 10 - 40 t [rad/s].

The target is 8.8 rad/s at 30 ms, with relative tolerance 0.01. The stored finite-step LAMMPS value is 8.7999300950861805 rad/s. At 50 ms the reference omega_z is 8.0000697292061798 rad/s. The YAML epsilon = 1e-10 is a separate same-code trajectory regression tolerance and is not the cross-code analytic acceptance tolerance. Report torque-law and angular-impulse checks as additional independent checks, not extra official pass criteria.

## Exact units and parameter conversion

For R = 0.0025 m, sphere mass m = (4 pi/3) rho R^3 = 0.0001767145867644259 kg; solid-sphere I = (2/5) m R^2 = 4.417864669110648e-10 kg m2. The fixed floor has no dynamics. The sphere-plane effective radius is nominal R, as assigned in contactSearch.cpp:287, not R minus overlap. The torsional resistance in contactFunctions.h:295-305 uses D = 2 R = 0.005 m twice: the spring pseudo-velocity is D omega_z and the physical torque is D times the pseudo-force.

LAMMPS SDS twisting inputs act directly on angular displacement: M_trial = -ktw theta - gtw omega_z and |M| <= mut |Fn|. Consequently mut carries units of length. FunDEM uses length-scaled spring and dashpot coefficients with a dimensionless torsional friction coefficient. Exact conversion is

- FunDEM torsional stiffness = ktw/D^2 = 40000 N/m.
- FunDEM torsional dashpot = gtw/D^2 = 0.04 N s/m.
- FunDEM torsional friction = mut/D = 0.02 (dimensionless).

Using mut/R = 0.04 would double the torque cap and be incorrect. Do not confuse D with the Hertz/Marshall contact-patch radius sqrt(R delta). This selected case is SDS, so no contact-patch-radius conversion is required.

Use explicit stiffness mode with normal 1e6, sliding 2.8e5, rolling 0, torsional 40000 N/m. Use explicit damping with normal, sliding and rolling zero and torsional 0.04 N s/m. Both materials retain restitution 1.0. Automatic restitution damping alone would set all dashpots to zero here and omit the official nonzero twisting dashpot. Young modulus and Poisson ratio required by the project schema are inactive placeholders in this explicit-linear case, not official or fitted material parameters. Sliding friction is 0.5; rolling friction and stiffness are zero. Cohesion is zero.

## Torque law and scope of equivalence

LAMMPS gran_sub_mod_twisting.cpp SDS accumulates a scalar angular history, computes the sum of elastic and viscous torque, truncates that sum to mut |Fn|, and adjusts angular history using both the cap and dashpot. FunDEM integrateTorsionalSpringForce accumulates a normal-axis vector length history, applies radial return to its elastic pseudo-force, then adds as much dashpot force as remains admissible inside the Coulomb disk; that second step does not change the elastic history. These algorithms do not establish general equivalence for unsaturated or reversing torsional histories.

For this benchmark omega_z stays positive near 8 to 10 rad/s. The dashpot torque alone is approximately 8e-6 to 1e-5 N m, hundreds of times above the approximately 2e-8 N m friction cap; the elastic increment at each nonzero time step also reaches the cap. The resulting law in both solvers is therefore the same saturated resistance M_z = -mut |Fn|. This benchmark validates the torque cap, its sign and units, spin decay, and torque-to-angular-momentum consistency. It does not separately identify the unsaturated twisting spring or dashpot response. Unlike the previous oblique-impact case, this is intrinsic normal-axis contact torque; the ordinary tangential-force arm R-delta contributes no torque here because the tangential force is zero.

LAMMPS normal-model Fncrit is abs(Fntot), agreeing with FunDEM's absolute normal-force contract. With restitution 1.0 there is no normal damping. Pure normal-axis spin produces no sliding or rolling velocity at the floor contact; force is along +z and torque along -z throughout this test.

## Normal load oscillation and independent reference

The initial overlap is delta0 = R-z0 approximately 2e-10 m, giving Fn0 approximately 2e-4 N, whereas mg = 1.767145867644259e-4 N. Thus the normal load oscillates about mg by about 13.177%, because the normal contact is undamped. A nearly linear angular-velocity envelope does not imply a constant instantaneous torque.

For a continuous Hooke oscillator with no loss of contact, let wn = sqrt(kn/m), A = kn delta0 - mg. Then

Fn(t) = mg + A cos(wn t),

omega_z(t) = omega0 - (mut/I) [mg t + A sin(wn t)/wn].

This provides an additional continuous-model reference. The official target is the constant-mean-load envelope; tiny non-monotone deviations from it with step refinement are possible because the oscillator phase changes. Report measured step-halving differences without claiming monotone convergence unless observed.

## Native CSV and integration audit conventions

Per-step trajectory output stores completed velocities after the final velocity-Verlet half-kick. The stored contact forces and torques are the force assembly used immediately before that kick. For law reconstruction at step n > 0, use omega_eval = omega_saved - dt M_saved/(2 I). Contact history uses that evaluation velocity.

The saved step-zero project/trajectory frame precedes initial force assembly and has zero displayed force/torque and no contact row. For derived impulse arrays only, reconstruct Fn0 = kn (R-z0) and M0 = -mut Fn0 from the unchanged initial state. The first integration step uses these initialized forces. Do not modify the retained native CSV. The per-particle force_z column is contact force and excludes gravity. Thus verify

m [vz(t)-vz0] = integral Fn dt - mg t,

I [omega_z(t)-omega0] = integral M_z dt,

using trapezoidal per-step integration, including the analytically reconstructed initial force/torque. Compare contact torque_z to the particle torque_z; no tangential moment or rolling/torsional decomposition ambiguity is present in this pure-twist case. The saved torsional spring should satisfy D*k_tors*spring_z = mut*Fn after each nonzero step while capped. Report all force, torque, momentum and history residuals in dem13-physics-audit.json when the native outputs complete.

## Completed independent native checks

All four completed native trajectories were independently checked, including project parameters, all 300000 contact rows, and all 300004 trajectory rows. All four pass the official 30 ms / 8.8 rad/s / 1% criterion. Native result retention is checked separately by archive-retention-audit.json.

| Run | omega_z at 30 ms (rad/s) | Relative error (%) | omega_z at 50 ms (rad/s) | Max torque-cap residual (N m) |
|---|---:|---:|---:|---:|
| cpu | 8.799930095273 | 0.000794372 | 8.000069706515 | 3.309e-24 |
| gpu | 8.799930095273 | 0.000794372 | 8.000069706515 | 3.309e-24 |
| cpu-dt-half | 8.799934025545 | 0.000749710 | 8.000058919050 | 3.309e-24 |
| gpu-dt-half | 8.799934025545 | 0.000749710 | 8.000058919050 | 3.309e-24 |

The largest angular-impulse residual is 2.9913e-22 N m s. The inferred LAMMPS twisting-friction length is 0.0001 m at every contact sample, and the effective native dimensionless coefficient is 0.02. The minimum dashpot torque divided by the cap exceeds 400, confirming saturation throughout the saved contact trajectory. Tangential force and nonaxial torque/spin remain zero.

The continuous oscillating-load solution has maximum full-trajectory spin errors of 6.0108237e-5 rad/s at the original step and 1.5498685e-5 rad/s at half step, a ratio of 0.2578463. An exact discrete velocity-Verlet diagnostic uses theta=acos(1-(wn dt)^2/2), Fn_n=mg+A cos(n theta), and omega_n=omega0-(mut/I)[mg n dt+A dt sin(n theta)/(2 tan(theta/2))]. Maximum native spin differences from this discrete closed form are 2.7548452e-10 and 8.1859852e-11 rad/s, respectively. The continuous-time load-phase difference accumulated over many undamped cycles is consistent with finite-step oscillator dispersion; instantaneous torque and force laws remain satisfied.
