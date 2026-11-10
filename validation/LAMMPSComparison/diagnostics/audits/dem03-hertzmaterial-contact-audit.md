# Contact and independent physics audit: dem03-hertzmaterial-wall

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. No solver source, fitted coefficient or official initial state is changed.

## Official reference and criterion

The [official YAML](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem03-hertzmaterial-wall-3d-si.yaml) is pinned to commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`, Git blob `2b1e846a6377f3218c0bc3756e3c430738d12a33`. This case uses one free sphere and one fixed plane at z=0: R=0.0025 m, density=7000 kg/m^3, initial center (0,0,0.00255) m, initial velocity (0,0,-3.9) m/s, and zero spin. Official material input is E=1e11 Pa, nu=0.25, restitution=1. The wall command is `wall/gran granular hertz/material E en nu tangential linear_nohistory 0.0 0.0 damping coeff_restitution zplane 0.0 NULL`. Thus no tangential elastic or damping/friction force is present. Gravity, extra force, rolling/twisting resistance and global damping are zero.

Official dt=10 ns, segments 1300+838+900=3038 steps, duration=30.38 us. The half-step run uses dt=5 ns and 6076 steps with identical physical duration. The analytic model `hertz_normal_impact` checks zero-based segment 1: cumulative step 2138, time 21.38 us. It compares E0=(1/2)*m*3.9^2 against (2/5)*|F_particle|*overlap, within 1%. YAML epsilon=1e-10 is a same-code regression setting and is not the cross-code criterion. The three stored official endpoint positions/velocities are compared separately. No local LAMMPS executable rerun is claimed.

## Exact material-to-wall mapping

The source Hertz-material normal implementation (retained locally; not included) calls `mix_stiffnessE_wall(Emod,poiss)` for wall contacts. Crucially, the wall mixing function (retained locally; not included), lines 109-112, returns E/[2*(1-nu^2)]. The official plane has infinite dynamic mass, but its elastic convention is the equal-material effective modulus. Replacing it with E/(1-nu^2) would double stiffness and change the reference problem.

FunDEM uses identical finite E=100 GPa, nu=0.25, restitution=1 for the free sphere and fixed plane material. `linearStiffnessEnabled=false` selects the material Hertz branch. Its contact-material resolver (retained locally; not included) computes the same series effective modulus. The exact coefficients are:

- E*=53.333333333333336 GPa; (4/3)E*=71.11111111111111 GPa.
- Effective radius R_eff=R=0.0025 m and effective mass m_eff=m=0.0004581489286485116 kg. The two-equal-sphere R/2 and m/2 formulas do not apply.
- K=(4/3)*E*sqrt(R)=3555555555.555556 N/m^(3/2), Fn=K*delta^(3/2).
- Hertz elastic energy U=(2/5)*K*delta^(5/2)=(2/5)*Fn*delta.
- Initial kinetic energy E0=0.003484222602371931 J; inertia I=(2/5)*m*R^2=1.1453723216212791e-9 kg m^2.

`dampingEnabled=false` selects restitution-derived damping. Both restitution coefficients are exactly 1, so every automatic dashpot is zero. Friction coefficients are zero. FunDEM material-derived tangential/rolling/torsional stiffness defaults can be nonzero, but with zero motion and zero limits these channels remain inactive. This normal-only case does not validate active tangential history, nonzero damping or `limitDampingWithinCoulombDisk`.

## Independent continuum and integration checks

Energy integration gives delta_max=[5*m*V^2/(4*K)]^(2/5)=22.680643558230896 um and Fmax=384.05244029192096 N. Continuum contact begins at 12.820512820512854 us, reaches peak compression at 21.37889427557118 us and separates at 29.937275730629513 us. The official sampled 21.38 us criterion is retained exactly. Peak compression, peak force, and maximum full-trace relative mechanical-energy error use a separately declared supplemental 1% tolerance.

The 30.38 us endpoint follows separation; its force/overlap/elastic energy must be zero. The outgoing-speed ratio to 3.9 m/s is a supplemental restitution diagnostic; this YAML does not define an additional restitution acceptance test. The full mechanical energy uses particle translation and rotation plus Hertz elastic energy. No wall kinetic energy is added because its state is fixed.

CPU and CUDA retain velocity-Verlet half-kick/drift/force/half-kick integration. With zero damping, force depends only on current overlap and needs no evaluation-velocity reconstruction. Initial gap means the raw t=0 force is physically zero. Every particle momentum increment is compared against trapezoidal wall force impulse, and every position increment against the drift formula. Particle momentum is not conserved during wall impact: the fixed support receives the equal opposite reaction impulse. The audit never applies an isolated-particle momentum conservation assertion.

Every contact checks point z=0, normal (0,0,1), overlap=R-z, positive signed repulsive force, effective mass/radius, area=pi*R*delta, particle/contact force agreement, all three tangential spring components, inactive elastic-energy channels, zero extra torque and zero particle spin/transverse motion. Native results retain every numerical step including initial state: 3039,3039,6077,6077 frames, total 18232.

## Completed four-run independent result

All four runs pass the official 1% sampled energy check and the separately declared supplemental 1% checks. At dt=10 ns, delta_max=22.680646567542744 um, Fmax=384.05251672717037 N, maximum relative total-energy drift=4.985384483013888e-7, and official sampled PE error=3.317049170892261e-7. At dt=5 ns, the corresponding values are 22.68064391819282 um, 384.052449434657 N, 1.2552065564719328e-7, and 3.967684292227602e-8. Halving dt reduces the maximum energy drift by a factor of 3.9717642146696073.

CPU/GPU position, force and spin traces coincide numerically; maximum velocity differences are 2.220446049250313e-16 m/s at base dt and 4.440892098500626e-16 m/s at half dt. Maximum wall-force impulse versus sphere-momentum-change residual across all runs is 1.474514954580286e-17 kg m/s. Native force-law residual is at most 2.644355845404789e-10 N and native normal-energy residual is at most 3.997670250388552e-15 J. Both endpoints are separated; outgoing speed ratios are 1.000000005095901 (base) and 1.0000000002985259 (half), consistent with elastic restitution within integration error.

## Reproducible evidence

- Source manifest and hashes (retained locally; not included): 11 full official snapshots verified by Git blob, 15 current native snapshots verified by SHA256. YAML, normal material law, wall mixing, wall fix and analytic checker were fetched live at the pin for this case; other same-pin snapshots were reused and rehashed.
- Independent physics script (retained locally; not included) and audit JSON (retained locally; not included).
- Root reporting analyzer (retained locally; not included) and numerical summary (retained locally; not included).
- Official checkpoint comparison (retained locally; not included) and complete energy/force history (retained locally; not included).
- Native frame retention audit (retained locally; not included).

The independent physics script does not import the root analyzer. Raw native archives remain authoritative in Windows; compact inputs, CSV, reports, analysis and audited source snapshots are mirrored into the WSL workspace. Previous case outputs and reports remain unchanged.
