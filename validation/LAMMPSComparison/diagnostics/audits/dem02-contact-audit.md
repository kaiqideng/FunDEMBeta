# dem02 independent contact and analytic audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. Source pin: `f9720ce2d176d615495d1d6c31349b1f71be1c73`.

## Official model and criterion

The official YAML prescribes two equal spheres of radius 0.01 m and density 2800 kg/m3, initially at x = +/- 0.0101 m, moving at +/- 10 m/s. The normal law is `hertz 3.33333e10 0.0`; the second argument is the damping parameter and is zero. The default damping submodel is `viscoelastic`, so the dissipative force is zero. This is an elastic collision and must not be configured with restitution zero.

The declared `xmu=0.35` variable is unused. The actual pair command is `tangential linear_nohistory 0.0 0.0`, specifying zero tangential damping multiplier and zero friction. No gravity or rolling/torsional resistance is applied.

Official runs are 500 + 1008 + 1000 = 2508 steps at 20 ns. The analytic checker uses segment index 1, cumulatively step 1508, time 30.16 us. It computes overlap from particle positions and normal-force magnitude from particle 1, then checks `(2/5)*F*overlap` against initial relative kinetic energy, with 1% relative tolerance. This is the only official analytic pass/fail quantity. Global peak overlap, global peak force, and full-history total energy are supplemental checks requested for FunDEM validation.

## Exact normal-law mapping

LAMMPS: `F = kn*sqrt(R_eff)*delta^(3/2)`. FunDEM: `F=(4/3)*E_eff*sqrt(R_eff)*delta^(3/2)`. Both use `R_eff=R1*R2/(R1+R2)=0.005 m` for two spheres. Therefore `E_eff=(3/4)*kn=2.4999975e10 Pa` is exact. For identical materials with Poisson ratio zero, `E_eff=E/2`, so use `E=4.999995e10 Pa`. These E/nu values are an algebraic representation of the official Hertz force coefficient, not extra independently prescribed material data or calibrated values.

Keep FunDEM's all-linear-stiffness switch false, preserving Hertz. With each sphere's restitution set to 1 and automatic damping selected, beta = -ln(1)/sqrt(pi^2+ln(1)^2)=0, so normal, sliding, rolling, and torsional dashpots are all zero. Friction coefficients are zero. Although Hertz/Mindlin tangential stiffnesses remain internally available, no transverse motion or nonzero tangential force is excited in this one-dimensional frictionless case.

The Hertz conservative energy is `U=(2/5)*F*delta=(2/5)*K*delta^(5/2)`, not the linear-spring formula `(1/2)*F*delta`. The shared FunDEM host/device response uses the correct 0.4 factor.

## Independent analytic values

Let `m=4*pi*r^3*rho/3`, `m_red=m/2`, `V=20 m/s`, and `K=kn*sqrt(R_eff)`.

- Sphere mass: 0.011728612573401897 kg.
- Reduced mass: 0.005864306286700948 kg.
- Hertz K: 2.3570202469325542e9 N/m^(3/2).
- Initial total and relative kinetic energy: 1.1728612573401898 J.
- Maximum overlap from energy, `(5*m_red*V^2/(4*K))^(2/5)`: 0.00027411281173443945 m.
- Peak normal force, `K*delta_max^(3/2)`: 10696.884705232753 N.
- First touch from initial gap / relative speed: 10 us.
- Compression duration, `(delta_max/V)*(2/5)*B(2/5,1/2)`: 20.16973591208614 us.
- Peak compression time: 30.16973591208614 us.
- Full separation time: 50.339471824172276 us.
- Official check time: 30.16 us, just before the continuous peak.
- Official endpoint: 50.16 us, about 0.17947 us BEFORE full separation.

The official endpoint therefore still contains a small positive overlap and elastic energy. It must not be described as completed free-flight restitution; the remaining PE belongs in the energy balance. A precise peak from all saved steps must be distinguished from the official 30.16 us checkpoint.

For velocity Verlet, positions and force correspond to the same saved full time; saved velocities include the final half kick. Force-vs-overlap and conservative total energy should be evaluated on these synchronized states. A force impulse computed using trapezoidal time integration should match `m*(v-v0)`; equal-and-opposite particle forces should cancel total momentum at every step.

## Audited source evidence

All seven official snapshots were fetched live from GitHub at the pin and independently checked using Git blob hashes before copying byte-identical local snapshots. Manifest: source-manifest.json (retained locally; not included). Ten read-only FunDEM snapshots are recorded there with SHA-256 hashes.

- [Official YAML](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem02-hertz-twosphere-3d-si.yaml).
- `gran_sub_mod_normal.cpp:179-191`: second Hertz coefficient is damping, normal force `k*a*delta`.
- `granular_model.cpp:216`: default viscoelastic damping; `409-447`: effective radius and `dR=delta*R_eff`.
- `gran_sub_mod_damping.cpp`: viscoelastic prefactor is `damp*m_eff*a`.
- `test_analytic_models.cpp:156-196`: official peak energy criterion, 0.4 force-overlap product.
- `test_dem_common.cpp`: cumulative segment driver and optional analytic invocation.
- FunDEM `contactMaterialFunctions.h:108-143`: Hertz material mixing, all damping paths, 0.4 Hertz elastic energy.
- FunDEM `contactDetection.h:193-223`: sphere-sphere distance, overlap and effective radius.

FunDEM `SphereDEM.cpp:196-217` executes CPU/GPU half kick, drift, interaction force assembly, final half kick. `motionIntegration.h` implements the corresponding shared integration formulas.

No solver implementation or historical case data was modified.

## Native-run independent review

All four completed native CSV pairs were independently audited with `audit_dem02_physics.py`; every official criterion passes. CPU/GPU values coincide to the listed precision. Base dt: overlap 0.0002741127562990858 m, force 10696.88146029747 N, maximum relative total-energy drift 3.562655536137882e-7, official checkpoint energy error 5.055887846161511e-7. Half dt: overlap 0.0002741128193712479 m, force 10696.885152257066 N, maximum relative total-energy drift 8.944617535533439e-8, official checkpoint energy error 7.181148652390945e-7. The global sampled peak changes from step 1508 (30.16 us) to step 3017 (30.17 us); the official check remains 30.16 us for both timesteps. The full-trace energy drift improves by about four times, while the fixed-time peak-energy-only error does not improve monotonically because the official check is slightly before the true peak and omits its small kinetic energy.

Force-law relative residual is 3.41e-16 or less and the saved normal elastic-energy residual is 6.67e-16 J or less. The largest force-impulse/momentum residual is 1.61e-15 kg m/s. Total momentum, action/reaction force residual, all transverse velocities, and spin remain zero. The endpoint overlap remains positive, about 3.5894e-6 m, with about 2.3013e-5 J of normal elastic energy. Root analysis code was independently reviewed for source constants, state synchronization, exact Hertz mapping, energy counting, checkpoint timing and all four variant schemas; its values agree with this separate calculation.
