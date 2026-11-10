# PFC linear restitution: selected full-dashpot subcase

> Historical diagnostic/coverage evidence from the 28–29 October 2026 campaign. Capability statements concern that configured source snapshot and application path, not an exhaustive current-core feature list. Migration changes links and artifact labels, not the recorded outcome.

## Official source and exact scope

Primary source: [PFC 7.0, Linear Contact Model: Calibrating the Normal Critical Damping Ratio](https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/linear_contact_model/cmlinear_restitution.html). The complete `restitution.p3dat` appears inline on this official page. Live connector text snapshots and an extracted datafile are retained under reference/pfc/linear-restitution-source (retained locally; not included). The source manifest (retained locally; not included) records SHA256, URLs and the extraction scope. These are web-connector text snapshots, not raw HTML, pinned Git blobs or a local PFC execution.

The original program has two rows of 16 noninteracting balls, with damping modes Md=0 and Md=1. This case selects only original ball 8 in the Md=0 row, beta_n=0.2, initially at (0.75,0.05,1.0) m. It does not claim the entire 32-ball sweep or no-tension Md=1. The selected ball's original physical coordinates are retained. The immobile floor is represented by an SDF plane; other balls and reflect-domain boundaries are inactive for this one-dimensional motion.

## Parameters and force law

The official input defines R=0.025 m, rho=2650 kg/m3, kn=50000 N/m, centroid height 1.0 m, gravity (0,0,-10) m/s2, no initial velocity or spin, and a fixed timestep 1e-5 s followed by `model solve age 1.0`. Thus the actual free-fall clearance is 0.975 m. The linear-model property table explicitly sets defaults ks=0, fric=0, reference gap rgap=0, absolute normal update lin_mode=0 and tangential critical damping dp_sratio=0. The [ball attribute documentation](https://docs.itascacg.com/pfc700/pfc/pfcmodule/doc/manual/ball_manual/ball_commands/cmd_ball.attribute.html) explicitly gives local damping default 0, so no unlisted body damping is introduced.

The [linear model](https://docs.itascacg.com/pfc700/common/contactmodel/linear/doc/manual/cmlinear.html), Eq. 10, specifies full normal dashpot load for Md=0 and reduced mass m_c=m for ball-facet contact. With upward compressive force positive, while overlap delta=R-z>0:

- m=4*pi*rho*R^3/3=0.1734420944169365 kg;
- beta_n=0.2;
- c_n=2*beta_n*sqrt(m*kn)=37.24965443243054 N s/m;
- signed F_n=kn*delta-c_n*v_z, without a no-tension cap;
- normal elastic energy U=kn*delta^2/2;
- zero tangential force, history, contact moments and spin.

FunDEM uses linear contact stiffness kn=50000 and automatic damping. The material input is e0=exp(-beta*pi/sqrt(1-beta^2))=0.526620599330303, which analytically recovers beta=0.2 and the same c_n. The manual-damping switch remains off and all stored manual coefficients remain zero. Shear/rolling/torsional stiffness, friction and damping are zero. No coefficient is fitted and no solver source is modified.

The native normal force must be audited with the force-assembly half-step velocity. `motionIntegration.h` adds gravity separately from the accumulated native force, so saved trajectory `force_*` excludes gravity. Saved velocity is after the final half-kick, giving v_half=v_saved-(dt/2)*(F_native/m+g_vector). The normal force's signed channel is used; its magnitude channel is absolute and cannot replace it. A negative unloading tail is correct for Md=0.

## Published comparison target

The official program finds the first height minimum and then the first maximum of each one-step height history. It computes e_h=sqrt(max(z_first_apex-R,0)/(height-R)). It then computes `err=max(abs(e_h-e_analytic))` across the Md=0 row. The page reports maximum difference below 1% for dt=1e-5 s; the code establishes that this is an absolute restitution-coefficient difference, 0.01, rather than a relative one-percent error. The program contains no automated acceptance assertion. This selected-subcase verification therefore reports whether abs(e_h-e0)<0.01, describing it as the published fixed-dt benchmark target, without presenting it as a preexisting PFC test assertion.

Gravity remains active during contact. Consequently the first-rebound coefficient derived from height differs slightly from the unforced damped-oscillator formula e0 even as dt tends to zero. As a supplemental reference, solve delta''+(c_n/m)*delta'+(kn/m)*delta=g from delta(0)=0 and delta'(0)=sqrt(2*g*(height-R)), locate the first return to zero overlap, then propagate ballistic rebound to the apex. This gravity-inclusive reference is derived independently, does not replace the official comparator and is not fitted to FunDEM output.

## Runs, retention and interpretation

Four native runs use CPU/GPU and dt=1e-5/5e-6 s, ending at 1 s. The native result output interval is 1 ms, with 1001 initial-inclusive frames in every run, totaling 4004 native frames. Every numerical step is retained in trajectory/contact CSV (100001 or 200001 trajectory states). Complete Windows result directories and ZIP are authoritative; compact configs, CSV, reports, references and audits are mirrored to WSL. A second impact may occur before 1 s and must be reported from the trace; the published metric remains the first rebound apex.

Energy is dissipated and sphere momentum changes due to the floor and gravity. Supplemental checks compare total force with contact force plus gravity, impulse with momentum change, spring energy and dashpot force law, mass and radius, inactive transverse/rotational channels, and CPU/GPU traces. They do not assert conservation of sphere momentum or mechanical energy. Retention/hash verification is separate from physical pass/fail.

Only FunDEM's native CPU/CUDA libraries are run here. There is no local PFC run and no comparison to unpublished raw PFC numerical trajectories. This is an exact active-equation mapping of one documented independent ball, with source-derived analytical comparisons.
