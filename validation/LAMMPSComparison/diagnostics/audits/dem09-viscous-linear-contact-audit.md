# LAMMPS dem09: free fall with linear viscous drag

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

The [official input](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem09-viscous-linear-3d-si.yaml) is pinned at the newer commit already used by dem13, **29dd3fd7642a6fe9e8b4dd99120df4068b221250**. Its Git blob is **c1948eb56cac61f1c6159fc060d7d191c1ab897c**. The source carries the official `granular unvetted` tag. This independently executed input is outside the historical24-input catalog, which remains unchanged.

The canonical input is reference/lammps/dem09-viscous-linear-3d-si.yaml (retained locally; not included). The source manifest (retained locally; not included) retains nine live-fetched, Git-blob-verified LAMMPS files: input, driver, shared trajectory runner, analytic checker, sphere-mass definition, gravity, velocity-Verlet integration and viscous-force implementation. Fifteen current native FunDEM source snapshots separately document force modules, gravity, mass and CPU/CUDA integration.

## Exact physical mapping

One sphere has diameter1e-4 m, density2000 kg/m3, initial center (0,0,0.008) m, zero velocity/spin, gravity (0,0,-9.81) m/s2 and linear drag F_drag=-gamma*v, gamma=1e-7 kg/s. The official `pair_style none` means there are no contacts, walls, spring forces or material-pair coefficients. Native material elasticity/restitution values are inactive schema placeholders and do not calibrate this force law.

The sphere mass is m=4*pi*rho*R^3/3=1.0471975511965978e-9 kg. FunDEM's `eventParticleDamping` translational parameter has units1/s and applies F=-dampingRate*m*v. Therefore setting dampingRate=gamma/m=95.4929658551372 1/s is an exact algebraic unit conversion. Rotational drag is zero, and no hydrodynamic/buoyancy module is added. No coefficient is fitted and no solver source changes are made.

In `nativeEventSolver.cpp`, CPU damping calls `particle.addForce((-rate/inverseMass)*velocity)`. `nativeEventKernel.cu` adds the identical value to the device force accumulator. `SphereDEM.cpp` clears forces, assembles contact forces and then adds external modules. `motionIntegration.h` adds gravity separately during velocity updates. Thus the trajectory force channel contains the actual drag contribution and excludes gravity. It is not a normal contact force. The empty contact CSV and absent contact/material-pair configuration verify that distinction.

## Official criterion and supplemental transient references

The original timestep is1e-4 s. Two run segments500 and300 steps give checkpoints0.05 and0.08 s. The official analytic model `terminal_velocity_linear` checks the final downward speed against the **asymptotic terminal speed**, v_inf=m*g/gamma=0.1027300797723863 m/s, with relative tolerance0.005 (0.5%). It does not compare against the finite-time exponential speed. The input's epsilon1e-10 applies to same-code regression states, not cross-code acceptance.

The exact continuous solution is used only as a supplementary transient reference, with lambda=gamma/m and tau=1/lambda:

    v_z(t) = -v_inf*(1-exp(-lambda*t))
    z(t) = z0-v_inf*[t-(1-exp(-lambda*t))/lambda]

At0.08 s the physical continuous solution still differs from the infinite-time terminal speed by exp(-lambda*T)=0.0004810991 relative. Refinement can reduce error to the finite-time exponential solution while slightly increasing error to the asymptotic target. Those are different comparisons and must remain separate.

## Native timing and force/impulse audits

FunDEM applies drag between the drift and final velocity half-kick. Saved velocity is after the final half-kick, so for steps n>=1 the force-assembly velocity is

    v_assembly = v_saved-(dt/2)*(F_drag/m+g_vector)

Use v_assembly when checking F_drag=-gamma*v and inferring gamma; using v_saved directly introduces a discretization-dependent apparent coefficient. The initial sample is special: it represents the specified zero velocity and zero drag before either half-kick, and must not be reconstructed by undoing a nonexistent final kick.

With positive downward assembly speed u_n, q=1-lambda*dt and initial v=0, the native discrete recurrence has an independent closed form:

    u_n = v_inf+(g*dt/2-v_inf)*q^(n-1), n>=1
    v_down,n = (1-lambda*dt/2)*u_n+g*dt/2
    z_n = z0-dt*[n*v_inf+(g*dt/2-v_inf)*(1-q^n)/(1-q)]

This exactly describes the declared integrator and is not fitted to the trace. Supplemental checks compare this recurrence, the continuous transient, direct force law, force-plus-gravity impulse, drift and half-kicks, mass/inertia and zero transverse/spin/torque channels. Mechanical energy decreases through viscous dissipation; conservation of particle momentum or mechanical energy is not an acceptance condition.

## Execution and retention

CPU/GPU runs use dt=1e-4 and5e-5 s, both ending at0.08 s. Every integration step is retained as native output and CSV:801/801/1601/1601 initial-inclusive frames, **4804 native frames total**. All four projects, raw result directories, complete native ZIP, CSV, source snapshots and English report are retained. Compact artifacts are mirrored to WSL; complete raw result directories and ZIP remain authoritative on Windows. This validates external translational drag and gravity integration; contact restitution, friction, rotational drag and fluid-flow coupling require their own cases.
