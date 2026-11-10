# Contact and independent physics audit: dem02-hertzmaterial-twosphere

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. No solver source or fitted parameter is changed for this case.

## Pinned official input and acceptance

The [official YAML](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem02-hertzmaterial-twosphere-3d-si.yaml) is pinned at commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`, Git blob `7a5d55b4492162756a80aad7c965316a8068a0ec`. Its eight official source snapshots were fetched live with the GitHub connector and verified byte-for-byte by Git blob SHA; ten current native source snapshots are retained by SHA256. See the source manifest (retained locally; not included).

Both spheres remain free. Radius=0.01 m, density=2800 kg/m^3, initial centers x=-/+0.0101 m, velocities vx=+/-10 m/s, zero transverse motion/spin/gravity/external force. The official material data are E=4.8e10 Pa and nu=0.20. `pair_coeff` uses `hertz/material E 0.0 nu tangential mindlin NULL 0.0 0.35`: the second normal argument is a zero damping prefactor, not restitution zero. Sliding friction 0.35 is ACTUALLY supplied to Mindlin, unlike the earlier `dem02-hertz-twosphere` file where the similarly named variable is unused. It is retained in FunDEM.

The official segments are 500+1008+1000=2508 steps at dt=20 ns (50.16 us). Analytic model `hertz_normal_impact`, zero-based segment 1, checks cumulative step 1508 (30.16 us). It compares initial relative kinetic energy with `(2/5)*|F_particle|*overlap` at that sampled instant with 1% relative tolerance. YAML regression epsilon=1e-10 is not used as the cross-code acceptance criterion. The stored YAML states at all three segment endpoints are compared separately. No local LAMMPS executable is claimed to have been rerun.

## Direct material mapping and force/energy laws

LAMMPS `gran_sub_mod.cpp` supplies the reciprocal-series E/G mixing formulas; `gran_sub_mod_normal.cpp` lines 198-236 implements `hertz/material`: for pair contacts it computes `k=(4/3)*mix_stiffnessE(E,E,nu,nu)`, then inherits `F=k*sqrt(Reff*delta)*delta` from Hertz. FunDEM `contactMaterialFunctions.h` lines 105-117 computes the same series effective modulus and secant normal stiffness. The exact mapping is therefore:

- E*=E/[2(1-nu^2)]=25e9 Pa, Reff=0.005 m;
- K=(4/3)*E*sqrt(Reff)=2357022603.955158 N/m^(3/2);
- Fn=K*delta^(3/2), U=(2/5)*K*delta^(5/2)=(2/5)*Fn*delta;
- mass=0.011728612573401897 kg per sphere, mred=0.005864306286700948 kg;
- initial total kinetic energy=1.1728612573401898 J.

FunDEM directly uses E=48 GPa and nu=0.20, with `linearStiffnessEnabled=false`. No substitute elastic constants or fitted stiffness are needed. `dampingEnabled=false` selects restitution-derived damping. Both restitution coefficients are 1, hence the dissipation factor and all four automatic damping coefficients are exactly zero. Explicit zero coefficient fields are inactive. Rolling/torsional friction is zero, while their Hertz/Mindlin stiffness defaults need not be zero: zero motion and friction leave them inactive.

LAMMPS `mindlin NULL` obtains tangential stiffness from material properties (8G*sqrt(Reff*delta)). FunDEM's nonlinear branch provides the same stiffness form. This central collision never activates tangential velocity, spring history, force, spin, rolling or torsional torque. Thus the run verifies the material-to-normal-Hertz branch and inactive tangential channels; it does not establish equivalence of active friction history, tangential damping, or the damping-within-Coulomb-disk algorithm.

## Independent continuum and integration checks

For Vrel=20 m/s, continuum energy gives delta_max=[5*mred*Vrel^2/(4*K)]^(2/5)=0.00027411270208928186 m and Fmax=10696.888983989631 N. Contact begins at 10 us, continuum peak is 30.169727844189347 us and separation 50.339455688378694 us. These are supplemental diagnostics, separate from the official sampled 30.16 us criterion.

The fixed official endpoint 50.16 us precedes full separation; its remaining Hertz elastic energy must be included in the energy balance. Outgoing speed ratio at that instant is not reported as separated-body restitution. The half-step runs use dt=10 ns for 5016 steps, preserving the same physical duration and checkpoints. Native binary results and CSV retain every numerical step including initial state: 2509,2509,5017,5017 frames, 15052 total.

FunDEM CPU and CUDA retain velocity-Verlet half-kick/drift/force/half-kick integration. With no damping, force depends only on current overlap, so no force-evaluation velocity correction is needed for the Hertz law. Independently, every complete-step momentum change is checked against trapezoidal force impulse. At t=0 the gap gives zero force; both-body action/reaction, total momentum, zero transverse state/spin and contact channels are checked directly.

## Reproducible evidence

- Independent physics script (retained locally; not included)
- Independent audit JSON (retained locally; not included)
- Root independent analyzer (retained locally; not included)
- Numerical summary (retained locally; not included)
- Official checkpoint comparisons (retained locally; not included)
- Full energy/force history (retained locally; not included)
- Complete native frame retention audit (retained locally; not included)

The independent script checks all four CSV trajectories without importing the reporting analyzer. It reconstructs overlap, Fn, normal elastic energy, kinetic energy of both bodies, total energy, continuum maxima, official peak-energy acceptance, momentum/impulse, effective mass/radius/area, signed repulsive force, and all zero tangential/rolling/torsional channels. It separately compares CPU/GPU traces and time-step energy convergence. Full raw results remain authoritative under the Windows artifact root; the WSL mirror holds compact configurations, CSV, analysis, reports and source snapshots.
