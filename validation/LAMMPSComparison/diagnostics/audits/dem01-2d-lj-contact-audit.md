# Independent contact-law and dimensional-mapping audit: dem01-hooke-2d-lj

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

## Official reference

The exact LAMMPS input is pinned to commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`, Git blob `9bc3f60cd4c5059210c74c24a8bb40a5b822ed17`. [Official YAML](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem01-hooke-2d-lj.yaml). The YAML and relevant mass/inertia sources were independently live fetched during this verification; thirteen full pinned LAMMPS snapshots and fifteen current FunDEM source snapshots are retained and hash verified in the source manifest (retained locally; not included).

The official input selects `units lj`, `dimension 2`, `atom_style sphere`, ordinary `set ... diameter 1 density 1`, `fix nve/sphere disc`, and `fix enforce2d`. Two particles begin at x=-/+0.55 with incoming vx=+/-1. Radius is0.5, normal Hooke stiffness10000, restitution0.5; linear-no-history tangential stiffness/damping and friction are zero. Gravity, initial spin and external forces are zero. The pair style is granular; this input contains no Lennard-Jones pair potential.

## Exact active-equation mapping and limits

Reduced-unit scales L*=1m, M*=1kg, T*=1s give energy scale1J, force scale1N, stiffness scale1N/m, density scale1kg/m3 and damping scale1N s/m. All official numerical values are preserved in the native project, including radius0.5m and density1kg/m3. This is an explicit dimensional realization of the official reduced-unit equations, not a claim that these numbers describe a particular real material.

Critically, the official mass is **not** pi*R^2*density. In pinned `set.cpp`, lines1295-1300 enable a disk-mass flag only for the separate keyword `density/disc`; lines1343-1345 use sphere-volume mass for ordinary `density`. Thus m=4*pi*R^3*rho/3=pi/6 and reduced mass=pi/12. The native sphere mass matches directly without a density adjustment.

Pinned `fix_nve_sphere.cpp`, lines47 and60-64, changes only the inertia prefactor from0.4 to0.5 with the `disc` option; translation still divides by the same rmass. The official disk inertia is0.5*m*R^2 and native sphere inertia is0.4*m*R^2. These differ and are explicitly retained as a scope limitation. In this central, frictionless collision, torque, spin, transverse velocity and out-of-plane motion remain zero, so the different inertias multiply inactive rotational equations. Native three-dimensional spheres therefore reproduce the active one-dimensional translational equations exactly. This does not validate a general two-dimensional solver, disk rotation, or a Lennard-Jones force model.

The unbounded native domain also reproduces the official trajectory over its finite interval: the particles remain well inside the official x/y box, and z=0 remains invariant. No periodic image or boundary contact participates.

## Damping and force timing

LAMMPS `hooke Kn en ... damping coeff_restitution` and FunDEM linear automatic damping yield the same coefficient:

    beta = -ln(e)/sqrt(pi^2+ln(e)^2)
    cn = 2*beta*sqrt(m_reduced*Kn)
    Fn = Kn*delta - cn*vn_evaluation

Native `linearStiffnessEnabled=true` with Kn=10000 and Kt=Kr=KT=0, `dampingEnabled=false` (manual override disabled), identical material restitution0.5 and all stored manual dashpots0 preserve automatic normal damping. All other friction coefficients are zero. Young modulus/Poisson ratio are inactive linear-branch placeholders. Unlike the previous nonzero-Kt history case, no explicit normal damping conversion is configured here.

Native saved velocities include the final half-kick. Before auditing force or inferring cn, reconstruct force-assembly velocity as `v_evaluation=v_saved-dt*F_saved/(2m)`. Contact master/slave orientation is resolved from each stored contact. Initial contact force is genuinely zero because the initial surface gap is0.1. Positive overlap retains a signed negative viscous unloading tail; force is not clipped at zero. This is a damping convention rather than adhesion.

## Time and acceptance

The official timestep is0.0001 reduced time units and segments are947+500 steps, so checkpoints are0.0947 and0.1447. The official `collision_restitution` analytic checker runs at zero-based segment1, the final checkpoint. Its two criteria are:

- Relative restitution: `e_measured=-(vx_left-vx_right)/(2*vin)`; `abs(e_measured-0.5)/0.5 <=0.02`.
- Momentum: `abs(m_left*vx_left+m_right*vx_right) <=0.02*(m_left+m_right)*vin`.

The half-step comparison uses0.00005 and2894 steps at the same end time. Each base run retains1448 native frames including the initial state, and each half-step run2895; four runs total8686. Every native state is also retained in CSV. YAML epsilon1e-10 is a same-code regression tolerance, not cross-code acceptance. No local LAMMPS execution is claimed; official saved states and the pinned analytic checker are the reference.

The continuous relative oscillator predicts first contact at0.05, geometric separation after pi/omega_d of contact, final separated velocities-/+0.5, restitution0.5, final system kinetic energy fraction e^2=0.25 and energy loss fraction0.75. Energy conservation is not an acceptance criterion for this dissipative case. Peak compression and viscous-power integration are supplemental diagnostics.

## Evidence

- Independent executable audit (retained locally; not included) and independent numeric audit (retained locally; not included).
- Pinned source snapshot generator (retained locally; not included).
- Root numeric analysis (retained locally; not included), native projects and raw results (retained locally; not included), [English report](../../cases/dem01-hooke-2d-lj/REPORT.md).

## Independently verified numerical outcome

All four actual native runs pass the official2% restitution and momentum checks. At dt=0.0001 the restitution is0.5004258213318077 on both CPU and GPU, with0.0851643% relative error to0.5. At dt=0.00005 it is0.501496394606861 on CPU and0.5014963946068609 on GPU, with0.299279% relative error. The half-step error is larger, although both grids pass. No monotonic restitution convergence is claimed from two timesteps. The finite-step contact begins at0.0501 and0.05005 respectively; both retain their last positive-overlap force at0.0664 before geometric contact loss on different following grid times. These sampled onset/separation events and the discontinuous signed damping force are consistent with phase-sensitive, nonmonotonic restitution errors.

Independent half-kick force reconstruction recovers cn=22.04795889476057 N s/m. The largest force-law residual is1.42e-14N and largest per-particle force-impulse/momentum closure residual is9.99e-16kg m/s. All total linear/angular momenta and action-reaction force residuals are zero in the symmetric traces. All non-normal force, sliding history, extra contact energy, resistance torque and spin channels are exactly zero. The base runs retain23 negative signed-force samples and the half-step runs45; signed unloading is preserved rather than clipped.

Initial system kinetic energy is0.5235987755982988J, with ideal separated energy0.1308996938995747J. Measured final energy fractions are0.250426002656 at the original step and0.251498633804 at the half step. Maximum compression changes from0.00760326137153m to0.00759709966504m, approaching the continuous0.00759114739482m. The final states are separated. All8686 native frames, every-step two-particle/contact CSV and four loadable configurations are preserved; the two methods independently read the new four runs rather than substituting an earlier case.

