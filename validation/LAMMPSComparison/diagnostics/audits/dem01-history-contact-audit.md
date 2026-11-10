# Independent contact-law audit: dem01-hookehistory-twosphere

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

## Official reference and acceptance

The exact official input is pinned to commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`, Git blob `6900a6efd5cf80f862da2241c2b48e7d4cbaf7f1`. [Official YAML](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem01-hookehistory-twosphere-3d-si.yaml). It was independently live fetched during this verification. Nine full pinned source snapshots are byte verified against their Git blob hashes; fifteen current FunDEM source snapshots are retained with SHA256 in the source manifest (retained locally; not included).

Two equal free spheres have radius 0.0025 m, density 2700 kg/m3, initial centers x=-/+0.00255 m and incoming vx=+/-2 m/s. All initial spin, gravity, imposed force and constraints are zero. Official parameters are normal Hooke stiffness Kn=1000000 N/m, linear-history tangential stiffness Kt=280000 N/m, tangential damping multiplier 0, sliding friction mu=0 and normal restitution e=0.5. No rolling or twisting resistance is active.

At the official timestep 30 ns, 2068+727 steps finish at 83.85 us. The `collision_restitution` analytic checker runs at zero-based segment 1. Its two independent criteria are:

- `e_measured = -(v_left,x - v_right,x)/(2*vin)` and `abs(e_measured - 0.5)/0.5 <= 0.02`.
- `abs(m_left*v_left,x + m_right*v_right,x) <= 0.02*(m_left+m_right)*vin`, with vin=2 m/s.

The half-step comparison uses 15 ns and 5590 steps at the same physical end time. Every integration state is retained natively and in CSV: 2796 frames per base run, 5591 per half-step run, 16774 across CPU and CUDA. YAML epsilon 1e-10 is a same-code regression setting, not cross-code acceptance. No local LAMMPS execution is claimed; pinned saved states and the pinned analytic checker form the reference.

## Exact damping mapping

LAMMPS `hooke Kn en ... damping coeff_restitution` interprets the second normal coefficient as e. Its Hooke restitution prefactor and force routine yield:

    beta = -ln(e)/sqrt(pi^2+ln(e)^2)
    cn = 2*beta*sqrt(m_reduced*Kn)
    m = 4*pi*R^3*rho/3 = 0.00017671458676442588 kg
    m_reduced = m/2 = 0.00008835729338221294 kg
    cn = 4.050468687151587 N s/m

The official history model reads `Kt=280000`, `xt=0`, `mu=0`; its tangential dashpot is `ct = xt*normal_damping_prefactor = 0`. FunDEM automatic damping would derive a nonzero tangential coefficient from this nonzero Kt and e=0.5. Therefore this exact configuration uses `linearStiffnessEnabled=true`, `dampingEnabled=true`, explicit cn from the analytic conversion above, and explicit ct=cr=cT=0. Rolling/twisting stiffnesses and all friction coefficients are zero. Restitution 0.5 is retained in the material but is inactive for actual damping while the explicit override is enabled. Young modulus and Poisson ratio are inactive placeholders in this linear stiffness branch.

This case checks the converted explicit normal dashpot and the zero-friction degeneration of a nonzero-stiffness history configuration. It does not provide a new automatic-damping check or validate active tangential friction, history rotation or the damping-in-Coulomb-disk limiter. The earlier normal Hooke case separately checks automatic normal damping with zero Kt. Its unchanged trajectory is used here only as a supplemental degeneracy comparison.

## Force timing and signed unloading

At positive overlap, signed normal force is `Fn = Kn*delta - cn*vn_evaluation`. Native completed-step velocities include the final half-kick. Reconstruct each force-assembly velocity as `v_evaluation = v_saved - dt*F_saved/(2*m)` before inferring cn. Initial force is genuinely zero because the spheres initially have a gap. Contact master/slave orientation can change, so the independent audit computes each signed force along the stored master normal and compares the corresponding master particle force.

Both contact models retain the signed negative unloading tail; no damping-force clipping is requested. `normal_force_magnitude` is absolute magnitude, while `normal_force_signed` preserves sign. The negative tail at positive overlap is a viscous convention, not adhesion. All tangential spring components, non-normal forces, extra contact energies, torque and spin must remain exactly zero.

## Continuum diagnostics

For relative overlap after first contact, `m_reduced*delta'' + cn*delta' + Kn*delta = 0`. Initial relative speed is 4 m/s and first contact is at 25 us. The underdamped oscillator predicts geometric separation at 55.2407516365 us, peak overlap 27.8916565193 um at 38.0300426838 us, separated velocities -/+1 m/s and restitution 0.5. These are supplemental continuous-model diagnostics; the official checks concern separated restitution and total x momentum at 83.85 us.

Initial system kinetic energy is 0.000706858347058 J. The ideal separated kinetic energy is 0.000176714586764 J, so the ratio is e^2=0.25 and the intended energy loss is 0.000530143760293 J. Energy conservation is not required for this dissipative collision. The independent audit compares geometric overlap, normal spring energy, signed force law, drift/kick timing, force impulse and particle momentum changes, total momentum, angular momentum, and viscous-power quadrature. Finite-step viscous-power quadrature is diagnostic rather than an exact equality requirement.

## Retained evidence

- Independent executable audit (retained locally; not included) and independent numeric audit (retained locally; not included).
- Preexisting-artifact baseline (retained locally; not included): thirteen old manifest records and fifty-one report/reference file hashes frozen before shared metadata changes.
- Root analysis (retained locally; not included), loadable projects and full native output (retained locally; not included), and [English report](../../cases/dem01-hookehistory-twosphere/REPORT.md).

## Independently verified numerical outcome

All four real native runs pass the official restitution and momentum checks. At 30 ns the measured restitution is 0.4996646896338126 on CPU and 0.4996646896338127 on GPU; at 15 ns it is 0.49984130980353825 on both. Relative errors to e=0.5 fall from 0.0670621% to 0.0317380%. Total momentum and angular momentum remain zero in these symmetric traces. Final kinetic-energy fractions are 0.249664802067 and 0.249841334986, approaching the ideal 0.25.

Inferred cn is 4.050468687151587 N s/m in the base runs and 4.050468687151586 N s/m in the half-step runs. The maximum signed-force reconstruction residual is 3.56e-15 N, and the largest per-particle momentum/force-impulse closure residual is 2.93e-18 kg m/s. There are 140 negative-force samples in each base run and 279 in each half-step run, with minima -8.10497386771 N and -8.09233080662 N. All tangential force/history, extra torque, spin and non-normal elastic energy channels remain exactly zero.

The supplemental comparison with the earlier zero-Kt automatic normal Hooke case reads all four old full trajectories. Across 36 numeric trajectory columns, the maximum absolute difference is 7.11e-15, coming from floating-point force arithmetic. The old reference files remain unchanged. This confirms the expected zero-friction normal degeneration without substituting old outputs for the four newly executed native runs. Maximum overlap is a supplemental phase-sensitive metric; monotonic improvement of every diagnostic is not asserted.
