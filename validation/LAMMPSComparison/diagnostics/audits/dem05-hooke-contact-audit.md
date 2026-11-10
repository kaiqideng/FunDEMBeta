# dem05 Hooke-history sliding-to-rolling source and physics audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 28 October 2026. Source inspection is read-only. Measured pass/fail comes from the four native trajectories and independent audit JSON.

## Official input and acceptance

The source is [LAMMPS dem05-hookehistory-slip-3d-si.yaml](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem05-hookehistory-slip-3d-si.yaml), pinned at f9720ce2d176d615495d1d6c31349b1f71be1c73, Git blob 9d4f2bbc7cf3aa5c10728a2cb2e90cbce024fa91. It labels LAMMPS 4 July 2026 and stored regressions generated 21 July 2026.

One solid sphere: radius 0.0025 m, density 2700 kg/m3, initial center (0,0,0.002499982) m, velocity (1,0,0) m/s and zero spin. Gravity is (0,0,-1) m/s2. Lower fixed plane z=0; upper plane z=1 is inactive. The contact is Hooke kn=10000 N/m, tangential linear_history kt=2800 N/m, friction 0.4, normal restitution 1, tangential damping multiplier 0. Rolling/torsional resistance is absent.

The original timestep is 10 us, cumulative segments 25000+25000+25000 steps, checkpoints 0.25, 0.50 and 0.75 s. Four native runs retain the duration: CPU/GPU at 10 us and CPU/GPU at 5 us. CSV contains every numerical step; native binary results contain the initial frame plus 1 ms output, 751 frames per run and 3004 total.

The [official analytic checker](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp#L82) checks only final vx and omega_y, each with relative tolerance 0.01:

- v_final = 5 u0/7 = 0.7142857142857143 m/s.
- omega_final = v_final/R = 285.7142857142857 rad/s.
- Ideal slip cessation time t*=2u0/(7 mu g)=0.7142857142857143 s is a supplementary diagnostic, not an extra official assertion.

The input epsilon=1e-10 is a same-code stored regression threshold, not the cross-code analytic acceptance tolerance. Compare stored checkpoint states separately. No local LAMMPS executable rerun is asserted.

## Stiffness and automatic zero damping

FunDEM sets linearStiffnessEnabled=true, kn=10000 and kt=2800 directly. E/nu remain inactive schema placeholders; they do not produce these stiffnesses. dampingEnabled=false selects automatic restitution-derived damping, rather than manual dashpot fields. With both materials e=1:

beta = -ln(e)/sqrt(pi^2+ln(e)^2) = 0, and cn=2 beta sqrt(m_eff kn)=0; all translational and rotational automatic dashpots are zero. LAMMPS coeff_restitution gives zero normal damping and its tangential multiplier is zero. Stored explicit dashpot fields are also zero. Set rolling/torsional stiffness/friction to zero, global damping zero, no constraints or external force modules.

Primary formulas: [LAMMPS Hooke normal](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_normal.cpp), [restitution damping](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_damping.cpp), and [linear-history tangential model](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_tangential.cpp).

The shared active laws are Fn=kn delta and Ft=-kt xi, with tangential history increment followed by Coulomb radial return |Ft|<=mu|Fn|. The zero-dashpot limit cannot identify, exercise, or validate the nonzero-damping scaling rule in limitDampingWithinCoulombDisk. Friction still dissipates energy; normal/tangential elastic oscillations can persist because dashpot damping is zero.

## Plane torque-arm conventions

LAMMPS uses nominal sphere radius R for rotational wall slip and friction torque: vt=vx-R omega_y, tau_y=-R Fx. Evidence: [granular_model.cpp wall geometry](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L471) and [wall force application](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/fix_wall_gran.cpp#L577).

FunDEM uses the actual shared plane point at z=0: arm a=R-delta=z_center, vt=vx-a omega_y and tau_y=-a Fx. Source snapshots: FunDEMBeta/src/execution/contactDetection.h, contactFunctions.h and cpu/particleFunctions.h, plus interaction/contact.cpp. Current source and observed force/torque identities establish the convention used in this native run; identity with the installed GUI executable is not asserted.

This Hooke variant is stiff enough that the difference is very small. Sphere mass m=1.767145867644259e-4 kg, I=(2/5)mR^2=4.417864669110648e-10 kg m2. Static overlap mg/kn=1.767145867644259e-8 m (17.671458676 nm), delta/R=7.06858347e-6. This is separate from the previously tested much softer Hertz-Mindlin case; do not alter either result or fit R/friction/stiffness to compensate.

A source-based constant-height diagnostic is v_geo=[m a^2/(I+m a^2)]u0, omega_geo=v_geo/a and t_geo=(u0-v_geo)/(mu g). It can explain a tiny offset but is not the official target. Finite tangential springs produce oscillations after the first slip reversal; a first zero crossing alone does not establish sustained zero slip.

## Independent every-step identities

The script audit_dem05_hooke_physics.py (retained locally; not included) reads all four complete CSV files and project snapshots. It checks input mapping, Hooke normal law, tangential spring force and radial return, Coulomb ratio, true plane point, particle-contact force equality, actual friction torque arm, native extra rolling/torsional torque zero, and conservative spring energies.

Native output velocities include the completed final half-kick. Recover force-assembly values using ve=v_saved-dt F/(2m) and omega_e=omega_saved-dt tau/(2I), then check xi_new=clip(xi_old+dt*(ve-a*omega_e),-mu Fn/kt,+mu Fn/kt). Vertical evaluation velocity would include gravity in the corresponding half-kick.

The raw initial saved frame precedes solver force assembly. Initialization calls assembleForceAndTorque(historyTimeStep=0), so reconstruct Fn0=kn*(R-z0), while Fx0=0 and tau0=0 remain correct for zero initial spring and zero damping. Only derived quadrature arrays use Fn0; raw CSV/native frames are never changed. Particle force excludes gravity. Thus Jx=m(vx-u0), Jz-mgt=m vz and integral(tau_y dt)=I omega_y must close under trapezoidal velocity-Verlet quadrature. Particle momentum itself is not conserved because the plane and gravity are external. Position drift and completed-step kick identities are also tested.

The initial normal load differs slightly from weight and oscillates. The undamped discrete vertical normal-force solution has phase theta=acos(1-(kn/m)dt^2/2), Fn[n]=mg+(Fn0-mg)cos(n theta). This separates ordinary time-discretization phase drift from a contact-law defect. Mechanical energy includes translational and rotational kinetic energy, gravity, normal spring and tangential spring energy; frictional loss is expected.

## Periodic representation and retained evidence

Official x is periodic in [-0.05,0.5), length 0.55 m; y is periodic over 0.1 m. An isolated 0.005 m sphere cannot contact its own image. A wide uniform stationary FunDEM plane with unwrapped x is physically equivalent here. For stored-position comparison, x_wrapped=-0.05+positive_mod(x_unwrapped+0.05,0.55). The final official wrapped x=0.087756546272827377 m corresponds to unwrapped x=0.6377565462728274 m. Raw FunDEM coordinates remain unwrapped. This equivalence is specific to the isolated sphere and uniform plane.

Source manifest (retained locally; not included) contains 11 official snapshots with Git blob verification and 12 current native source snapshots with SHA256. The official YAML, analytic checker, wall code and tangential model were independently fetched live at the pinned commit. Other immutable same-pin snapshots were reused and byte-verified. Independent physics JSON (retained locally; not included) contains observed numerical metrics and per-step checks. Full run archives and frame hashes are retained separately; native execution evidence determines completion, not this source inspection alone.

## Observed four-run independent audit

All four runs pass the official final 1% vx and omega_y checks. This independent script audited every numerical CSV step and all source/input mappings.

| Run | Final vx (m/s) | Final omega_y (rad/s) | Relative vx error | Relative omega error |
|---|---:|---:|---:|---:|
| cpu | 0.714236805511 | 285.761174702955 | 0.006847228% | 0.016411146% |
| gpu | 0.714236805494 | 285.761174719443 | 0.006847231% | 0.016411152% |
| cpu-dt-half | 0.714237445206 | 285.760535127898 | 0.006757671% | 0.016187295% |
| gpu-dt-half | 0.714237445203 | 285.760535130940 | 0.006757672% | 0.016187296% |

Across all four runs, maximum tangential-history return residual is 1.11008e-21 m; maximum horizontal/vertical impulse residual is 1.2719e-17 kg m/s; maximum angular impulse residual is 5.9557e-21 kg m2/s. The observed Coulomb ratio never exceeds 0.4 beyond floating-point roundoff.

CPU dt first geometric-slip zero is 0.714296411716 s. Final geometric slip is -0.000161025043156 m/s, confirming the endpoint is a finite-spring oscillatory state. Initial-to-final total mechanical energy loss is 2.52451950677e-05 J; conservation is not the acceptance target for this frictional problem.
