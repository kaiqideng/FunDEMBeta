# Independent contact-law audit: dem01-hooke-twosphere

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

## Reference and acceptance

The official LAMMPS input is pinned to commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`, input Git blob `57c2041792b5e3744e978e4a1834e46ceaaadb7b`. [Official input](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem01-hooke-twosphere-3d-si.yaml).

The model is two equal, free spheres: diameter 0.005 m, density 2700 kg/m3, centers x = -/+0.00255 m and initial x velocity +/-2 m/s, zero angular velocity. Normal Hooke stiffness is 1e6 N/m and material restitution is 0.5. All tangential force, rolling resistance and torsion parameters are zero. There is no gravity, constraint, cohesion, global damping or applied force.

Official timestep is 3e-8 s. The 2068 + 727 steps give a final physical time 8.385e-5 s. The analytic test `collision_restitution` runs at zero-based segment 1, which is the final time after both segments. It checks:

- Relative restitution error: abs(e_measured - 0.5) / 0.5 <= 0.02.
- e_measured = -(v_left,x - v_right,x)/(2*2 m/s).
- abs(m_left*v_left,x + m_right*v_right,x) <= 0.02*(m_left + m_right)*2 m/s.

The YAML epsilon 1e-10 belongs to LAMMPS same-code regression and is not a cross-code acceptance tolerance. No local LAMMPS rerun is claimed: pinned saved states and the pinned analytic checker supply the reference.

## Automatic damping mapping

LAMMPS `hooke kn en damping coeff_restitution` interprets the second normal coefficient as restitution e. Its Hooke prefactor is -2 ln(e)/sqrt(pi^2 + ln(e)^2); the inherited force routine multiplies it by sqrt(m_reduced*Fn_elastic/overlap). With Fn_elastic = kn*overlap, this is the constant dashpot:

    beta = -ln(e)/sqrt(pi^2 + ln(e)^2)
    cn = 2*beta*sqrt(m_reduced*kn)

FunDEM uses `linearStiffnessEnabled=true` and `dampingEnabled=false`. The false flag means restitution-derived damping is active. Both materials have e=0.5; the harmonic mean is therefore also 0.5. The complete sphere-pair linear stiffness group must explicitly set kt=kr=ktors=0. The corresponding automatic tangential, rolling and torsional dashpots then vanish.

Calculated constants:

| Quantity | Value |
|---|---:|
| Sphere mass | 0.00017671458676442588 kg |
| Reduced mass | 8.835729338221294e-5 kg |
| beta | 0.2154537619662468 |
| Automatic cn | 4.050468687151587 N s/m |
| Effective sphere-pair radius | 0.00125 m |

Young modulus and Poisson ratio are schema placeholders in this explicit linear stiffness branch; they do not define this official kn. Stored explicit dashpot values, if zero, are inactive while dampingEnabled=false; they must not be presented as the actual physical damping.

## Force timing and signed unloading

Both implementations retain the total signed normal force without `limit_damping` clipping. At positive overlap,

    signed Fn = kn*overlap - cn*v_normal,evaluation.

FunDEM's contact force is assembled after drift at the half-step velocity. The saved completed-step velocity includes the final half-kick, so reconstruct:

    v_evaluation = v_saved - dt*F_saved/(2*m).

Using saved full-step velocity directly would bias the inferred dashpot coefficient. The initial gap makes t=0 force exactly zero; no zero-time reconstruction is needed.

The normal output magnitude is abs(signed Fn). A negative tail during unloading is the retained dashpot convention, not a cohesive constitutive model. Friction is zero, so the tail cannot drive a friction-cap discrepancy. The native contact master can change; independent auditing canonicalizes normal and force to the left sphere before comparison.

## Continuum diagnostics

For relative overlap delta after first touch, m_reduced*delta'' + cn*delta' + kn*delta = 0. The relative incoming speed is 4 m/s. First touch is 25 microseconds. With a=cn/(2*m_reduced) and wd=sqrt(kn/m_reduced-a^2),

    delta(tau) = 4/wd * exp(-a*tau) * sin(wd*tau)
    contact duration = pi/wd.

Predicted geometric separation is 55.2407516365 microseconds, maximum overlap is 27.8916565193 micrometers at 38.0300426838 microseconds. At separation, expected sphere velocities are -/+1 m/s, giving e=0.5.

Initial total kinetic energy is 0.000706858347058 J. Expected separated kinetic energy is 0.000176714586764 J, a quarter of the initial value; energy loss is 0.000530143760293 J. Mechanical energy is intentionally dissipated. Compare K + 0.5*kn*delta^2 and integrated viscous power as diagnostics; do not require mechanical energy conservation.

The per-particle momentum impulse has magnitude 0.000530143760293 N s. Total linear momentum and total orbital plus spin angular momentum remain zero. No spin is expected.

## Why dem08 DMT was not selected

The official dem08-dmt-twosphere input is a run-0 pull-off force check at nearly zero overlap, requiring abs(F)=4*pi*gamma*R_effective=0.00628318530718 N for gamma=0.1 and R_effective=0.005 m. Current FunDEM sphere material/contact response has no surface energy, DMT pull-off term or adhesive normal-force model. Its negative viscous tail requires separating velocity, so it cannot reproduce static DMT adhesion. Replacing the DMT force with a fitted dashpot or applied load would no longer verify the official contact law.

## Evidence

All four completed native runs pass the official restitution and momentum checks. Original-step CPU/GPU restitution is 0.4996646896338126 / 0.4996646896338127; half-step CPU/GPU restitution is 0.49984130980353814 / 0.49984130980353825. Relative restitution errors are about 0.0670621% and 0.0317380%. Every-contact damping inference gives mean cn=4.050468687151586 N s/m, and signed force reconstruction residual is at most 7.11e-15 N. Total momentum and angular momentum stay exactly zero in these symmetric traces; the largest particle impulse closure residual is 2.93e-18 kg m/s.

The sampled maximum overlap changes from 27.8858964901 to 27.8823971880 micrometers at half-step, compared with the continuous 27.8916565193 micrometers. Its small relative diagnostic error increases from 0.0206514% to 0.0331975%; monotonic improvement of every metric is not claimed. Contact activation is phase dependent and the viscous force has a finite jump at first positive overlap. The official acceptance concerns restitution and momentum. Independent root-analysis crosschecks agree on 32 core values and four inferred damping means.

- Official source snapshot manifest (retained locally; not included): nine complete official snapshots verified against their pinned Git blob hashes.
- Twelve native source snapshots with SHA256 are listed in that same manifest under fundem_source_snapshots.
- Independent executable audit (retained locally; not included) verifies every output step and each stored contact, including force timing, automatic damping, impulse closure and the exact official criterion.
- Independent numeric audit (retained locally; not included) contains the four completed-run checks once native CSVs are available.
- Pre-existing artifact baseline (retained locally; not included) retains the eight previous manifest records and report/reference hashes before metadata changes.

