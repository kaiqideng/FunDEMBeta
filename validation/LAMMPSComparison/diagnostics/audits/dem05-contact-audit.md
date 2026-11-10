# dem05 contact-model audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 28 October 2026. Source and physics audit only; actual native trajectories determine measured pass/fail.

## Official setup and acceptance

Reference: [LAMMPS dem05-mindlin-slip-3d-si.yaml](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem05-mindlin-slip-3d-si.yaml), commit f9720ce2d176d615495d1d6c31349b1f71be1c73, version label 4 Jul 2026, regression generated 21 July 2026.

One homogeneous solid sphere: diameter 0.005 m, density 2700 kg/m3, center z = 0.0024393 m, vx = 1 m/s, zero initial spin. Gravity is exactly 1 m/s2 downward. The stationary lower plane is z = 0. The upper plane at z = 1 never touches the sphere. Contact: hertz/material E = 10000 Pa, e = 1, nu = 0.33; tangential mindlin NULL 1.0 0.4; damping coeff_restitution. No rolling/twisting resistance. Timestep 1e-5 s. Cumulative checkpoints: 25000, 50000, 75000 steps = 0.25, 0.50, 0.75 s.

The YAML epsilon = 1e-10 is internal LAMMPS regression precision. The separate analytic slip_cessation tolerance is 0.01, applied only to vx and omega_y at the final segment. Do not transfer epsilon to cross-code verification.

## Elastic mapping

The LAMMPS wall is infinite in mass, but its hertz/material mapping does not take an infinitely stiff elastic half-space limit. Pinned source formulas are:

- E* = E / [2 (1 - nu^2)].
- G* = E / [4 (2 - nu) (1 + nu)].
- Fn = (4/3) E* sqrt(R delta) delta.
- Kt = 8 G* sqrt(R delta).

Assign BOTH the FunDEM sphere and wall E = 10000 Pa, nu = 0.33. Their series-compliance combination reproduces both E* and G* exactly. Do not double E or make the wall infinitely stiff. The wall is kinematically fixed.

| Quantity | Value |
|---|---:|
| Sphere mass m | 1.767145867644259e-4 kg |
| Radius R | 0.0025 m |
| Inertia I = (2/5) m R^2 | 4.417864669110648e-10 kg m2 |
| Effective E* | 5611.04253170239 Pa |
| Effective G* | 1125.5684120480844 Pa |
| Static overlap from Fn = mg | 6.065687964282859e-5 m |
| Static center height h = R - delta | 0.0024393431203571713 m |
| h/R | 0.9757372481428686 |

Evidence:

- [gran_sub_mod.cpp lines 87-122](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod.cpp#L87): pair/wall E and G mixing.
- [gran_sub_mod_normal.cpp lines 209-220](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_normal.cpp#L209): wall Hertz uses wall-specific mixing.
- [gran_sub_mod_tangential.cpp lines 299-317](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_tangential.cpp#L299): NULL Mindlin stiffness uses wall-specific G*.
- Local FunDEMBeta/src/execution/contactMaterialFunctions.h lines 107-117: series-compliance E/G and Hertz/Mindlin stiffness.

## Damping and friction

LAMMPS coeff_restitution uses log(e), so e = 1 gives zero normal damping. The tangential multiplier xt = 1 multiplies the zero normal damping prefactor, giving zero tangential damping. The numeral 1.0 after mindlin NULL is not a dimensional dashpot.

Use both FunDEM materials' e = 1 with dampingEnabled = false (automatic restitution-derived damping, no manual override). beta = -log(e)/sqrt(pi^2 + log(e)^2) is zero, making the automatic dashpots zero. Set sliding friction 0.4, rolling friction 0, torsional friction 0, global damping 0. Frictional energy loss still occurs.

With zero dashpots, both codes use tangential displacement history, current Mindlin stiffness, and Coulomb radial return |Ft| <= mu |Fn|. This plain mindlin variant has no explicit unloading-history rescale. Damping-limiter differences are inactive here. Normal and tangential elastic oscillations can persist, so the final instant need not be an exactly steady force-free rolling state.

Evidence:

- [gran_sub_mod_damping.cpp lines 163-183](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_damping.cpp#L163).
- [gran_sub_mod_tangential.cpp lines 339-435](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/gran_sub_mod_tangential.cpp#L339).
- FunDEMBeta/src/execution/contactFunctions.h lines 89-102, 163-195; contactMaterialFunctions.h lines 120-141.

## Contact lever arm: a genuine model-convention difference

The elastic mapping above is exact. A separate kinematic convention can still cause the official 1% final check to fail. No material/friction parameter fitting is appropriate.

LAMMPS wall convention: slip = vx - R omega_y, torque tau_y = -R Fx. Its rotational contact velocity explicitly uses radi = R. Pair-contact torque has a penetration adjustment, but the wall branch explicitly scales torque by radi.

Primary evidence: [granular_model.cpp lines 471-510](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L471); [fix_wall_gran.cpp lines 577-617](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/fix_wall_gran.cpp#L577).

FunDEM sphere/level-set convention: contact point = position - (R - overlap) normal. On the lower plane this point is on z = 0, so a = R - overlap = z_center. Both rigid-body contact velocity and moment of force use this point: slip = vx - a omega_y, tau_y = -a Fx. Hertz stiffness still uses physical R. Here delta/R is approximately 2.4263%, material relative to a 1% criterion.

Local source evidence:

- FunDEMBeta/src/execution/contactDetection.h lines 181-185: sphereLevelSetContactPoint.
- FunDEMBeta/src/interaction/contactSearch.cpp lines 281-287: contact point and physical effective radius.
- FunDEMBeta/src/execution/contactFunctions.h lines 70-85: contact-point velocity.
- FunDEMBeta/src/interaction/contact.cpp lines 50-62: feeds that velocity into tangential force.
- FunDEMBeta/src/execution/cpu/particleFunctions.h lines 113-115, 147-164: cross(point - position, force) enters particle torque.

## Analytical targets

Official solid-sphere ideal:

- Gross sliding: vx = u0 - mu g t; omega_y = [5 mu g / (2 R)] t.
- Slip cessation ts = 2 u0 / (7 mu g) = 0.7142857142857143 s.
- Rolling velocity = 5 u0 / 7 = 0.7142857142857143 m/s.
- Rolling spin = velocity/R = 285.7142857142857 rad/s.
- At t = 0.75 s, BOTH relative errors in vx and omega_y must be <= 1%.

[Official analytic implementation](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp#L82) asserts only those two final components. Slip residual and cessation time are supplementary diagnostics. The three YAML regression checkpoint vectors are additional comparison data.

For FunDEM's convention, a constant-height solid-sphere approximation instead gives:

v_roll = u0 / [1 + (2/5) R^2/a^2],
omega_roll = v_roll/a,
ts = u0 / {mu g [1 + (5/2) a^2/R^2]}.

At a = static Hertz height:

| Source-based expectation, not measured | Value | Difference from official ideal |
|---|---:|---:|
| v_roll | 0.7041558362936668 m/s | -1.4181829189% |
| omega_roll | 288.66617017394566 rad/s | +1.0331595609% |
| ts | 0.7396104092658325 s | about +3.55% |

These expectations MUST NOT replace official acceptance targets. Slight vertical motion, finite tangential stiffness, and a short post-cessation window affect native final values.

Recommended direct trajectory checks:

1. Gross sliding: |Fx|/Fn = 0.4. Measure tau_y/(-Fx), compare with z_center and R.
2. Plot both vx - R omega_y (official convention) and vx - z_center omega_y (FunDEM convention).
3. At 0.25/0.50 s, vx should be near 0.9/0.8 m/s. Smaller FunDEM torque arm predicts omega_y near 97.5737/195.1474 rad/s, versus official 100/200.
4. At constant height, vx + [I/(m a)] omega_y = u0. For LAMMPS a = R this is vx + 0.001 omega_y = 1. Variable height makes the corrected expression approximate.
5. Compare CPU/GPU and dt/2 at matching timestamps. A persistent offset after dt reduction supports a convention diagnosis, without proving all other solver behavior.

## Periodic positions and an unwrapped equivalent

LAMMPS has periodic x in [-0.05, 0.5) and y in [-0.05, 0.05). One sphere cannot contact its own image since periodic lengths greatly exceed its diameter. The plane is uniform. Wrapping therefore changes only position representation. A FunDEM unwrapped trajectory on a wide uniform stationary plane is physically equivalent for this test.

For position regression: x_wrapped = -0.05 + positive_mod(x_unwrapped + 0.05, 0.55). The third official checkpoint x_wrapped = 0.08774648294987510 m corresponds to x_unwrapped = 0.6377464829498751 m. Velocities/spins compare directly. Record this representation mapping and preserve raw unwrapped coordinates.

This equivalence does not generalize to multiple spheres or short/nonuniform walls. The finite FunDEM plane must cover the full unwrapped path. A plane has linear signed distance, so interior trilinear level-set interpolation can reproduce its plane exactly if the sphere stays in the valid grid.

## Reproducible source evidence

Six pinned official sources are saved under tools/dem05-contact-source/. manifest.json records original paths, commit, Git blob hashes, URLs, and snapshot paths. Current FunDEM sources were inspected read-only. Native contact forces/torques should independently confirm the convention of the executed library. Current source alone cannot establish byte identity with a prebuilt library or installed GUI executable.

No solver code, benchmark input, friction value, or material parameter was changed by this audit.
