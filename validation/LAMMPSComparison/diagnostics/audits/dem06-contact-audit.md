# DEM06 Hooke-history spinning-sphere impact: independent contact audit

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 2026-10-28. Scope: exact official input, source mapping, and later independent force/impulse reconstruction. No solver source is modified.

## Official source and acceptance

The pinned [LAMMPS input](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem06-hookehistory-spin-3d-si.yaml) has Git blob d1cc7b040decc016ec05e95b5e0dc2e283406434. It records version 4 Jul 2026, generated 21 Jul 2026. The two particles are finite, solid spheres and both are integrated by fix nve/sphere. The large sphere is a very heavy moving collision partner, not a fixed wall.

| Quantity | Official value |
|---|---:|
| Small diameter, density | 0.2 m, 2700 kg/m3 |
| Large diameter, density | 1.0 m, 2.7e6 kg/m3 |
| Initial centers | small (0,0,0.60001) m; large (0,0,0) m |
| Small initial velocity | (0,0,-0.2) m/s |
| Small initial spin | (0,22.86,0) rad/s |
| Large initial velocity and spin | zero |
| Normal Hooke stiffness | 1.2e8 N/m |
| Tangential linear-history stiffness | 3.4e7 N/m |
| Normal restitution input | 0.5 |
| Sliding coefficient | 0.40 |
| Tangential damping multiplier | 0.0 |
| Rolling, twisting, gravity, cohesion | none |
| Time step and run segments | 1 us; 600 + 700 + 500 steps |

The official analytic model is spin_impact at zero-based segment 1: cumulative step 1300, time 1.3 ms. All three small-sphere outputs must independently satisfy relative error <= 0.07 against vx = 0.12 m/s, vz = 0.1 m/s, omega_y = 19.86 rad/s. The half-step run must check step 2600 at the same physical time. The complete trajectory runs to 1.8 ms (1800 or 3600 steps).

The [analytic source](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp#L197-L220) uses the heavy-partner gross-sliding approximation:

- Delta vx = mu (1+e) vin = 0.12 m/s.
- vz_out = e vin = 0.1 m/s.
- omega_y_out = omega0 - 2.5 Delta vx / Rsmall = 19.86 rad/s.

The YAML epsilon 1e-10 is a same-code stored-trajectory regression tolerance, not the cross-code analytic tolerance. Stored LAMMPS output at 1.3 ms is vx = 0.12706938882175314 m/s, vz = 0.099865657293162532 m/s, omega_y = 19.683950332657247 rad/s. The stored vx is already 5.89% above the ideal approximation. A successful cross-code check should retain that distinction.

## Explicit material mapping

Masses msmall = 11.309733552923257 kg, mlarge = 1413716.6941154068 kg. Solid-sphere inertias are Ismall = 0.04523893421169304 kg m2 and Ilarge = 141371.6694115407 kg m2. The [LAMMPS pair code](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/pair_granular.cpp#L241-L252) uses the finite reduced mass

mred = msmall mlarge / (msmall+mlarge) = 11.309643075778652 kg.

For the Hooke model and damping coeff_restitution:

beta = -ln(e) / sqrt(pi^2 + ln(e)^2) = 0.2154537619662468,

cn = 2 beta sqrt(mred kn) = 15874.466906486981 N s/m.

There is no Hertz sqrt(5/6) factor. Since the official tangential multiplier is zero, ct is exactly zero. The FunDEM configuration must use all-linear stiffness override and explicit damping override: kn = 1.2e8 N/m, kt = 3.4e7 N/m, cn above, ct = 0, rolling/torsional stiffness/dashpots/friction = 0, and global damping = 0. E and Poisson ratio required by the schema are inactive placeholders in this mode, not extra official material parameters.

FunDEM automatic damping would also derive nonzero tangential damping. Therefore this exact official mapping deliberately selects explicit coefficients. This case verifies the official Hooke friction/spin behavior rather than retesting the automatic-damping branch.

## Signed normal force and rotating contact geometry

Let n point from large center to small center, delta = Rsmall + Rlarge - |xsmall-xlarge|. For positive geometric overlap,

Fn = kn delta - cn dot(vsmall_eval-vlarge_eval,n).

The source defaults limit_damping = 0 and the input does not enable it. Both solvers retain the negative signed normal unloading tail and use mu |Fn| for the tangential Coulomb limit. Do not replace signed Fn by its stored magnitude in normal impulse integration, and do not assume integral |Fn| dt equals integral Fn dt. With e = 0.5 the attractive tail is material to the comparison with the ideal gross-sliding formula.

FunDEM uses one shared sphere/sphere point halfway between the facing undeformed surfaces:

point = xlarge + (Rlarge-delta/2) n = xsmall - (Rsmall-delta/2) n.

Hence its force-generated torques are tau_small = (point-xsmall) cross Fsmall and tau_large = (point-xlarge) cross (-Fsmall). The extra contact.torque field is zero because rolling and torsion are disabled.

The [LAMMPS pair torque](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/src/GRANULAR/granular_model.cpp#L498-L508) uses the same Ri-delta/2 arms. The wall-case R versus R-delta torque difference does not apply here. There are nevertheless two source-level algorithm differences:

1. LAMMPS contact slip uses nominal Ri omega_i + Rj omega_j in its rotational term. FunDEM evaluates both rigid-body velocities at the shared point, shortening each rotational arm by delta/2.
2. LAMMPS projects and rescales its old tangential history when rotating the contact normal. FunDEM projects without norm rescaling before adding the current slip increment. Unlike a flat-wall test, this contact normal rotates slightly as the small sphere acquires x velocity.

Both schemes use a Coulomb radial return with ct = 0. The native law must be checked against its actual common-point half-step slip and projected spring history. These differences should be disclosed without fitting parameters or assuming an unmeasured attribution.

## Timing and independent checks

FunDEM uses velocity Verlet. Completed native frames include the final velocity/angular-velocity half kick, but force and torque were evaluated immediately before it. For each particle j:

v_eval_j = v_saved_j - dt F_saved_j/(2 mj),

omega_eval_j = omega_saved_j - dt tau_saved_j/(2 Ij).

The saved force must use these reconstructed velocities. The per-step finite-difference checks are Delta v = dt (F_previous+F_current)/(2m) and Delta omega = dt (tau_previous+tau_current)/(2I). Trapezoidal integrals of each particle's force and torque should match its linear and intrinsic angular momentum changes.

Additional diagnostics:

- Equal and opposite pair forces, Coulomb bound |Ft| <= mu |Fn|, shared-point torque closure, and zero extra resistance torque.
- Integral signed Fn and integral |Fn|; tensile tail contribution; cap fraction and slip sign.
- Total linear momentum P = sum m v, with large-sphere recoil explicitly retained.
- Total angular momentum about a fixed origin L = sum (x cross mv + I omega). Shared-point internal forces give zero total internal moment; finite integration may retain a small discretization residual that should be quantified rather than asserted exactly zero.
- Actual contact point, normal/overlap geometry, explicit cn and ct, both masses/inertias, and inactive transverse components.
- CPU/GPU trajectory agreement and dt/2 differences at common physical times, without assuming refinement must remove the ideal-approximation discrepancy.

## Source retention

reference/lammps/dem06-source contains nine complete live-fetched official files at the pinned commit. source-manifest.json records official and independently calculated Git blob hashes, SHA-256 hashes, and sizes. All nine official Git blobs match. Relevant native sources are copied without editing into reference/fundem/dem06-source; the manifest records their SHA-256 hashes. These document inspected source, while native-audit/build evidence documents actual linked runtime provenance. This does not establish byte identity with an installed GUI executable.

The official f/f/f box is a neighbor/deletion domain; it does not add physical walls. Both particles remain inside it by center coordinates. The authored native project must have no external walls, prescribed motion, gravity or periodic force modules.


## Independent completed-run results

All four runs were independently read by tools/audit_dem06_physics.py. The machine-readable evidence is tools/dem06-physics-audit.json. The official small-sphere checks all pass; every one of vx, vz, omega_y is within the official 7% tolerance at 1.3 ms.

| Run | vx (m/s) | vz (m/s) | omega_y (rad/s) |
|---|---:|---:|---:|
| CPU dt | 0.12706941711922226 | 0.09986565350999287 | 19.683949625161674 |
| GPU dt | 0.12706941711921863 | 0.09986565350999180 | 19.683949625161823 |
| CPU dt/2 | 0.12705825036481969 | 0.10004190299640328 | 19.684230091611113 |
| GPU dt/2 | 0.12705825036482080 | 0.10004190299640738 | 19.684230091611120 |

The base-CPU analytic relative errors are 5.89118093% in vx, 0.13434649% in vz, and 0.88645707% in spin. The time-step refinement changes final vx by -1.11667544e-5 m/s, vz by +1.76249486e-4 m/s, and spin by +2.80466449e-4 rad/s. Refinement reduces each analytic discrepancy here, but it does not remove the approximate gross-sliding vx gap.

Base CPU has 988 contact samples, 970 Coulomb-capped and 18 uncapped. Of these samples, 137 retain a negative normal force with minimum -1584.198852277977 N. Signed scalar normal impulse is 3.391438289861633 N s, absolute normal impulse is 3.6195026567997686 N s, and tensile impulse magnitude is 0.1140321834690677 N s. Thus the absolute impulse exceeds the signed impulse by about 6.725%. The tangential impulse x component is 1.4370547329028547 N s, slightly below mu times the absolute normal impulse because of the initial spring buildup. A small normal x impulse of 6.651744080819267e-5 N s raises the total particle x impulse to 1.4371212503436634 N s. This measured decomposition explains the scale of the ideal-approximation gap without parameter changes.

The reconstructed slip stays negative in x throughout contact; its magnitude remains at least 1.84153 m/s. There is no sticking or slip reversal. The native normal force and the projected/capped tangential history laws reconstruct to floating-point precision. Per-step particle forces are exactly equal and opposite; extra rolling/torsional resistance torque remains zero. Each force/torque impulse reproduces its own particle momentum change.

Across four runs the total linear momentum drift is below 1.89e-14 kg m/s, and total angular momentum drift about a fixed origin is below 1.91e-14 kg m2/s. The finite large sphere recoils and spins as required; it has not been replaced by a fixed boundary. The approximately 125000:1 mass ratio supports the official heavy-partner idealization, while full pair momentum is still used in the native verification.

Results preserve the two source-level slip/history conventions described above. Passing the official 7% analytic criterion is not a claim of byte-exact agreement with all stored LAMMPS trajectory values, nor a test of automatic damping, nonzero tangential dashpots, sticking, or energy conservation for an inelastic frictional impact.
