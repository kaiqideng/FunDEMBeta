# Hooke/history oblique wall impact: LAMMPS benchmark verification

**Case:** `dem04-hookehistory-oblique`  
**Evidence date:** 28 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines frictional impulse and spin generation in oblique Hooke/history impact on a plane using pinned LAMMPS inputs and analytic assertions. At 28.57 µs, the base-step v_x, v_z and spin errors are 0.001658%, 0.000687% and 0.270217%, respectively; all satisfy the official 1% limits. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `f9720ce2d176d615495d1d6c31349b1f71be1c73` and the input Git blob is `6ab54dcfe7fcef2d4dcc16eb430a29ba8ebe8f3a`. LAMMPS was not rerun locally. The pinned upstream inputs and their original qualifications are retained.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / FunDEM configuration |
|---|---|
| Sphere diameter; density | 0.005 m; 2700 kg/m^3 |
| Solid-sphere mass; inertia | 0.000176714586764 kg; 4.41786466911e-10 kg m^2 |
| Initial center | (0, 0, 0.002550) m |
| Initial translation; spin | (2.5069, 0, -2.9876) m/s; omega = 0 |
| Wall; gravity | Fixed plane z = 0; gravity = 0 |
| Normal / tangential contact | Hooke / linear_history; kn = 2.0e7 N/m; kt = 5.7e6 N/m |
| Restitution; sliding friction | e = 0.98; mu = 0.092 |
| Explicit damping mapping | cn = 0.764595808533 N s/m; ct = cr = ctorsion = 0 |
| Rolling / torsional / global damping | Rolling and torsional stiffness/friction = 0; global damping = 0 |
| Time step; official segments | dt = 1e-8 s; 1724 + 1133 + 900 = 3757 steps |
| Refined run | dt = 5e-9 s; 7514 steps; same 37.57 us endpoint |
| Native saved output | Every integration step, including t = 0 |
| Native material encoding | E = 10000 Pa, nu = 0.33, encoded restitution = 0.98, density = 2700 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

The official gross-sliding impulse approximation for a solid sphere with zero initial spin is [2]
$$
v_z^+=e v_{z,\rm in},\qquad \Delta v_x=\mu(1+e)v_{z,\rm in}.
\tag{1a}
$$

$$
v_x^+=v_{x,\rm in}-\Delta v_x,\qquad\omega_y^+=\frac{5\Delta v_x}{2R}.
\tag{1b}
$$

Here \(v_{z,\rm in}\) is the positive incoming normal speed. These are component targets, not an assertion that the source and native contact-point kinematics coincide.

## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

1% vx, vz and omega_y at 28.57 us.

The three independent analytic targets are \(v_x^+=1.962678784\,\mathrm{m/s}\), \(v_z^+=2.927848\,\mathrm{m/s}\) and \(\omega_y^+=544.221216\,\mathrm{rad/s}\). Each must satisfy the 0.01 relative limit; passing only a subset is insufficient.

For a native quantity \(q\) and reference \(q_{\rm ref}\), the error definitions are

$$
\varepsilon_{\rm rel}=\frac{|q-q_{\rm ref}|}{|q_{\rm ref}|},\qquad\varepsilon_{\rm abs}=|q-q_{\rm ref}|.
\tag{2}
$$

An error is accepted when it does not exceed its stated limit. Relative errors and limits are dimensionless; percentage columns equal \(100\varepsilon_{\rm rel}\). Absolute errors carry the units of the quantity. Zero-reference assertions use their stated absolute bounds rather than dividing by zero. The analytic tolerances in [1,2], or separately labeled supplemental diagnostic limits, determine acceptance. The YAML same-code regression epsilon is retained as source metadata and is not adopted as a cross-software acceptance tolerance.

## 3. Results and discussion

**Table 2.** Retained measurements for the two backends and two time steps. Numerical values, error types, limits and decisions are reproduced from the campaign record without refitting. `dt` denotes \(\Delta t\), and `dt/2` denotes the refined step; `um` denotes µm where used.

| Run | dt (s) | Steps / native frames | vx (m/s) | vx error (%) | vz (m/s) | vz error (%) | omega_y (rad/s) | Spin error (%) | All 1% checks |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| CPU dt | 1e-08 | 3,757 / 3,758 | 1.962646252 | 0.001658 | 2.927827885 | 0.000687 | 542.75063903 | 0.270217 | PASS |
| GPU dt | 1e-08 | 3,757 / 3,758 | 1.962646252 | 0.001658 | 2.927827885 | 0.000687 | 542.75063903 | 0.270217 | PASS |
| CPU dt/2 | 5e-09 | 7,514 / 7,515 | 1.962651011 | 0.001415 | 2.927892057 | 0.001505 | 542.74583807 | 0.271099 | PASS |
| GPU dt/2 | 5e-09 | 7,514 / 7,515 | 1.962651011 | 0.001415 | 2.927892057 | 0.001505 | 542.74583807 | 0.271099 | PASS |

Refinement changes the spin error from 0.270217% to 0.271099%, so the spin offset is not eliminated by a smaller time step. The source wall uses nominal R in spin-slip and torque, while the native shared contact point uses R − δ. The explicit normal dashpot preserves the supplied restitution and zero tangential damping. All three component criteria pass; this must not be interpreted as pointwise contact-kinematic identity.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem04-hookehistory-oblique. Panels (a) and (b) show velocity (x component) and angular velocity (y component) of the selected body 'sphere'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows the signed residual of (b). Both residual panels use linear scales, with zero shown at its true ordinate. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

Selected rebound velocities and spin are tested; nominal-radius and shared-point wall lever arms are not identical.

This report reanalyzes the retained 28 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem04-hookehistory-oblique-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem04-hookehistory-oblique-3d-si.yaml), commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem04_hookehistory_oblique --parallel
build/validation/LAMMPSComparison/lammps_dem04_hookehistory_oblique --mode cpu --output build/comparisons/dem04-hookehistory-oblique-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem04-hookehistory-oblique.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_checkpoint_comparison.csv](checkpoints/official_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
