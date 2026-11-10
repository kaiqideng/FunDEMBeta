# Hertz oblique wall impact: LAMMPS benchmark verification

**Case:** `dem04-hertz-oblique`  
**Evidence date:** 29 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines the normal and tangential impulses and generated spin of an oblique Hertz sphere–wall collision on a common Coulomb-capped path using pinned LAMMPS inputs and analytic assertions. At the base step, the relative rebound errors in v_z, v_x and ω_y are 1.94271769117e-06, 2.2751099776e-05 and 0.00450781769376, respectively, within the 0.01 limits. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `29dd3fd7642a6fe9e8b4dd99120df4068b221250` and the input Git blob is `5b07ff13041a03aff78974ac92de06db51eafcc0`. LAMMPS was not rerun locally. Upstream qualification: `granular`.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / native representation |
|---|---|
| sphere: initial state | R=0.0025 m; position=(0, 0, 0.00255) m; velocity=(2.5069, 0, -2.9876) m/s; spin=(0, 0, 0) rad/s. |
| sphere: mass and inertia | mass=0.000176714586764 kg; isotropic inertia=4.41786466911e-10 kg m^2. |
| Normal elastic law | Hertz Fn_elastic=(4/3)*E_eff*sqrt(R_eff)*delta^1.5; native E_eff=3.9277275e+10 Pa. Stored linear kn is inactive. |
| Damping convention | Automatic damping uses the harmonic material restitution and the active normal/tangential stiffness; stored manual dashpots are inactive. |
| Tangential parameters | kt=0 N/m (explicit in linear mode), mu=0.092 |
| Rolling and twisting parameters | kr=0, mu_r=0, ktwist=0, mu_twist=0 |
| Gravity | (0, 0, 0)m/s^2 |
| Base / refined time step (s) | 1e-08/5e-09 |
| Base / refined integration steps | 4217/8434 |
| Native material encoding | E = 73644890625 Pa, nu = 0.25, encoded restitution = 0.98, density = 2700 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

Source wall spin-slip and torque use R; native uses the shared contact point at R-delta. Maximum source/native contact-speed difference is 0.00512712 m/s; torque-arm difference is 1.40351e-05 m. Official selected-time acceptance does not imply pointwise kinematic identity. Source constant kt=1.5e10 N/m and ct=0 differ from native overlap-dependent Hertz/Mindlin stiffness and automatic damping. Their total Ft agrees only under the proved common saturation. Halving dt cannot remove the R versus R-delta torque convention.


- The four runs are accepted only for their independently proved common-Coulomb-cap path. Unsaturated tangential stiffness, damping and history operators are not validated.
- Direct Hertz conversion uses E_eff=0.75kn. Native nu=0.25 is an explicit unfitted encoding choice because the source direct-kn law supplies no Poisson ratio.
- The source constant kt=1.5e10 N/m and ct=0 are replayed separately; native Hertz/Mindlin kt(delta) and automatic ct(delta) differ. All actual contact states must yield the same capped total force vector.
- Source wall torque uses R; native shared-point torque uses R-delta. The persistent spin offset is quantified and cannot be eliminated by timestep refinement.

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

Model `oblique_impact` is evaluated at original step 3317 and time 3.317e-05s. Half-step runs use twice the native check step at the same physical time. The exact reference, relative/absolute error type and unchanged limit are listed separately for each quantity. Whole histories and source golden states supplement these selected-time assertions.

For a native quantity \(q\) and reference \(q_{\rm ref}\), the error definitions are

$$
\varepsilon_{\rm rel}=\frac{|q-q_{\rm ref}|}{|q_{\rm ref}|},\qquad\varepsilon_{\rm abs}=|q-q_{\rm ref}|.
\tag{2}
$$

An error is accepted when it does not exceed its stated limit. Relative errors and limits are dimensionless; percentage columns equal \(100\varepsilon_{\rm rel}\). Absolute errors carry the units of the quantity. Zero-reference assertions use their stated absolute bounds rather than dividing by zero. The analytic tolerances in [1,2], or separately labeled supplemental diagnostic limits, determine acceptance. The YAML same-code regression epsilon is retained as source metadata and is not adopted as a cross-software acceptance tolerance.

## 3. Results and discussion

**Table 2.** Retained measurements for the two backends and two time steps. Numerical values, error types, limits and decisions are reproduced from the campaign record without refitting. `dt` denotes \(\Delta t\), and `dt/2` denotes the refined step; `um` denotes µm where used.

| Run | Quantity | Reference | Native value | Error type | Error | Limit | Result |
|---|---|---:|---:|---|---:|---:|---|
| CPU dt | Rebound vz (m/s) | 2.927848 | 2.92784231202 | relative | 1.94271769117e-06 | 0.01 | PASS |
| CPU dt | Rebound vx (m/s) | 1.962678784 | 1.9626341309 | relative | 2.2751099776e-05 | 0.01 | PASS |
| CPU dt | Rebound omega_y (rad/s) | 544.221216 | 541.767965973 | relative | 0.00450781769376 | 0.01 | PASS |
| GPU dt | Rebound vz (m/s) | 2.927848 | 2.92784231202 | relative | 1.94271769117e-06 | 0.01 | PASS |
| GPU dt | Rebound vx (m/s) | 1.962678784 | 1.9626341309 | relative | 2.27510997765e-05 | 0.01 | PASS |
| GPU dt | Rebound omega_y (rad/s) | 544.221216 | 541.767965973 | relative | 0.00450781769377 | 0.01 | PASS |
| CPU dt/2 | Rebound vz (m/s) | 2.927848 | 2.92784934954 | relative | 4.60932160063e-07 | 0.01 | PASS |
| CPU dt/2 | Rebound vx (m/s) | 1.962678784 | 1.96263407953 | relative | 2.27772711529e-05 | 0.01 | PASS |
| CPU dt/2 | Rebound omega_y (rad/s) | 544.221216 | 541.768004547 | relative | 0.00450774681542 | 0.01 | PASS |
| GPU dt/2 | Rebound vz (m/s) | 2.927848 | 2.92784934954 | relative | 4.60932159911e-07 | 0.01 | PASS |
| GPU dt/2 | Rebound vx (m/s) | 1.962678784 | 1.96263407953 | relative | 2.27772711534e-05 | 0.01 | PASS |
| GPU dt/2 | Rebound omega_y (rad/s) | 544.221216 | 541.768004547 | relative | 0.00450774681542 | 0.01 | PASS |

The spin discrepancy remains about 0.4508% when the time step is halved, while the normal-velocity error decreases. This persistent offset is consistent with the distinct torque arms: the source uses R and the native wall contact uses R − δ. The source constant k_t and zero c_t differ from native overlap-dependent stiffness and automatic damping. Independent replay proves agreement only for the total tangential force along the common saturated path. Neither the selected-time PASS nor refinement establishes equivalence of the unsaturated spring, damping or history operators.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem04-hertz-oblique. Panels (a) and (b) show velocity (z component) and angular velocity (y component) of the selected body 'sphere'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows recomputed error divided by its stated tolerance; the dashed vertical line is the acceptance boundary ε/τ = 1. Absolute errors are |actual−expected|; relative errors divide this by |expected|. Rebound vz (m/s): relative tolerance 0.01; Rebound vx (m/s): relative tolerance 0.01; Rebound omega_y (rad/s): relative tolerance 0.01. Exactly 0 zero ratios are displayed at the explicit logarithmic floor 1.15233039978e-05; every positive ratio is unmodified. The four run markers have deterministic vertical offsets within each categorical criterion row solely to distinguish coincident values. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

Conditional PASS only: replay proves equal total tangential force on a shared Coulomb cap at every contact. Unsaturated stiffness, damping and histories differ.

This report reanalyzes the retained 29 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem04-hertz-oblique-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem04-hertz-oblique-3d-si.yaml), commit `29dd3fd7642a6fe9e8b4dd99120df4068b221250`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem04_hertz_oblique --parallel
build/validation/LAMMPSComparison/lammps_dem04_hertz_oblique --mode cpu --output build/comparisons/dem04-hertz-oblique-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem04-hertz-oblique.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_metrics_comparison.csv](checkpoints/official_metrics_comparison.csv)
- [Checkpoint: source_checkpoint_comparison.csv](checkpoints/source_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
