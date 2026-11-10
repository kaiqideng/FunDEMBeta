# Transient sliding and spin: LAMMPS benchmark verification

**Case:** `dem05-sliptime`  
**Evidence date:** 29 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines the translational and rotational transient while a sphere is still undergoing Coulomb sliding using pinned LAMMPS inputs and analytic assertions. At 0.4 s, the base-step sliding velocity is 1.24000196289 m/s and spin is 159.996905954 rad/s, within the official 0.01 relative limits. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `29dd3fd7642a6fe9e8b4dd99120df4068b221250` and the input Git blob is `20f0a60a7ead5d669a420f48008aae1fd26fe928`. LAMMPS was not rerun locally. Upstream qualification: `unvetted`.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / native representation |
|---|---|
| sphere: initial state | R=0.0025 m; position=(0, 0, 0.002499982) m; velocity=(1.4, 0, 0) m/s; spin=(0, 0, 0) rad/s. |
| sphere: mass and inertia | mass=0.000176714586764 kg; isotropic inertia=4.41786466911e-10 kg m^2. |
| Normal elastic law | Linear Hooke Fn_elastic=kn*delta, kn=10000 N/m. |
| Damping convention | Automatic damping uses the harmonic material restitution and the active normal/tangential stiffness; stored manual dashpots are inactive. |
| Tangential parameters | kt=2800 N/m (explicit in linear mode), mu=0.4 |
| Rolling and twisting parameters | kr=0, mu_r=0, ktwist=0, mu_twist=0 |
| Gravity | (0, 0, -1)m/s^2 |
| Base / refined time step (s) | 1e-05/5e-06 |
| Base / refined integration steps | 100000/200000 |
| Native material encoding | E = 10000 Pa, nu = 0.25, encoded restitution = 1, density = 2700 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

Source wall spin-slip and torque use R; native uses the shared contact point at R-delta. Maximum source/native contact-speed difference is 7.198e-06 m/s; torque-arm difference is 1.8e-08 m. Official selected-time acceptance does not imply pointwise kinematic identity.


- Source periodic x/y do not affect contact dynamics. Map saved x coordinates modulo original periodic box length0.55m for checkpoint comparison. Source ceiling atz=1m remains inactive but is retained.

Before the ideal slip-cessation time, the official Coulomb transient is [2]
$$
v_x(t)=u_0-\mu gt,\qquad\omega_y(t)=\frac{5\mu gt}{2R},\qquad t<t_s=\frac{2u_0}{7\mu g}.
\tag{1}
$$


## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

Model `slip_transient` is evaluated at original step 40000 and time 0.4s. Half-step runs use twice the native check step at the same physical time. The exact reference, relative/absolute error type and unchanged limit are listed separately for each quantity. Whole histories and source golden states supplement these selected-time assertions.

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
| CPU dt | Sliding vx (m/s) | 1.24 | 1.24000196289 | relative | 1.58297653426e-06 | 0.01 | PASS |
| CPU dt | Sliding omega_y (rad/s) | 160 | 159.996905954 | relative | 1.93377878995e-05 | 0.01 | PASS |
| GPU dt | Sliding vx (m/s) | 1.24 | 1.24000196289 | relative | 1.58297716781e-06 | 0.01 | PASS |
| GPU dt | Sliding omega_y (rad/s) | 160 | 159.996905953 | relative | 1.93377928463e-05 | 0.01 | PASS |
| CPU dt/2 | Sliding vx (m/s) | 1.24 | 1.24000145496 | relative | 1.1733511648e-06 | 0.01 | PASS |
| CPU dt/2 | Sliding omega_y (rad/s) | 160 | 159.997413889 | relative | 1.61631921143e-05 | 0.01 | PASS |
| GPU dt/2 | Sliding vx (m/s) | 1.24 | 1.24000145496 | relative | 1.17335294151e-06 | 0.01 | PASS |
| GPU dt/2 | Sliding omega_y (rad/s) | 160 | 159.997413887 | relative | 1.6163205679e-05 | 0.01 | PASS |

The reference values 1.24 m/s and 160 rad/s are evaluated before the ideal slip-cessation time, so the test concerns the sliding transient rather than final pure rolling. Refinement reduces the relative errors from 1.58297653426e-06 to 1.1733511648e-06 in velocity and from 1.93377878995e-05 to 1.61631921143e-05 in spin. Periodic x wrapping is applied only when comparing saved coordinates; it does not alter isolated-particle motion on the uniform plane.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem05-sliptime. Panels (a) and (b) show velocity (x component) and angular velocity (y component) of the selected body 'sphere'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows recomputed error divided by its stated tolerance; the dashed vertical line is the acceptance boundary ε/τ = 1. Absolute errors are |actual−expected|; relative errors divide this by |expected|. Sliding vx (m/s): relative tolerance 0.01; Sliding omega_y (rad/s): relative tolerance 0.01. Exactly 0 zero ratios are displayed at the explicit logarithmic floor 2.93337791199e-05; every positive ratio is unmodified. The four run markers have deterministic vertical offsets within each categorical criterion row solely to distinguish coincident values. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

Selected-time response before rolling is tested. A single isolated sphere permits an unwrapped equivalent; this does not generalize to many-body periodic flow.

This report reanalyzes the retained 29 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem05-sliptime-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem05-sliptime-3d-si.yaml), commit `29dd3fd7642a6fe9e8b4dd99120df4068b221250`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem05_sliptime --parallel
build/validation/LAMMPSComparison/lammps_dem05_sliptime --mode cpu --output build/comparisons/dem05-sliptime-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem05-sliptime.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_metrics_comparison.csv](checkpoints/official_metrics_comparison.csv)
- [Checkpoint: source_checkpoint_comparison.csv](checkpoints/source_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
