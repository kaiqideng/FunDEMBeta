# Curved-wall rebound: LAMMPS benchmark verification

**Case:** `dem11-regionwall`  
**Evidence date:** 29 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines restitution against a discretized inner cylindrical wall with an analytically preserved Hertz coefficient using pinned LAMMPS inputs and analytic assertions. The base-step rebound velocity is -0.999807429792 m/s, giving relative error 0.000192570208063 against -1 m/s and the official 0.01 limit. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `29dd3fd7642a6fe9e8b4dd99120df4068b221250` and the input Git blob is `e90c87ee07bed0d26ab25bb04320a43f73ef1b4c`. LAMMPS was not rerun locally. Upstream qualification: `unvetted`.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / native representation |
|---|---|
| sphere: initial state | R=0.0025 m; position=(0.0072, 0, 0) m; velocity=(2, 0, 0) m/s; spin=(0, 0, 0) rad/s. |
| sphere: mass and inertia | mass=0.000458148928649 kg; isotropic inertia=1.14537232162e-09 kg m^2. |
| Normal elastic law | Hertz Fn_elastic=(4/3)*E_eff*sqrt(R_eff)*delta^1.5; native E_eff=5.70157227e+10 Pa. Stored linear kn is inactive. |
| Damping convention | Automatic damping uses the harmonic material restitution and the active normal/tangential stiffness; stored manual dashpots are inactive. |
| Tangential parameters | kt=0 N/m (explicit in linear mode), mu=0 |
| Rolling and twisting parameters | kr=0, mu_r=0, ktwist=0, mu_twist=0 |
| Gravity | (0, 0, 0)m/s^2 |
| Base / refined time step (s) | 1e-08/5e-09 |
| Base / refined integration steps | 20500/41000 |
| Native material encoding | E = 106904480061 Pa, nu = 0.25, encoded restitution = 0.5, density = 7000 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

- The physical cylinder radius is 0.01 m with endcaps at z=+/-0.05 m. Sphere pose, velocity, mass, restitution and duration are retained; periodic Z is inactive at z=vz=0.
- The analytic conversion E_eff=(3/4)*kn*sqrt(8/7) maps the source inner-cylinder effective radius Rp*8/7 to native Rp, preserving elastic and automatic damping coefficients without fitting.
- Sliding friction is zero. The dormant native tangential stiffness differs from source kt=0 but contributes no force in this invariant normal motion.
- The prepared SDF grid has a symmetric transverse phase. It is the authoritative saved geometry; regenerating the built-in geometry would change that phase.
- Fine/coarse signed-distance interpolation errors are measured explicitly. This is a discretized-wall verification, not a claim of exact analytic cylinder geometry.
- A metadata-only URL typo was corrected to unittest/granular/tests. Original configured input and result bytes remain preserved with the correction record.

For incoming motion along positive \(x\), the official wall-restitution target is [2]
$$
v_x^+=-e v_{x,\rm in}.
\tag{1}
$$


## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

Model `wall_restitution` is evaluated at original step 20500 and time 0.000205s. Half-step runs use twice the native check step at the same physical time. The exact reference, relative/absolute error type and unchanged limit are listed separately for each quantity. Whole histories and source golden states supplement these selected-time assertions.

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
| CPU dt | Rebound vx (m/s) | -1 | -0.999807429792 | relative | 0.000192570208063 | 0.01 | PASS |
| GPU dt | Rebound vx (m/s) | -1 | -0.999807429792 | relative | 0.000192570208063 | 0.01 | PASS |
| CPU dt/2 | Rebound vx (m/s) | -1 | -0.999881839183 | relative | 0.000118160816632 | 0.01 | PASS |
| GPU dt/2 | Rebound vx (m/s) | -1 | -0.999881839183 | relative | 0.000118160816632 | 0.01 | PASS |

The refined rebound error is 0.000118160816632, and CPU/GPU values agree at the displayed precision. Independent geometry diagnostics report maximum actual-contact axis SDF errors of 1.0417437 µm on the 0.00025 m grid and 4.1679848 µm on the 0.0005 m grid. These are a separate spatial-refinement comparison, not a time-convergence result. The effective-modulus conversion preserves the source cylindrical Hertz coefficient algebraically. The retained symmetric SDF grid phase is part of the benchmark geometry; regenerating a nominal cylinder can change it.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem11-regionwall. Panels (a) and (b) show position (x component) and velocity (x component) of the selected body 'sphere'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows recomputed error divided by its stated tolerance; the dashed vertical line is the acceptance boundary ε/τ = 1. Absolute errors are |actual−expected|; relative errors divide this by |expected|. Rebound vx (m/s): relative tolerance 0.01. Exactly 0 zero ratios are displayed at the explicit logarithmic floor 0.0029540204158; every positive ratio is unmodified. The four run markers have deterministic vertical offsets within each categorical criterion row solely to distinguish coincident values. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

Conditional on the saved symmetric prepared SDF grid and algebraic curvature mapping. Rebuilt wall geometry is not byte-equivalent; geometric interpolation error is measured separately.

This report reanalyzes the retained 29 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem11-regionwall-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem11-regionwall-3d-si.yaml), commit `29dd3fd7642a6fe9e8b4dd99120df4068b221250`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem11_regionwall --parallel
build/validation/LAMMPSComparison/lammps_dem11_regionwall --mode cpu --output build/comparisons/dem11-regionwall-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem11-regionwall.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_metrics_comparison.csv](checkpoints/official_metrics_comparison.csv)
- [Checkpoint: source_checkpoint_comparison.csv](checkpoints/source_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
