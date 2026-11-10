# Hertz compression and peak force: LAMMPS benchmark verification

**Case:** `dem02-hertz-peak`  
**Evidence date:** 29 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines separate overlap and contact-force assertions at the prescribed near-peak state of an elastic two-sphere Hertz impact using pinned LAMMPS inputs and analytic assertions. At 30.16 µs, the base-step overlap and force relative errors are 2.6497570337e-07 and 7.70564282514e-07, respectively, below the 0.01 limits. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `29dd3fd7642a6fe9e8b4dd99120df4068b221250` and the input Git blob is `276bff16ac211f00347db6b5e20d23a25210e4ed`. LAMMPS was not rerun locally. Upstream qualification: `unvetted`.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / native representation |
|---|---|
| left: initial state | R=0.01 m; position=(-0.0101, 0, 0) m; velocity=(10, 0, 0) m/s; spin=(0, 0, 0) rad/s. |
| right: initial state | R=0.01 m; position=(0.0101, 0, 0) m; velocity=(-10, 0, 0) m/s; spin=(0, 0, 0) rad/s. |
| left: mass and inertia | mass=0.0117286125734 kg; isotropic inertia=4.69144502936e-07 kg m^2. |
| right: mass and inertia | mass=0.0117286125734 kg; isotropic inertia=4.69144502936e-07 kg m^2. |
| Normal elastic law | Hertz Fn_elastic=(4/3)*E_eff*sqrt(R_eff)*delta^1.5; native E_eff=2.4999975e+10 Pa. Stored linear kn is inactive. |
| Damping convention | Automatic damping uses the harmonic material restitution and the active normal/tangential stiffness; stored manual dashpots are inactive. |
| Tangential parameters | kt=0 N/m (explicit in linear mode), mu=0 |
| Rolling and twisting parameters | kr=0, mu_r=0, ktwist=0, mu_twist=0 |
| Gravity | (0, 0, 0)m/s^2 |
| Base / refined time step (s) | 2e-08/1e-08 |
| Base / refined integration steps | 2508/5016 |
| Native Hertz force coefficient | \(K_{\rm native}\) = 2357020246.932554 N/m^(3/2), derived from the retained contact coefficient |
| Official analytic force coefficient | \(K_{\rm check}\) = 2357023000 N/m^(3/2), the independently rounded YAML `kfac` |
| Native material encoding | E = 46874953125 Pa, nu = 0.25, encoded restitution = 1, density = 2800 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

- Direct Hertz kn mapped algebraically using E*=3kn/4 and equal-material E=2(1-nu^2)E*. Only normal channel active.

The Hertz relations below use a generic force coefficient \(K\), reduced mass \(m^*\) and approach speed \(V\). Native force and stored energy use \(K_{\rm native}=(4/3)E^*\sqrt{R^*}=k_n\sqrt{R^*}\). The official `hertz_peak` checker substitutes the separately rounded YAML `kfac`, denoted \(K_{\rm check}\), into Eq. (1b) to generate the reference overlap and force in Table 2 [1,2]. Both coefficients are retained independently; neither is fitted or silently replaced by the other.
$$
F_n=K\delta^{3/2},\qquad U(\delta)=\frac{2}{5}K\delta^{5/2}=\frac{2}{5}F_n\delta.
\tag{1a}
$$

$$
\delta_{\max}=\left(\frac{5m^*V^2}{4K}\right)^{2/5},\qquad F_{\max}=K\delta_{\max}^{3/2}.
\tag{1b}
$$


## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

Model `hertz_peak` is evaluated at original step 1508 and time 3.016e-05s. Half-step runs use twice the native check step at the same physical time. The exact reference, relative/absolute error type and unchanged limit are listed separately for each quantity. Whole histories and source golden states supplement these selected-time assertions.

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
| CPU dt | Overlap at official check time (m) | 0.000274112683666 | 0.000274112756299 | relative | 2.6497570337e-07 | 0.01 | PASS |
| CPU dt | Force at official check time (N) | 10696.8897029 | 10696.8814603 | relative | 7.70564282514e-07 | 0.01 | PASS |
| GPU dt | Overlap at official check time (m) | 0.000274112683666 | 0.000274112756299 | relative | 2.6497570337e-07 | 0.01 | PASS |
| GPU dt | Force at official check time (N) | 10696.8897029 | 10696.8814603 | relative | 7.70564282514e-07 | 0.01 | PASS |
| CPU dt/2 | Overlap at official check time (m) | 0.000274112683666 | 0.000274112732997 | relative | 1.79965200119e-07 | 0.01 | PASS |
| CPU dt/2 | Force at official check time (N) | 10696.8897029 | 10696.8800963 | relative | 8.98079902594e-07 | 0.01 | PASS |
| GPU dt/2 | Overlap at official check time (m) | 0.000274112683666 | 0.000274112732997 | relative | 1.79965200119e-07 | 0.01 | PASS |
| GPU dt/2 | Force at official check time (N) | 10696.8897029 | 10696.8800963 | relative | 8.98079902594e-07 | 0.01 | PASS |

The refined overlap error is 1.79965200119e-07, while the force error is 8.98079902594e-07. Refinement therefore improves one sampled quantity and slightly worsens the other. These are values at the official sample time, not independently located global maxima. CPU/GPU values coincide at the reported precision. The result tests only the elastic normal channel; all damping and frictional contributions are inactive.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem02-hertz-peak. Panels (a) and (b) show position (x component) and velocity (x component) of the selected body 'left'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows recomputed error divided by its stated tolerance; the dashed vertical line is the acceptance boundary ε/τ = 1. Absolute errors are |actual−expected|; relative errors divide this by |expected|. Overlap at official check time (m): relative tolerance 0.01; Force at official check time (N): relative tolerance 0.01. Exactly 0 zero ratios are displayed at the explicit logarithmic floor 4.49913000297e-06; every positive ratio is unmodified. The four run markers have deterministic vertical offsets within each categorical criterion row solely to distinguish coincident values. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

Normal elastic impact. Source coefficient conversion is algebraic, not fitted; tangential response is inactive.

This report reanalyzes the retained 29 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem02-hertz-peak-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem02-hertz-peak-3d-si.yaml), commit `29dd3fd7642a6fe9e8b4dd99120df4068b221250`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem02_hertz_peak --parallel
build/validation/LAMMPSComparison/lammps_dem02_hertz_peak --mode cpu --output build/comparisons/dem02-hertz-peak-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem02-hertz-peak.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_metrics_comparison.csv](checkpoints/official_metrics_comparison.csv)
- [Checkpoint: source_checkpoint_comparison.csv](checkpoints/source_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
