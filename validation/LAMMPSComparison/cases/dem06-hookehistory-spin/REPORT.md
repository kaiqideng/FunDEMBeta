# Spin-induced oblique collision: LAMMPS benchmark verification

**Case:** `dem06-hookehistory-spin`  
**Evidence date:** 28 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines spin-induced friction impulse and translation in a free two-sphere impact with strongly unequal masses using pinned LAMMPS inputs and analytic assertions. At 1.3 ms, base-step errors in v_x, v_z and ω_y are 5.89118%, 0.13435% and 0.88646%, all within the official 7% limits. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `f9720ce2d176d615495d1d6c31349b1f71be1c73` and the input Git blob is `d1cc7b040decc016ec05e95b5e0dc2e283406434`. LAMMPS was not rerun locally. The pinned upstream inputs and their original qualifications are retained.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / FunDEM configuration |
|---|---|
| Small sphere | Radius 0.1 m; density 2700 kg/m^3 |
| Large sphere | Radius 0.5 m; density 2.7e6 kg/m^3; mass ratio 125000 |
| Initial centers | Small: (0, 0, 0.60001) m; large: (0, 0, 0) m |
| Initial translation | Small: (0, 0, -0.2) m/s; large: zero |
| Initial spin | Small omega = (0, 22.86, 0) rad/s; large: zero |
| Normal / tangential model | Hooke / linear_history; kn = 1.2e8 N/m; kt = 3.4e7 N/m |
| Restitution / sliding friction | e = 0.5; mu = 0.40 |
| Small / large mass | 11.3097335529 / 1413716.69412 kg |
| Reduced mass | 11.3096430758 kg |
| Explicit damping | cn = 15874.4669065 N s/m; ct = cr = ctorsion = 0 |
| Other forces / constraints | No gravity, wall or global damping; both spheres free; rolling/torsional resistance = 0 |
| Base step / official segments | dt = 1e-6 s; 600 + 700 + 500 = 1800 steps |
| Refined run | dt = 5e-7 s; 3600 steps; same 1.8 ms endpoint |
| Output | Every integration step, including t = 0 |
| Native material encoding | material 1: E = 10000 Pa, nu = 0.33, encoded restitution = 0.5, density = 2700 kg/m^3; material 2: E = 10000 Pa, nu = 0.33, encoded restitution = 0.5, density = 2700000 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

For a solid sphere impacting a much heavier partner in gross sliding, the official approximate component targets are [2]
$$
\Delta v_x=\mu(1+e)v_{\rm in},\qquad v_z^+=ev_{\rm in},\qquad\omega_y^+=\omega_0-\frac{5\Delta v_x}{2R}.
\tag{1}
$$

The actual second sphere has finite mass and is integrated freely. Its recoil and spin are retained; the heavy-partner formula is the official approximate target, not a replacement for the two-body dynamics.

## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

7% vx, vz and omega_y at 1.3 ms.

The independent analytic targets are \(v_x^+=0.12\,\mathrm{m/s}\), \(v_z^+=0.10\,\mathrm{m/s}\) and \(\omega_y^+=19.86\,\mathrm{rad/s}\). Each uses relative limit 0.07, evaluated at cumulative step 1300. The end of the full 1.8 ms run is supplementary.

For a native quantity \(q\) and reference \(q_{\rm ref}\), the error definitions are

$$
\varepsilon_{\rm rel}=\frac{|q-q_{\rm ref}|}{|q_{\rm ref}|},\qquad\varepsilon_{\rm abs}=|q-q_{\rm ref}|.
\tag{2}
$$

An error is accepted when it does not exceed its stated limit. Relative errors and limits are dimensionless; percentage columns equal \(100\varepsilon_{\rm rel}\). Absolute errors carry the units of the quantity. Zero-reference assertions use their stated absolute bounds rather than dividing by zero. The analytic tolerances in [1,2], or separately labeled supplemental diagnostic limits, determine acceptance. The YAML same-code regression epsilon is retained as source metadata and is not adopted as a cross-software acceptance tolerance.

## 3. Results and discussion

**Table 2.** Retained measurements for the two backends and two time steps. Numerical values, error types, limits and decisions are reproduced from the campaign record without refitting. `dt` denotes \(\Delta t\), and `dt/2` denotes the refined step; `um` denotes µm where used.

| Run | vx (m/s) | vx error (%) | vz (m/s) | vz error (%) | omega_y (rad/s) | Spin error (%) | All official 7% |
|---|---:|---:|---:|---:|---:|---:|---|
| CPU dt | 0.12706942 | 5.89118 | 0.09986565 | 0.13435 | 19.6839496 | 0.88646 | PASS |
| GPU dt | 0.12706942 | 5.89118 | 0.09986565 | 0.13435 | 19.6839496 | 0.88646 | PASS |
| CPU dt/2 | 0.12705825 | 5.88188 | 0.10004190 | 0.04190 | 19.6842301 | 0.88504 | PASS |
| GPU dt/2 | 0.12705825 | 5.88188 | 0.10004190 | 0.04190 | 19.6842301 | 0.88504 | PASS |

The large partner remains finite and free, with mass ratio 125000; its recoil and rotation are included in momentum audits. The approximately 5.9% v_x offset from the simplified heavy-partner formula persists after refinement. The retained impulse audit separates the signed normal impulse from the absolute load used by the Coulomb cap: the base-step values are 3.391438290 and 3.619502657 N s. Their difference reflects the signed unloading tail, while early history buildup and changing contact directions also affect the tangential impulse. Total linear and angular momentum drifts remain at most 1.88e-14 kg m/s and 1.89e-14 kg m²/s in the archived audit. The official analytic and saved-state comparisons remain separate.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem06-hookehistory-spin. Panels (a) and (b) show velocity (x component) and angular velocity (y component) of the selected body 'small'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows the signed residual of (b). Both residual panels use linear scales, with zero shown at its true ordinate. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

A finite freely moving large partner is retained. The test uses explicit damping, not automatic tangential damping, and is not torsion about the normal.

This report reanalyzes the retained 28 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem06-hookehistory-spin-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem06-hookehistory-spin-3d-si.yaml), commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem06_hookehistory_spin --parallel
build/validation/LAMMPSComparison/lammps_dem06_hookehistory_spin --mode cpu --output build/comparisons/dem06-hookehistory-spin-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem06-hookehistory-spin.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_checkpoint_comparison.csv](checkpoints/official_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
