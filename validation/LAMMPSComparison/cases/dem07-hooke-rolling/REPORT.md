# Rolling resistance: LAMMPS benchmark verification

**Case:** `dem07-hooke-rolling`  
**Evidence date:** 28 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines Coulomb-capped rolling-resistance torque and spin decay about an axis tangent to a plane using pinned LAMMPS inputs and analytic assertions. At 7 ms, the base-step and refined spins are 0.6499413004 and 0.6499285592 rad/s against the official 0.65 rad/s target, within 1%. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `f9720ce2d176d615495d1d6c31349b1f71be1c73` and the input Git blob is `70c8d8237edc01e7b207a9a1724b416118e179b0`. LAMMPS was not rerun locally. The pinned upstream inputs and their original qualifications are retained.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official input / FunDEM configuration |
|---|---|
| Sphere diameter; density | 0.005 m; 2700 kg/m^3 |
| Mass; solid-sphere inertia | 0.000176714586764 kg; 4.41786466911e-10 kg m^2 |
| Initial center; translation | (0, 0, 0.002499982) m; zero translation |
| Initial angular velocity | (0, 1, 0) rad/s |
| Wall; gravity | Fixed plane z = 0; g = 1 m/s^2 downward |
| Normal model; stiffness | Hooke; kn = 10000 N/m |
| Sliding channel | Zero sliding stiffness, damping and friction |
| Rolling channel | Spring-dashpot-slider; kr = 10000 N/m, damping = 0, mu_r = 0.05 |
| Restitution; torsion; global damping | e = 1; torsional resistance = 0; global damping = 0 |
| Time step; duration | dt = 1e-6 s; 10000 steps; 0.01 s |
| Official segments | 3000 + 4000 + 3000 steps; unchanged continuous motion |
| Refined run | dt = 5e-7 s; 20000 steps; same 0.01 s endpoint |
| Native saved result | Every 10 us including both endpoints; 1001 frames per run |
| Detailed trajectory | Every integration step in CSV |
| Native material encoding | E = 10000 Pa, nu = 0.33, encoded restitution = 1, density = 2700 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

In the saturated rolling regime, the torque opposes spin. For the ideal constant normal load \(F_n=mg\), the official target is [2]
$$
M_y=-\mu_r RF_n,\qquad I=\frac{2}{5}mR^2,\qquad\omega_y(t)=\omega_0-\frac{5\mu_rg}{2R}t.
\tag{1}
$$

The actual torque follows the instantaneous contact load; the ideal constant-load line and the native load oscillation are distinct comparisons.

## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

1% omega_y at 7 ms.

The official `rolling_decay` target is \(\omega_y(7\,\mathrm{ms})=0.65\,\mathrm{rad/s}\) with relative limit 0.01, at cumulative step 7000. The 10 ms endpoint is a separate diagnostic and does not replace that check.

For a native quantity \(q\) and reference \(q_{\rm ref}\), the error definitions are

$$
\varepsilon_{\rm rel}=\frac{|q-q_{\rm ref}|}{|q_{\rm ref}|},\qquad\varepsilon_{\rm abs}=|q-q_{\rm ref}|.
\tag{2}
$$

An error is accepted when it does not exceed its stated limit. Relative errors and limits are dimensionless; percentage columns equal \(100\varepsilon_{\rm rel}\). Absolute errors carry the units of the quantity. Zero-reference assertions use their stated absolute bounds rather than dividing by zero. The analytic tolerances in [1,2], or separately labeled supplemental diagnostic limits, determine acceptance. The YAML same-code regression epsilon is retained as source metadata and is not adopted as a cross-software acceptance tolerance.

## 3. Results and discussion

**Table 2.** Retained measurements for the two backends and two time steps. Numerical values, error types, limits and decisions are reproduced from the campaign record without refitting. `dt` denotes \(\Delta t\), and `dt/2` denotes the refined step; `um` denotes µm where used.

| Run | dt (s) | Steps / native frames | omega_y at 7 ms (rad/s) | Official error (%) | omega_y at 10 ms (rad/s) | Official check |
|---|---:|---:|---:|---:|---:|---|
| CPU dt | 1e-06 | 10,000 / 1,001 | 0.6499413004 | 0.009031 | 0.5000467081 | PASS |
| GPU dt | 1e-06 | 10,000 / 1,001 | 0.6499413004 | 0.009031 | 0.5000467081 | PASS |
| CPU dt/2 | 5e-07 | 20,000 / 1,001 | 0.6499285592 | 0.010991 | 0.5000339920 | PASS |
| GPU dt/2 | 5e-07 | 20,000 / 1,001 | 0.6499285592 | 0.010991 | 0.5000339920 | PASS |

The official errors are 0.009031% and 0.010991%. Their slight increase under refinement reflects the comparison with a constant-load idealization, not failure of the 1% assertion. The initial overlap differs slightly from mg/k_n, leaving an undamped normal-load oscillation of about 1.86%. The archived torque audit verifies the instantaneous cap M_y = −μ_r R F_n, with maximum coefficient error 6.94e-18 and spin–angular-impulse residual 1.1e-16 rad/s. The 10 ms spin is a supplementary diagnostic. Sliding force is zero, so this test isolates rolling resistance rather than sphere–wall tangential-force torque.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem07-hooke-rolling. Panels (a) and (b) show angular velocity (y component) and torque (y component) of the selected body 'sphere'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows the signed residual of (b). Both residual panels use linear scales, with zero shown at its true ordinate. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

Capped rolling torque and spin decay are tested. Sliding friction is zero; torsion and translating rolling motion are outside scope.

This report reanalyzes the retained 28 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem07-hooke-rolling-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem07-hooke-rolling-3d-si.yaml), commit `f9720ce2d176d615495d1d6c31349b1f71be1c73`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem07_hooke_rolling --parallel
build/validation/LAMMPSComparison/lammps_dem07_hooke_rolling --mode cpu --output build/comparisons/dem07-hooke-rolling-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem07-hooke-rolling.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_checkpoint_comparison.csv](checkpoints/official_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
