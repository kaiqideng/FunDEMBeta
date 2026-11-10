# Preloaded elastic stack: LAMMPS benchmark verification

**Case:** `dem10-stack`  
**Evidence date:** 29 October 2026  
**Historical outcome:** **PASS**

## Abstract

This study examines mechanical-energy conservation for two elastic spheres confined between a floor and ceiling using pinned LAMMPS inputs and analytic assertions. At 20 ms, relative energy errors are 0.0035325% and 0.0008593% at the base and refined steps, below the official 0.05% limit. CPU and GPU calculations at the original and halved time steps give the historical outcome **PASS**. The measurements are retained campaign evidence; this revision changes their scientific presentation and does not certify a new run of the current core.

## 1. Benchmark definition

The benchmark is defined by the official input and saved reference states [1]; its analytic quantities and acceptance procedure are specified by the official checker [2]. The source commit is `29dd3fd7642a6fe9e8b4dd99120df4068b221250` and the input Git blob is `405583bd3e27294ab46506f6ba242485a756de33`. LAMMPS was not rerun locally. The pinned upstream input is tagged `unvetted`; the comparison does not remove that qualification.

**Table 1.** Physical configuration, coefficient mapping and integration settings. Values are in SI units except the explicitly mapped reduced-unit input where stated.

| Quantity | Official / native mapping |
|---|---|
| Spheres | R=0.0005 m; lower/upper density=20000/10000 kg/m^3 |
| Finite masses | Lower=1.0471975512e-05 kg; upper=5.23598775598e-06 kg |
| Initial centers | Official y=0.00045/0.00135 m -> native z; initial velocities/spins zero |
| Fixed planes | Official y=0/0.0018 m -> native z=0/0.0018 m |
| Contact law | Hooke kn=1000 N/m for all three contacts; restitution 1 gives zero automatic damping |
| Other channels | Zero tangential stiffness/friction, resistance and damping; no body damping; gravity 9.81 m/s^2 downward |
| Time | Official 2 us, five 2000-step segments to 20 ms; half-step 1 us, 20000 steps; continuous native events |
| Full output | Every integration step including initial frame: 10001/10001/20001/20001, total 60004 native frames |
| Native material encoding | material 1: E = 10000 Pa, nu = 0.33, encoded restitution = 1, density = 20000 kg/m^3; material 2: E = 10000 Pa, nu = 0.33, encoded restitution = 1, density = 10000 kg/m^3 |

The native material encoding is reported for reproducibility of the dedicated C++ input; it is not an additional fitted reference parameter. Under explicit Hooke stiffness, E and nu do not set the active spring constants. When manual contact damping is enabled, the encoded material restitution does not set the prescribed dashpot.

With the source vertical coordinate relabeled as \(z\), the official energy invariant includes both particles and each active contact once [2]:
$$
E=\sum_{i=1}^{2}\left(\frac{1}{2}m_iv_{iz}^2+m_igz_i\right)+\frac{1}{2}k_n\left(\delta_{\rm floor}^2+\delta_{12}^2+\delta_{\rm ceiling}^2\right).
\tag{1}
$$

Each overlap is set to zero when its contact is open. There is no viscous or frictional dissipation in this configuration.

## 2. Numerical method and verification procedure

The retained campaign comprises CPU and GPU calculations at \(\Delta t\) and \(\Delta t/2\). The refined runs use twice as many steps, with unchanged initial conditions and duration. Comparisons are made at equal physical times, using twice the base-step index in the refined runs. The native contact solver supplies the reported motion; the analytic reference is used for comparison, not to generate the native response.

0.05% relative total mechanical-energy error at 20 ms (official tag unvetted).

The official `stack_energy` assertion uses \(|E(20\,\mathrm{ms})-E_0|/|E_0|\leq0.0005\). The maximum whole-trace error in Table 2 is an additional diagnostic, not the endpoint assertion.

For a native quantity \(q\) and reference \(q_{\rm ref}\), the error definitions are

$$
\varepsilon_{\rm rel}=\frac{|q-q_{\rm ref}|}{|q_{\rm ref}|},\qquad\varepsilon_{\rm abs}=|q-q_{\rm ref}|.
\tag{2}
$$

An error is accepted when it does not exceed its stated limit. Relative errors and limits are dimensionless; percentage columns equal \(100\varepsilon_{\rm rel}\). Absolute errors carry the units of the quantity. Zero-reference assertions use their stated absolute bounds rather than dividing by zero. The analytic tolerances in [1,2], or separately labeled supplemental diagnostic limits, determine acceptance. The YAML same-code regression epsilon is retained as source metadata and is not adopted as a cross-software acceptance tolerance.

## 3. Results and discussion

**Table 2.** Retained measurements for the two backends and two time steps. Numerical values, error types, limits and decisions are reproduced from the campaign record without refitting. `dt` denotes \(\Delta t\), and `dt/2` denotes the refined step; `um` denotes µm where used.

| Run | Final E (microjoules) | Final relative error (%) | Maximum whole-trace error (%) | Official0.05% |
|---|---:|---:|---:|---|
| CPU dt | 7.615302319 | 0.0035325 | 0.0047021 | PASS |
| GPU dt | 7.615302319 | 0.0035325 | 0.0047021 | PASS |
| CPU dt/2 | 7.615505900 | 0.0008593 | 0.0011755 | PASS |
| GPU dt/2 | 7.615505900 | 0.0008593 | 0.0011755 | PASS |

The maximum whole-trace errors, 0.0047021% and 0.0011755%, are separate diagnostics and also decrease on refinement. Energy includes both particle kinetic energies, gravitational potential and each of the three contact spring energies once. The floor and ceiling supply external impulse; total particle momentum is not an invariant of this problem. No damping is present, and the upstream `unvetted` qualification is retained.

![Figure 1. Retained benchmark measurements](figures/comparison.svg)

**Figure 1.** Historical checkpoint comparison for dem10-stack. Panels (a) and (b) show position (z component) and velocity (z component) of the selected body 'lower'; (c) shows the signed FunDEM-minus-reference residual of (a). Black open diamonds denote pinned LAMMPS saved-state references; blue open circles, orange crosses, green open squares and purple plus signs denote CPU at Δt, GPU at Δt, CPU at Δt/2 and GPU at Δt/2, respectively. Every run uses its retained physical checkpoint times; no time shift, jitter, averaging, fitting or connecting trajectory is used. Coincident symbols are retained with different shapes and sizes. Axes state physical units and explicit powers of ten. Panel (d) shows the signed residual of (b). Both residual panels use linear scales, with zero shown at its true ordinate. The source's vertical y coordinate maps to native z; the selected body's source and native ordinates retain this documented mapping. These sparse historical checkpoints do not establish a new validation outcome for the current source or a continuous-time error bound.

## 4. Conclusions and limitations

All four retained runs satisfy the stated official analytic criteria for the active trajectory. The historical PASS is restricted to those quantities, sample times and model assumptions; it does not certify inactive contact channels or general many-particle behavior.

The source's initial precompression/bootstrap is retained. The energy criterion does not establish position, velocity or modal-phase accuracy.

This report reanalyzes the retained 29 October 2026 campaign; no new physical acceptance run is claimed. Distributed CSV checkpoint tables support the reported samples and scalar diagnostics. Full trajectories and native archives remain outside this repository, so the compact tables and figure must not be interpreted as a newly reconstructed complete time history. CPU/GPU agreement and a two-level time-step study are numerical consistency evidence, not independent proof of model equivalence or an asymptotic convergence order.

## References

[1] LAMMPS Development Team. [Official `dem10-stack-3d-si.yaml` input and saved reference states](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-stack-3d-si.yaml), commit `29dd3fd7642a6fe9e8b4dd99120df4068b221250`.

[2] LAMMPS Development Team. [Analytic DEM verification models, `test_analytic_models.cpp`](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/test_analytic_models.cpp), same pinned source revision.

## Appendix A. Reproducibility and evidence

Run from the directory containing `FunDEMBeta/`, after the [family build instructions](../../README.md#build-and-run):

```bash
cmake --build build --target lammps_dem10_stack --parallel
build/validation/LAMMPSComparison/lammps_dem10_stack --mode cpu --output build/comparisons/dem10-stack-cpu
```

Reserve a new output directory for each execution. Repeat with `--mode gpu`, `--half-step`, and their combination to obtain CPU/GPU results at both time steps; the half-step option doubles the step count to preserve duration. GPU execution requires a CUDA-enabled build and compatible device. The current dedicated C++ program provides a reproducible setup, but its completion is not a new acceptance certification. Compare its numerical output independently with the physical quantities, times and limits stated above.

- [Native C++ case](dem10-stack.cpp)
- [Exact historical numerical fields](HISTORICAL-RESULTS.csv) · [record guide](HISTORICAL-RESULTS.md)
- [Pinned sources and historical hashes](PROVENANCE.md)
- [Checkpoint: official_checkpoint_comparison.csv](checkpoints/official_checkpoint_comparison.csv)
- [Complete historical report](HISTORICAL-REPORT.md): unchanged campaign calculations, diagnostics, scope and retention detail.
