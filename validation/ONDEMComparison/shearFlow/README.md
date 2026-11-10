# Granular shear flow between bumpy walls: model definition and unresolved reference discrepancy

## Abstract

This report describes a CPU adaptation of the ON-DEM benchmark for monodisperse spheres sheared between two periodic bumpy walls [1,2]. The implementation uses native sphere storage, mass/inertia, velocity-Verlet integration, solver time and VTU output, with a case-specific external contact-force component. Historical transient comparisons from 17 October 2026 show mobile-particle kinetic-energy differences of −51.077875% and −64.970898% using the published coefficients at matched sample times. A separately identified diagnostic coefficient set reduces those two differences to −0.323035% and +0.487018%, but does not reproduce trajectories or instantaneous contact pairs. The full 200-million-step benchmark, time-step study and steady continuum-field comparison have **not** been completed. No validation PASS or general GPU/Hybrid capability is claimed.

## 1. Benchmark definition

The benchmark specifies a constant-volume, gravity-free assembly of identical spheres driven by two bumpy boundaries moving in opposite streamwise directions [1]. Quantities are dimensionless under $\rho_p=d=V=1$; they are not calibrated SI material properties. The original source IDs, particle order and coordinates are preserved when the input is explicitly retrieved.

**Table 1. Published benchmark inputs.** These remain the runner defaults. The diagnostic set in Table 2 is separate and does not replace them.

| Quantity | Symbol | Published/default value |
| --- | --- | --- |
| Mobile spheres | $N$ | 3000 |
| Prescribed wall spheres | $N_w$ | 240 bottom and 240 top |
| Diameter; density | $d$; $\rho_p$ | 1; 1 |
| Domain | $(L_x,L_y,L_z)$ | $(20,10.3923048454133,22)$ |
| Periodicity | — | X/Y periodic; Z nonperiodic |
| Bottom/top wall-sphere center | $z_\mathrm{bot}$; $z_\mathrm{top}$ | 0.5; 21.5 |
| Clear wall gap | $H$ | 20 |
| Bottom/top prescribed velocity | $\mathbf{V}_\mathrm{bot}$; $\mathbf{V}_\mathrm{top}$ | $(−1,0,0)$; $(1,0,0)$ |
| Gravity | $\mathbf{g}$ | Zero |
| Pair stiffness | $k_n$; $k_t$ | $2\times10^5$; $5.7143\times10^4$ |
| Direct dashpot coefficients | $c_n$; $c_t$ | 50; 14.2857 |
| Sliding friction | $\mu$ | 0.3 |
| Time step | $\Delta t$ | $7.2320\times10^{-5}$ |
| Full step count | $N_t$ | 200000000 |
| Output interval | $N_s$ | 50000 steps |

The external component owns all five pair coefficients. Nominal native material data supplies mass and inertia, not an additional contact law or restitution-to-dashpot conversion. For an overlap $\delta$, minimum-image contact normal $\mathbf{n}$ and normal relative contact velocity $v_n$, the adopted direct-force law is

$$
F_n=\max(0,k_n\delta-c_n v_n),\qquad \mathbf{F}_{t,\mathrm{trial}}=-k_t\boldsymbol{\xi}_t-c_t\mathbf{v}_t.\tag{1}
$$

The spring is transported into the current tangent plane, integrated once per history advance and returned after capping the **total** elastic-plus-damping tangential force:

$$
\|\mathbf{F}_t\|\leq\mu F_n.\tag{2}
$$

This explicit return differs from the default core law, which limits damping separately. Rolling and torsional resistance are absent. The page does not uniquely specify all force/history conventions; noncohesive normal truncation, spring transport and Coulomb return must be compared with the author's actual implementation before claiming exact cross-code equivalence.

Equation (1) is this adapter's default `clampNonAttractive` release convention,
not a description of the current native contact implementation. The latter
currently leaves its normal-force truncation disabled. The adapter's optional
`retainUntilSeparation` policy is a separate C++ diagnostic setting; neither
the published/default command nor the retained transient outcomes are silently
switched to that policy by a native-core change.

## 2. Numerical method and verification procedure

One linked-cell periodic search is rebuilt each step, evaluating each unordered pair once. Minimum images determine separation, normals, contact lever arms, relative velocities and torques. Tangential histories use stable particle-pair indices across boundary wrapping. Open contacts lose their previous history. Wall spheres retain prescribed velocities and disabled accelerations; wall–wall forces are omitted.

`periodicContactForce.h/.cpp` adds loads through the external-force hook without clearing unrelated loads, moving particles or advancing solver time. Initialization advances history by zero; each subsequent velocity-Verlet step advances it once by the current integration step. Native nonperiodic contact evaluation is not also run. The isolated CPU adaptation is not a general GPU/Hybrid solver benchmark.

The retained implementation checks cover periodic mechanics, Mercury-format IO and additive external-force assembly. A 120-particle linked-cell/brute-force comparison is included in the historical checks. The frozen 3480-particle input passed initialization and a 100-step smoke run. A separate 1000-step extraction comparison preserved particle/contact files byte-for-byte at steps 0, 500 and 1000. These checks establish specified implementation behavior, **not** agreement with the full ON-DEM benchmark.

Two 100000-step transient comparisons were executed on 17 October 2026. Both used the snapshot-inferred step $7.192545032\times10^{-5}$ to match reference timestamps at steps 50000 and 100000. Thus the published-coefficient transient comparison is not a run at the page's default time step. Its initial coordinates, radii and translational velocities match the reference; negligible initial angular velocities differ by at most $3.970502315\times10^{-12}$.

For mobile particles only, the instantaneous kinetic energy and signed comparison error are

$$
K_m(t)=\sum_{i\in\mathrm{mobile}}\left(\frac{1}{2}m_i\|\mathbf{v}_i\|^2+\frac{1}{2}\boldsymbol{\omega}_i^{\mathsf{T}}\mathbf{I}_i\boldsymbol{\omega}_i\right),\qquad \varepsilon_K(t)=100\frac{K_m^{\mathrm{FunDEM}}(t)-K_m^{\mathrm{reference}}(t)}{K_m^{\mathrm{reference}}(t)}\%.\tag{3}
$$

Prescribed wall grains are excluded. Matching this scalar cannot establish trajectory or stress-profile equivalence. No official quantitative acceptance threshold is supplied here, and no new tolerance is introduced to convert these comparisons into a PASS.

**Table 2. Published defaults versus explicitly inferred/diagnostic settings.** The inferred set is not an author-confirmed correction; tangential history/damping cannot be uniquely recovered from the available snapshots.

| Quantity | Published/default | Snapshot-inferred or diagnostic | Evidential status |
| --- | ---: | ---: | --- |
| $\Delta t$ | $7.2320\times10^{-5}$ | $7.192545032\times10^{-5}$ | Inferred from example time / 50000 |
| $c_n$ | 50 | 15 | Conditional force audit, assuming velocity Verlet |
| $c_t$ | 14.2857 | $30/7$ | Diagnostic choice; not uniquely inferred |
| $k_t$ | 57143 | $400000/7$ | Diagnostic exact ratio; not author-confirmed |
| $k_n$; $\mu$ | 200000; 0.3 | Unchanged | Preserved in both runs |

## 3. Results and discussion

**Table 3. Historical instantaneous mobile-particle kinetic-energy comparison at matched reference times (17 October 2026).** The published-coefficient and diagnostic runs both use the snapshot-inferred time step. Differences are signed relative to the author's example, not errors against an experimental observation or completed steady benchmark.

| Normalized time | Reference $K_m$ | Published-coefficient FunDEM $K_m$ | Published-coefficient $\varepsilon_K$ | Diagnostic $\varepsilon_K$ |
| --- | ---: | ---: | ---: | ---: |
| 3.596272516 | 5.20357556 | 2.54569974 | −51.077875% | −0.323035% |
| 7.192545032 | 8.47848212 | 2.96993614 | −64.970898% | +0.487018% |

The unresolved differences using published coefficients are approximately 51% and 65%. Smaller diagnostic energy differences do not provide an exact snapshot reproduction: individual trajectories and instantaneous contact pairs still differ. The defaults remain unchanged. Neither run is a steady-state validation or calibrated pass, and no new curve is interpolated from these two scalar observations.

The supplied sample files contain three complete 3480-particle snapshots at times 0, 3.596272516 and 7.192545032. Independent geometric search finds zero, one and three physical contacts, corresponding exactly to zero, two and five exported contact records. Their small contact files are therefore not evidence of truncation.

Assuming the reference uses velocity Verlet, the force-evaluation velocity is recovered by subtracting the final half kick, $\mathbf{v}_\mathrm{assembly}=\mathbf{v}_\mathrm{output}-\Delta t\mathbf{F}/(2m)$, for mobile grains. With all exported forces included, four unique contacts imply

$$
c_n=\frac{k_n\delta-F_n}{v_{n,\mathrm{assembly}}}=15.000000\quad\text{to within }6.5\times10^{-8}.\tag{4}
$$

This conditional result and the snapshot-inferred step support a discrepancy between the example snapshots and the published parameter table. The missing author simulation script prevents confirmation. The tangential dashpot and history update cannot be uniquely identified from these scalar elongation samples. With $m_\mathrm{eff}=\pi/12$ and $k_n=2\times10^5$, the collision-time estimate divided by 50 is approximately $7.2319727693\times10^{-5}$ for $c_n=50$, consistent with the page's rounded default; this direct-force convention is a documented inference, not an author-confirmed interpretation.

No result figure is supplied for this unresolved case. The retained numerical comparison is Table 3; a fitted history, late-time continuum profile or visual validation result cannot be inferred from it.

## 4. Conclusions and limitations

The isolated CPU implementation retains reproducible periodic mechanics, native integration and discrete-output checks. Published-coefficient transients do not agree with the supplied example energies, and a conditional force audit identifies a likely normal-damping inconsistency. The diagnostic run reduces two scalar energy differences without resolving model provenance or reproducing individual trajectories.

The **full 200-million-step comparison has not been run**. Time-step convergence and the benchmark's late-time averaged streamwise velocity, concentration, granular temperature, pressure and shear-stress profiles remain unverified. The required coarse-graining definition and 100 vertical slices must be followed. The page's phrase about X/Z averaging conflicts with a vertical profile; the explicit interpretation in the retained notes is averaging in the periodic horizontal X/Y directions. Instantaneous particle-center bins are not MercuryCG fields and cannot substitute for the required late-time analysis.

The pinned upstream repository supplies no redistribution license or simulation script. Original coordinates are therefore **not bundled**; retrieval remains explicit, TLS-verified, commit-pinned and SHA256-checked. Successful compilation, initialization, a short solve or exported visualization is insufficient to claim full validation.

## References

1. ON-DEM. *Granular shear flows between bumpy walls case*. Public benchmark definition, lead Dalila Vescovi. [Benchmark page](https://on-dem.atlassian.net/wiki/spaces/Index/pages/244744204/Granular+shear+flows+between+bumpy+walls+case). Accessed 31 October 2026.
2. Vescovi, D. *GranularShearFlows*: initial geometry and example outputs, commit `56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b`. [Pinned author repository](https://github.com/dalilavescovi/GranularShearFlows/tree/56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b). No redistribution license or solver script is supplied in the retained source inventory.

## Appendix A. Reproducibility and evidence

Keep the terminal in the directory containing `FunDEMBeta/`. Build the CPU adaptation, run its three small checks and explicitly retrieve the frozen initial geometry:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_BUILD_VALIDATIONS=ON
cmake --build build --target ondemShearFlow checkOndemShearFlow --parallel
cmake -DDESTINATION=build/validation/ondemShearFlow_input/GranularShearFlow.data -P FunDEMBeta/validation/ONDEMComparison/shearFlow/fetchInitialConfiguration.cmake
build/validation/ondemShearFlow --check
build/validation/ondemShearFlow --steps 100 --output-interval 50 --vtu --output build/ondem-short-run
```

No steps run by default. `checkOndemShearFlow` downloads nothing and does not start the full benchmark. Nonempty output directories are refused. Multi-configuration builds require `--config Release`, the executable's configuration subdirectory and explicit `--input FILE` if the input is outside its default location.

Reproduce the historical matched-time and diagnostic settings in separate new directories:

```bash
OMP_NUM_THREADS=2 build/validation/ondemShearFlow --steps 100000 --output-interval 50000 --time-step 7.192545032e-5 --output build/ondem-example-time-run
OMP_NUM_THREADS=2 build/validation/ondemShearFlow --steps 100000 --output-interval 50000 --time-step 7.192545032e-5 --normal-damping 15 --sliding-damping 4.285714285714286 --sliding-stiffness 57142.857142857145 --output build/ondem-snapshot-coefficients-run
```

The full duration is opt-in: `--steps 200000000 --output-interval 50000`. These instructions do not assert that it has been completed. `--vtu` adds visualization files, not another force calculation; `--no-output` disables frame writing for timing measurements.

**Table A1. Discrete outputs and comparison tools.** None of these outputs alone establishes steady-state validation.

| Output or tool | Meaning |
| --- | --- |
| `ShearFlow.data.NNNN` | Ordered mobile/prescribed particle states, stable source IDs |
| `ShearFlow.fstat.NNNN` | Actual contact records; mobile/mobile written in both directions, mobile/wall once, no fixed/fixed or dummy records |
| `--vtu` | Optional native particle visualization |
| `compareReference.py` | Time/row-checked transient metrics and instantaneous center-bin profiles |
| `auditReference.py` | Pinned-snapshot periodic-contact and conditional force audit |

`compareReference.py` requires numpy and matplotlib; the C++ program does not. After obtaining the six pinned reference particle/contact files, the commands and SHA256 checks in [HISTORICAL-NOTES.md](HISTORICAL-NOTES.md#sources-and-frozen-input) reproduce the retained audits. The notes retain original numerical details, execution timing, field conventions and coarse-graining caveats; this report changes presentation and evidential clarity, not those records.
