# Shared comparison utilities

Each comparison executable defines its own model in C++. This library supplies
CLI/output, explicit check evaluation and prepared native geometry only; it does
not load desktop projects or select a model by ID. It links `FunDEM::FunDEM`
without Qt. The ON-DEM adapter has its own CLI, documented in its
[case guide](../ONDEMComparison/shearFlow/README.md).

## Run one comparison

From the directory containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_VALIDATIONS=ON
cmake --build build --target lammps_dem01_hertz_twosphere --parallel
build/validation/LAMMPSComparison/lammps_dem01_hertz_twosphere --mode cpu --output build/comparisons/dem01-cpu
```

| Option | Meaning |
| --- | --- |
| `--mode cpu\|gpu\|hybrid` | Native backend; unavailable CUDA modes fail explicitly |
| `--half-step` | Δt/2, twice the default steps and twice the checkpoint indices |
| `--steps N` | Positive step cap; a shortened run is INCOMPLETE |
| `--output NEW_DIR` | Fresh directory; any existing directory is rejected |
| `--help` | Usage only; no geometry building or solve |

With no arguments, the executable runs its own case into a fresh
`<id>-results/` under the current directory. There is no `--case`, `--list` or
input-path argument. For multi-config Windows builds use `.exe` and the
configuration subdirectory. A serial MSVC CPU build can use
`-DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_ENABLE_OPENMP=OFF`; classic MSVC OpenMP does
not support every grid reduction. `lammpsComparison` and `pfcComparison` build
their families; they never start a suite.

## Output and interpretation

| File | Contents |
| --- | --- |
| `trajectory.csv` | Every-step state, force, torque and translational kinetic energy |
| `contacts.csv` | Family/indices, overlap, signed normal force, force and elastic energy |
| `comparison.csv` | Expected/actual checkpoints, error definitions and tolerances |
| `summary.txt` | Steps, fresh verdict and separately labeled historical outcome |

Translational energy excludes rotation. Signed normal-force projections are
distinct from the aggregate normal-force-magnitude checkpoint. Use each
report's stated observable and regime, not just similar column names.

Eighteen of the 38 LAMMPS/PFC cases configure 38 automatic source metrics; the
other twenty report **NOT EVALUATED**, not PASS. Short runs are **INCOMPLETE**.
Historical PASS/FAIL is never substituted for the fresh verdict.

| Exit code | Meaning |
| --- | --- |
| 0 | Completion or intentional short run; not necessarily a comparison PASS |
| 1 | Input, runtime or IO error |
| 2 | Full run failed a configured acceptance gate |

Always compare a new run independently and retain it outside the source tree.

## Documentation checks

Reports follow the [publication specification](../PUBLICATION.md): abstract,
benchmark definition, numerical method, results and discussion, conclusions and
limitations, references, and a reproducibility appendix. Exact numerical tables
and archived reports remain unchanged. PDF and SVG figures preserve vectors;
600 dpi PNG exports use the same 180 mm final width.

From the directory containing `FunDEMBeta/`:

```bash
python3 FunDEMBeta/validation/common/auditDocumentation.py
```

The standard-library audit checks local document/image links, preserved CSV
hashes, archived and current report hashes, case counts/outcomes, publication
sections, vector PDF/static SVG, PNG dimensions and resolution, figure provenance
and current checksum inventories. It does not run a solver. An intentional
documentation or figure update can refresh current ledgers with
`--write-manifests`; historical ledgers and numerical data are never rewritten.
