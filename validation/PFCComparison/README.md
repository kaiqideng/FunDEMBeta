# PFC comparisons

Three native C++ cases retain **3 historical PASS** outcomes from 28–29 October
2026. **PFC was not executed locally.** PFC is proprietary Itasca software;
the targets here come from its public documentation and analytical relations.
These records are not a fresh acceptance result for the current core.

| Case | Scope | Historical outcome |
| --- | --- | --- |
| [Damping ratio 0.2](cases/pfc-linear-restitution-md0-beta02/REPORT.md) | First rebound, signed unloading and impulse | PASS |
| [Damping ratio 0.6](cases/pfc-linear-restitution-md0-beta06/REPORT.md) | First rebound and static settling | PASS |
| [Sixteen-ratio sweep](cases/pfc-linear-restitution-md0-sweep/REPORT.md) | Complete Md=0 row; 64 ball/run checks | PASS |

The criterion is absolute first-apex restitution error **<0.01**, not a universal
relative 1% tolerance. The sweep does not cover the Md=1 no-tension row or the
complete 32-ball example. [Coverage](COVERAGE.md) documents all fifteen official
verification-index entries and the excluded operators.

## Build and run

Keep the terminal in the directory containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_VALIDATIONS=ON
cmake --build build --target pfc_pfc_linear_restitution_md0_beta02 --parallel
build/validation/PFCComparison/pfc_pfc_linear_restitution_md0_beta02 --mode cpu --output build/comparisons/pfc-beta02-cpu
```

`pfcComparison` builds all three executables; it does not run them. Each case is
defined in `cases/<id>/<id>.cpp` and links the public core without Qt, project
loading or a PFC installation. GPU/Hybrid require a CUDA build and device.
Use `--half-step` for Δt/2, or `--steps N` for an explicitly INCOMPLETE short run.
Always choose a new output directory. See [shared options](../common/README.md).

## Reports and evidence

Each `REPORT.md` has an abstract, benchmark definition, numerical method,
results and discussion, conclusions and limitations, references, and a
reproducibility appendix, following the [publication specification](../PUBLICATION.md);
`HISTORICAL-REPORT.md` retains the detailed original report. Checkpoint and
`HISTORICAL-RESULTS.csv` values are unchanged. `PROVENANCE.md` retains original
source and archive hashes; its archived-report hash identifies the detailed
report, not a rewritten presentation. The
[historical index](HISTORICAL-INDEX.md) preserves case identities and outcomes.
Full native histories, proprietary programs and vendor data are not included.

Figures are generated only from the retained checkpoint tables, with no fitted
or synthetic trajectories. They are supplied as vector PDF, editable SVG and
600 dpi PNG at 180 mm width, with embedded scientific type, physical units and
self-contained captions. `FIGURES.json` records exact CSV sources, selected
rows/fields, dimensions and artifact hashes.
Recreate them with Python and matplotlib:

```bash
python3 FunDEMBeta/validation/PFCComparison/plotCheckpoints.py
python3 FunDEMBeta/validation/PFCComparison/plotCheckpoints.py --check
```

Reformatting, regenerating a figure, successful compilation or completing a new
solve never upgrades the historical outcome. See
[evidence integrity](EVIDENCE.md) for the documentation-only revision policy.
