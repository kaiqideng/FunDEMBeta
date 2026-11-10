# LAMMPS comparisons

35 native C++ cases retain **33 historical PASS and 2 historical FAIL** outcomes
from 28–29 October 2026. They are not a new acceptance result for the current core.

Pinned official LAMMPS inputs, saved states and analytical criteria are cited
in each report. Most comparisons did not rerun LAMMPS locally. PASS applies
only to the stated observable and regime; the two Mindlin sliding-to-rolling
failures remain FAIL. Region-wall geometry and Hertz-oblique saturation
restrictions remain as documented in their reports.

## Read a case

Every `REPORT.md` has an abstract, benchmark definition, numerical method,
results and discussion, conclusions and limitations, references, and a
reproducibility appendix. Numbered equations, tables and self-contained figure
captions follow the [publication specification](../PUBLICATION.md). Figures
show retained campaign measurements; they do not imply new executions. The unchanged original
analysis is available separately as `HISTORICAL-REPORT.md`.

| File | Purpose |
| --- | --- |
| `REPORT.md` | Publication-format comparison, equations, criterion, figure and scope |
| `HISTORICAL-REPORT.md` | Byte-preserved report before presentation cleanup |
| `HISTORICAL-RESULTS.csv` | Exact numerical fields and original pass flags |
| `PROVENANCE.md` | Source pins, campaign hashes and execution inventory |
| `checkpoints/*.csv` | Available compact saved-state comparisons |
| `<case>.cpp` | Standalone native-core reproduction input |

The [61-input coverage matrix](diagnostics/reference/lammps/20261029-expanded-suite/coverage-matrix.md)
separates completed comparisons from unexecuted or excluded inputs. The earlier
24-input inventory is a frozen source snapshot, not the current coverage count.

## Build and run

Run from the directory containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_VALIDATIONS=ON
cmake --build build --target lammps_dem01_hertz_twosphere --parallel
build/validation/LAMMPSComparison/lammps_dem01_hertz_twosphere --mode cpu --output build/comparisons/dem01-hertz-twosphere-cpu
```

`lammpsComparison` builds all family targets; it is not an executable.
Each case is fully defined in its C++ source. No Qt, desktop project/config
loader or reference-program installation is needed. GPU/Hybrid require a
CUDA-enabled build and compatible device. Use `--half-step` for dt/2 and
preserved default duration; a positive `--steps` limit gives an INCOMPLETE
short run. Use a new output directory for every invocation. See
[shared options and output](../common/README.md).

## Case index

| Case | Historical outcome | Report | C++ source |
| --- | --- | --- | --- |
| `dem01-hertz-tsuji` | **PASS** | [Report](cases/dem01-hertz-tsuji/REPORT.md) | [Source](cases/dem01-hertz-tsuji/dem01-hertz-tsuji.cpp) |
| `dem01-hertz-twosphere` | **PASS** | [Report](cases/dem01-hertz-twosphere/REPORT.md) | [Source](cases/dem01-hertz-twosphere/dem01-hertz-twosphere.cpp) |
| `dem01-hooke-2d-lj` | **PASS** | [Report](cases/dem01-hooke-2d-lj/REPORT.md) | [Source](cases/dem01-hooke-2d-lj/dem01-hooke-2d-lj.cpp) |
| `dem01-hooke-massvel` | **PASS** | [Report](cases/dem01-hooke-massvel/REPORT.md) | [Source](cases/dem01-hooke-massvel/dem01-hooke-massvel.cpp) |
| `dem01-hooke-twosphere` | **PASS** | [Report](cases/dem01-hooke-twosphere/REPORT.md) | [Source](cases/dem01-hooke-twosphere/dem01-hooke-twosphere.cpp) |
| `dem01-hooke-twosphere-e02` | **PASS** | [Report](cases/dem01-hooke-twosphere-e02/REPORT.md) | [Source](cases/dem01-hooke-twosphere-e02/dem01-hooke-twosphere-e02.cpp) |
| `dem01-hooke-twosphere-e08` | **PASS** | [Report](cases/dem01-hooke-twosphere-e08/REPORT.md) | [Source](cases/dem01-hooke-twosphere-e08/dem01-hooke-twosphere-e08.cpp) |
| `dem01-hooke-veldamp` | **PASS** | [Report](cases/dem01-hooke-veldamp/REPORT.md) | [Source](cases/dem01-hooke-veldamp/dem01-hooke-veldamp.cpp) |
| `dem01-hookehistory-twosphere` | **PASS** | [Report](cases/dem01-hookehistory-twosphere/REPORT.md) | [Source](cases/dem01-hookehistory-twosphere/dem01-hookehistory-twosphere.cpp) |
| `dem01-legacy-hooke-twosphere` | **PASS** | [Report](cases/dem01-legacy-hooke-twosphere/REPORT.md) | [Source](cases/dem01-legacy-hooke-twosphere/dem01-legacy-hooke-twosphere.cpp) |
| `dem01-legacy-hookehistory-twosphere` | **PASS** | [Report](cases/dem01-legacy-hookehistory-twosphere/REPORT.md) | [Source](cases/dem01-legacy-hookehistory-twosphere/dem01-legacy-hookehistory-twosphere.cpp) |
| `dem02-hertz-peak` | **PASS** | [Report](cases/dem02-hertz-peak/REPORT.md) | [Source](cases/dem02-hertz-peak/dem02-hertz-peak.cpp) |
| `dem02-hertz-twosphere` | **PASS** | [Report](cases/dem02-hertz-twosphere/REPORT.md) | [Source](cases/dem02-hertz-twosphere/dem02-hertz-twosphere.cpp) |
| `dem02-hertzmaterial-twosphere` | **PASS** | [Report](cases/dem02-hertzmaterial-twosphere/REPORT.md) | [Source](cases/dem02-hertzmaterial-twosphere/dem02-hertzmaterial-twosphere.cpp) |
| `dem02-legacy-hertzhistory-twosphere` | **PASS** | [Report](cases/dem02-legacy-hertzhistory-twosphere/REPORT.md) | [Source](cases/dem02-legacy-hertzhistory-twosphere/dem02-legacy-hertzhistory-twosphere.cpp) |
| `dem03-hertz-peak-wall` | **PASS** | [Report](cases/dem03-hertz-peak-wall/REPORT.md) | [Source](cases/dem03-hertz-peak-wall/dem03-hertz-peak-wall.cpp) |
| `dem03-hertz-tsuji-wall` | **PASS** | [Report](cases/dem03-hertz-tsuji-wall/REPORT.md) | [Source](cases/dem03-hertz-tsuji-wall/dem03-hertz-tsuji-wall.cpp) |
| `dem03-hertz-wall` | **PASS** | [Report](cases/dem03-hertz-wall/REPORT.md) | [Source](cases/dem03-hertz-wall/dem03-hertz-wall.cpp) |
| `dem03-hertzmaterial-wall` | **PASS** | [Report](cases/dem03-hertzmaterial-wall/REPORT.md) | [Source](cases/dem03-hertzmaterial-wall/dem03-hertzmaterial-wall.cpp) |
| `dem03-legacy-hertzhistory-wall` | **PASS** | [Report](cases/dem03-legacy-hertzhistory-wall/REPORT.md) | [Source](cases/dem03-legacy-hertzhistory-wall/dem03-legacy-hertzhistory-wall.cpp) |
| `dem04-hertz-oblique` | **PASS** | [Report](cases/dem04-hertz-oblique/REPORT.md) | [Source](cases/dem04-hertz-oblique/dem04-hertz-oblique.cpp) |
| `dem04-hookehistory-oblique` | **PASS** | [Report](cases/dem04-hookehistory-oblique/REPORT.md) | [Source](cases/dem04-hookehistory-oblique/dem04-hookehistory-oblique.cpp) |
| `dem05-hookehistory-slip` | **PASS** | [Report](cases/dem05-hookehistory-slip/REPORT.md) | [Source](cases/dem05-hookehistory-slip/dem05-hookehistory-slip.cpp) |
| `dem05-mindlin-slip` | **FAIL** | [Report](cases/dem05-mindlin-slip/REPORT.md) | [Source](cases/dem05-mindlin-slip/dem05-mindlin-slip.cpp) |
| `dem05-slip-mu02` | **FAIL** | [Report](cases/dem05-slip-mu02/REPORT.md) | [Source](cases/dem05-slip-mu02/dem05-slip-mu02.cpp) |
| `dem05-sliptime` | **PASS** | [Report](cases/dem05-sliptime/REPORT.md) | [Source](cases/dem05-sliptime/dem05-sliptime.cpp) |
| `dem06-hookehistory-spin` | **PASS** | [Report](cases/dem06-hookehistory-spin/REPORT.md) | [Source](cases/dem06-hookehistory-spin/dem06-hookehistory-spin.cpp) |
| `dem06-spin-identical` | **PASS** | [Report](cases/dem06-spin-identical/REPORT.md) | [Source](cases/dem06-spin-identical/dem06-spin-identical.cpp) |
| `dem07-hooke-rolling` | **PASS** | [Report](cases/dem07-hooke-rolling/REPORT.md) | [Source](cases/dem07-hooke-rolling/dem07-hooke-rolling.cpp) |
| `dem09-viscous-linear` | **PASS** | [Report](cases/dem09-viscous-linear/REPORT.md) | [Source](cases/dem09-viscous-linear/dem09-viscous-linear.cpp) |
| `dem10-freefall` | **PASS** | [Report](cases/dem10-freefall/REPORT.md) | [Source](cases/dem10-freefall/dem10-freefall.cpp) |
| `dem10-stack` | **PASS** | [Report](cases/dem10-stack/REPORT.md) | [Source](cases/dem10-stack/dem10-stack.cpp) |
| `dem11-regionwall` | **PASS** | [Report](cases/dem11-regionwall/REPORT.md) | [Source](cases/dem11-regionwall/dem11-regionwall.cpp) |
| `dem13-twist-sds` | **PASS** | [Report](cases/dem13-twist-sds/REPORT.md) | [Source](cases/dem13-twist-sds/dem13-twist-sds.cpp) |
| `dem15-hooke-oblique` | **PASS** | [Report](cases/dem15-hooke-oblique/REPORT.md) | [Source](cases/dem15-hooke-oblique/dem15-hooke-oblique.cpp) |

## Historical evidence

Each case keeps its concise report, unchanged historical report, numerical
record, provenance and available compact checkpoint tables.
[HISTORICAL-INDEX.md](HISTORICAL-INDEX.md) and
[HISTORICAL-INDEX.csv](HISTORICAL-INDEX.csv) retain the case identities/outcomes.
Original numeric tokens, tolerances, reference citations and evidence hashes
are preserved. Private artifact locations and application identities are omitted.
Original JSON and result images are backed up outside the core source tree.
No historical report is upgraded by successful compilation or stepping.

`EVIDENCE-SHA256SUMS.txt` describes the current distributed evidence.
`HISTORICAL-EVIDENCE-SHA256SUMS.txt` preserves the earlier checksum ledger with
its original path labels. Its CSV checksums remain valid; its 35 report
checksums were already stale before this cleanup. Archived report hashes are
identified separately in each `PROVENANCE.md`. See [evidence integrity](EVIDENCE.md)
for the distinction; the old digests have not been rewritten.

## Figures and integrity

Every report includes one multi-panel checkpoint figure at 180 mm width,
supplied as vector PDF, editable SVG and a 600 dpi PNG fallback. Embedded fonts,
physical units, distinct markers and self-contained captions support publication
use. `FIGURES.json` records sources, selected rows and fields, dimensions and
artifact hashes. Figures use only retained measurements. Sparse
checkpoint markers are not interpolated into full trajectories. Error plots
divide each recorded error by its own stated tolerance; the boundary is 1.

From this directory, with Python 3 and matplotlib available:

```bash
python3 tools/plotCheckpoints.py --write
python3 tools/plotCheckpoints.py --check
sha256sum -c EVIDENCE-SHA256SUMS.txt
```

The drawing script does not execute LAMMPS or FunDEM, fit parameters, or change
historical outcomes. After an intentional presentation update, refresh the
current ledger with the [shared documentation audit](../common/README.md#documentation-checks).

Excluded mappings and source-law diagnostics remain in [diagnostics](diagnostics/LAMMPS-DIAGNOSTICS.md).
These are distinct from the 35-case count and its two historical official FAILs.
