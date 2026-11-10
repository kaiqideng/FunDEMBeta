# Reference comparisons

| Family | Included scope | Historical status |
| --- | --- | --- |
| [LAMMPS](LAMMPSComparison/README.md) | 35 independent native C++ reference cases | 33 PASS; 2 FAIL retained |
| [PFC](PFCComparison/README.md) | Three Md=0 restitution cases against public analytical targets | 3 PASS; no local PFC execution |
| [ON-DEM](ONDEMComparison/README.md) | Isolated periodic bumpy-wall shear adaptation | Transient discrepancies unresolved; full benchmark not run |

Historical outcomes concern the recorded campaign, not a fresh acceptance
result for the current source. Successful compilation, a completed solve or a
readable frame is not automatically a physical-comparison PASS. Excluded models
and unavailable operators remain documented rather than replaced by tuned ones.

## Current-source scope

The retained campaigns predate subsequent contact and Bond changes. In the
current native contact implementation, the normal-force truncation line is
disabled; the resultant can remain signed while the contact overlaps. The
ON-DEM adapter owns a separate release convention and retains its non-attractive
default. Its contact equation must not be changed to match the native law.

Bond damage/energy changes and the removal of persistent endpoint-torque arrays
are not validated by the retained contact-comparison outcomes. These families
do not supply current Bond-fracture acceptance evidence. Use the current
constitutive tests for regression checks and a separately recorded comparison
campaign for new physical-validation claims. Do not relabel the stored results
as a run of the latest source.

## Build and run

Keep the terminal in the directory containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_VALIDATIONS=ON
cmake --build build --target lammpsComparison pfcComparison --parallel
```

These are build-only aggregate targets, not executables or implicit full-suite
runs. Each case defines its model in `cases/<id>/<id>.cpp` and links the public
core without Qt, a desktop project loader or a reference-program installation.
Use the family README for an executable command and
[shared options](common/README.md) for output and verdict semantics.
`FUNDEM_BUILD_VALIDATIONS=OFF` omits comparison executables.
ON-DEM's input retrieval and long run are opt-in; restricted coordinates are
not distributed with the source.

## Read and preserve evidence

Start with a family's README, then a case's `REPORT.md`. The English reports
use an abstract, benchmark definition, numerical method, results and discussion,
conclusions and limitations, references, and a reproducibility appendix.
[Publication specifications and PDF generation](PUBLICATION.md) describe the
numbered equations and tables, vector PDF/SVG figures and 600 dpi raster exports.
Figures are derived from retained CSV checkpoints; they do not replace exact
numerical records or demonstrate a new run.

Keep new results outside source in a new directory for each invocation.
`HISTORICAL-*` documents retain detailed prior reports, exact results and hashes;
they are not rewritten to match a later run. Full native archives and large
histories are retained separately. Never convert a historical FAIL into PASS
because a new input compiles or exports successfully.
