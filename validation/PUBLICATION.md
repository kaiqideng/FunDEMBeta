# Publication-format verification reports

This edition presents the retained comparison campaign as an English technical
report series. It is a reproducible publication format, not a claim of peer
review, experimental validation or fresh acceptance of the current solver.
The [family index](README.md) identifies all included cases and their outcomes.

## Report structure

Every LAMMPS/PFC case has its own editable `REPORT.md`, with a title and dated
outcome, Abstract, Benchmark definition, Numerical method and verification
procedure, Results and discussion, Conclusions and limitations, References,
and a reproducibility appendix. Tables, figures and equations are numbered
locally within each case. The original tables and their numerical precision
are preserved. References identify the actual authority: pinned LAMMPS inputs
and assertions, or PFC's public analytical targets. PFC was not run locally.

The bound PDF uses A4 pages, 180 mm text width, embedded STIX General type,
10.5 pt body text, 9 pt reproducibility appendices, numbered equations, ruled tables without vertical borders,
a linked contents list and case bookmarks. The ON-DEM adaptation is an
additional scope report; its unresolved transient discrepancies and missing
full benchmark run remain explicit. No unsupported ON-DEM result plot is added.

## Figure specification

Each of the 38 LAMMPS/PFC cases supplies the same figure in three formats:

| File | Intended use |
| --- | --- |
| `figures/comparison.pdf` | Vector figure with embedded TrueType fonts |
| `figures/comparison.svg` | Editable vector artwork |
| `figures/comparison.png` | 600 dpi raster fallback at 180 mm width |

Figures use 8.5 pt scientific type, panel letters, physical units, explicit
axis scale factors and distinguishable marker shapes and colors. Captions
define the selected body, observable, reference, backend, time step and error
normalization. Zero errors displayed on logarithmic axes use a stated drawing
floor only; source values and acceptance decisions are unchanged. Sparse saved
states are plotted as checkpoints without fitted or interpolated histories.
The two historical LAMMPS FAILs remain visible. PFC's absolute coefficient
criterion does not conceal the large relative errors near critical damping.

Each family's `FIGURES.json` records exact source CSV hashes, fields, selected
rows, panel definitions, captions, physical dimensions, artifact hashes and
layout checks. `HISTORICAL-*` files and the old checksum ledgers remain frozen.
Current presentation checksums are refreshed separately.

## Reproduce the presentation

Use Python 3.10 or later with `matplotlib`, `reportlab` and `pypdf` installed.
Run from the directory containing `FunDEMBeta/`. Choose an output directory
outside the core source tree. These commands do not run either solver.

```bash
python3 FunDEMBeta/validation/LAMMPSComparison/tools/plotCheckpoints.py --write
python3 FunDEMBeta/validation/PFCComparison/plotCheckpoints.py
python3 FunDEMBeta/validation/LAMMPSComparison/tools/plotCheckpoints.py --check
python3 FunDEMBeta/validation/PFCComparison/plotCheckpoints.py --check
python3 FunDEMBeta/validation/common/buildPublication.py --output publication-output/FunDEM_Verification_Report.pdf
python3 FunDEMBeta/validation/common/auditDocumentation.py
```

The PDF builder rejects an existing output path, embeds the original vector
figures, and records source hashes, page count and vector placements beside
the PDF. Its scratch assets stay outside the source tree. It does not modify
case reports or numerical records. A publication delivery also requires visual
inspection of rendered PDF pages and figures; an integrity audit alone cannot
establish correct scientific content or readable layout.

For a portable distribution, place a copy of this `validation/` directory
beside the bound PDF. Its relative evidence hyperlinks then resolve without
the original workstation path. Running native C++ inputs still requires the
FunDEM core checkout and the build instructions in the family guides.

After an intentional report or figure change, update each affected report's
current hash in `PROVENANCE.md`, then explicitly refresh the current family
ledgers with `auditDocumentation.py --write-manifests`. Do not refresh archived
report hashes or historical CSV hashes. The shared report builder, drawing
style and documentation audit are in [common](common/README.md).
