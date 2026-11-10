# Evidence integrity

This presentation cleanup does not rerun the 2026-10-28–29 campaign, change
acceptance criteria, or replace its 33 PASS and two FAIL outcomes.

| Record | Meaning |
| --- | --- |
| `REPORT.md` | Concise presentation of the historical comparison |
| `HISTORICAL-REPORT.md` | Byte-identical report present immediately before this cleanup |
| `HISTORICAL-RESULTS.csv`, `checkpoints/*.csv` | Unchanged historical numerical records |
| `PROVENANCE.md` | Source pins and archive hashes, including the corresponding `HISTORICAL-REPORT.md` digest |
| `FIGURES.json` | Exact checkpoint and generated image hashes, selected columns, and figure method |
| `EVIDENCE-SHA256SUMS.txt` | Current distributed files and figures |
| `HISTORICAL-EVIDENCE-SHA256SUMS.txt` | Unchanged earlier checksum ledger, with its original path labels |

## Historical ledger limitation

The old ledger's 89 CSV checksums still match the retained numerical tables.
Its 35 `REPORT.md` checksums were already out of date before this cleanup:
they do not identify the report bytes present at the start of this task.
Those old digests are preserved, not replaced or presented as successful
checks of newer documents.

Each archived report matches both the pre-cleanup backup and the report digest
previously recorded in its `PROVENANCE.md`. That field is now labeled
`Archived report SHA256 (HISTORICAL-REPORT.md)`. This is distinct from the
older campaign report hashes also retained in the provenance record.

For the current distributed tree, run from `LAMMPSComparison/`:

```bash
sha256sum -c EVIDENCE-SHA256SUMS.txt
```

The [shared audit](../common/README.md#documentation-checks) independently
checks historical CSVs, the archived reports, current files and local links.
It lists the old ledger/report discrepancies explicitly, rather than altering
historical evidence to manufacture a clean history.

## Figures

The drawing script reads the retained official-metric or checkpoint CSV for
each case. It recomputes the stored metric error where defined, checks recorded
criterion flags, and keeps physical units on checkpoint axes. Exact zero errors
on log axes are marked at the stated plotting bound, not replaced in the data.
Plots contain recorded markers and reference guides, with no smoothing,
trajectory interpolation, parameter fitting or new PASS/FAIL assignment.

Original native archives, full trajectories and vendor inputs remain outside
this source distribution. Their recorded hashes do not imply that the files
are included here. LAMMPS was not rerun locally for most comparisons; each
report retains its own reference-execution and model-mapping restrictions.
