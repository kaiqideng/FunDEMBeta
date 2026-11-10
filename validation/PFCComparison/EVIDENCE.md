# Evidence integrity

The 30 October 2026 presentation cleanup changes reading order and figure style,
not the historical experiment, acceptance criteria or numerical evidence.

| Record | Meaning |
| --- | --- |
| `REPORT.md` | Concise current presentation of the historical comparison |
| `HISTORICAL-REPORT.md` | Byte-identical report before this presentation cleanup |
| `HISTORICAL-COVERAGE.md` | Byte-identical full coverage inventory before this cleanup |
| `HISTORICAL-RESULTS.csv` | Unchanged exact historical numerical fields |
| `checkpoints/*.csv` | Unchanged compact numerical inputs for the figures |
| `PROVENANCE.md` | Historical source/archive hashes; archived-report hash identifies `HISTORICAL-REPORT.md` |
| `EVIDENCE-SHA256SUMS.txt` | Current distributed documents, evidence and figures |
| `HISTORICAL-EVIDENCE-SHA256SUMS.txt` | Unchanged pre-cleanup manifest; paths describe that archived snapshot |

The old manifest is not a checksum list for rewritten current documents. For
example, its `cases/<id>/REPORT.md` digest now identifies that case's
`HISTORICAL-REPORT.md`, and its `COVERAGE.md` digest identifies
`HISTORICAL-COVERAGE.md`. Its original README digest refers to the separately
retained pre-cleanup snapshot. Do not replace those old digests with new values
or interpret the old manifest against a different directory tree.

For the current tree, run from `PFCComparison/`:

```bash
sha256sum -c EVIDENCE-SHA256SUMS.txt
```

`plotCheckpoints.py` reads the three retained checkpoint CSV files and writes
vector SVG/PDF figures and 600 dpi PNG exports. It checks that all expected
historical rows are present, that coefficient errors agree with their recorded
targets, and that input hashes remain unchanged. It does not interpolate a new
physical trajectory, run PFC or FunDEM, or produce a new PASS/FAIL verdict.

Detailed native archives, original JSON records and per-step histories remain
outside this source distribution. Their recorded hashes are identifiers, not
claims that those archives are supplied here. PFC itself was not executed
locally. The historical three PASS outcomes retain their original scope.
