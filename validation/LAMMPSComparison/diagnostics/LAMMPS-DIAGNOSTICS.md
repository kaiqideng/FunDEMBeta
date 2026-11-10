# LAMMPS diagnostic exclusions

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 29 October 2026. Five completed excluded attempts retain 22,325 frames; startup preflight contributes zero executions/frames. These are not official analytic FAIL outcomes or four-run validations.

- [dem04-oblique-grazing](dem04-oblique-grazing/DIAGNOSTIC-REPORT.md): one CPU-only source-history mapping counterexample, 4,863 native frames.
- [dem06-spin-wall](dem06-spin-wall/DIAGNOSTIC-REPORT.md): one CPU-only source-history mapping counterexample, 1,659 native frames.
- [dem06-mindlin-spin-heavy](dem06-mindlin-spin-heavy/DIAGNOSTIC-REPORT.md): one CPU-only source-history mapping counterexample, 1,801 native frames.

The wrong E-times2 Tsuji wall mapping contributes 2 superseded base runs / 14,002 frames; corrected canonical runs are separate. Both source and native E_eff equal E/[2(1-nu^2)]. The cylinder padding rejection occurred before any integration.

The full excluded-attempt ledger and modulus-correction record are retained in
the original campaign archive, not distributed here.

## Source audits and coverage

These documents explain model mapping; they do not add accepted cases or
change either historical official FAIL outcome.

- [Full pinned 61-input coverage](reference/lammps/20261029-expanded-suite/coverage-matrix.md)
- [Conditional cap counterexamples](reference/lammps/20261029-expanded-suite/cap-degeneracy-diagnostics.md)
- [Active-history exclusions](reference/lammps/20261029-expanded-suite/conditional-history-probe-results.md)
- [Conditional Hertz-oblique equivalence](reference/lammps/20261029-expanded-suite/hertz-oblique-saturation-equivalence.md)
- [Family reports and reproduction](../README.md)
