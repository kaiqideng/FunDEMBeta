# LAMMPS diagnostic: dem06-spin-wall

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 29 October 2026. **CPU-only mapping probe; excluded from verified comparisons.** This is not a four-run benchmark and is not an official analytic FAIL.

The pinned [official input](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem06-spin-wall-3d-si.yaml) uses commit 29dd3fd7642a6fe9e8b4dd99120df4068b221250, Git blob 336074ed3cc40fff1a647ef9961f0e648d920125. LAMMPS was not rerun locally. Its Mindlin history and source ct=0.5cn were replayed on every saved native contact state. Native normal forces agree within 1.8189894e-12 N, but source-history tangential-force equivalence **FAILS**: not every step shares the same total Coulomb saturation.

The first counterexample is step 51 at 5.1e-05 s. Source trial magnitude=736491.0838 N; cap=353.329957106 N; source capped=True. Source/native force residual=706.659914212 N, versus declared replay tolerance 9.51370025157e-08 N. Source/native ct=2204.61457892/3949.63383181 N s/m. Maximum history replay residual=1810.11226529 N.

This contact-law mismatch prevents a faithful unchanged-parameter reproduction. No coefficients were tuned and no GPU/half-step runs were claimed. Source undeformed-radius slip and native shared-contact-point slip are also distinguished in the proof.

One actual CPU execution retains 1,659 native frames, original project snapshots, every-step trajectory/contact CSV and logs. Complete native result (retained locally; not included) | Source-history proof (retained locally; not included) | All-contact replay CSV (retained locally; not included) | Retention audit (retained locally; not included) | Full excluded ZIP (retained locally; not included) | [Combined English diagnostic report](../LAMMPS-DIAGNOSTICS.md).
