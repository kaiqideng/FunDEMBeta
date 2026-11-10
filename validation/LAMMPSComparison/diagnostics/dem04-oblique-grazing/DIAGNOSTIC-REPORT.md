# LAMMPS diagnostic: dem04-oblique-grazing

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

Date: 29 October 2026. **CPU-only mapping probe; excluded from verified comparisons.** This is not a four-run benchmark and is not an official analytic FAIL.

The pinned [official input](https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem04-oblique-grazing-3d-si.yaml) uses commit 29dd3fd7642a6fe9e8b4dd99120df4068b221250, Git blob b818942ecd2238aaccd3a9c42423fcb9c306ba6a. LAMMPS was not rerun locally. Its Mindlin history and source ct=0.5cn were replayed on every saved native contact state. Native normal forces agree within 2.8421709e-14 N, but source-history tangential-force equivalence **FAILS**: not every step shares the same total Coulomb saturation.

The first counterexample is step 2237 at 2.237e-05 s. Source trial magnitude=0.00453817992281 N; cap=0.0302570775012 N; source capped=False. Source/native force residual=0.0257188975784 N, versus declared replay tolerance 1e-10 N. Source/native ct=0.0689265761495/0.12348405008 N s/m. Maximum history replay residual=0.034005393656 N.

This contact-law mismatch prevents a faithful unchanged-parameter reproduction. No coefficients were tuned and no GPU/half-step runs were claimed. Source undeformed-radius slip and native shared-contact-point slip are also distinguished in the proof.

One actual CPU execution retains 4,863 native frames, original project snapshots, every-step trajectory/contact CSV and logs. Complete native result (retained locally; not included) | Source-history proof (retained locally; not included) | All-contact replay CSV (retained locally; not included) | Retention audit (retained locally; not included) | Full excluded ZIP (retained locally; not included) | [Combined English diagnostic report](../LAMMPS-DIAGNOSTICS.md).
