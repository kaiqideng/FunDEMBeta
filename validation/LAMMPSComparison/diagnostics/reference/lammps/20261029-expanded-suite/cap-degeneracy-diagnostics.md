# Conditional-cap source diagnostics

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

These are source calculations, not executed validation cases. The original input, material Hertz normal law, initial pose and timestep are preserved. The pinned source explicitly defines wall effective E as E/[2(1-nu^2)] in gran_sub_mod.cpp, so the two native materials retain E_source without doubling.

## Hertz rolling: decisive first-step counterexample

The initial normal preload is assembled before the first drift. With zero initial rolling history, the first rolling displacement is `R*omega0*dt`. Native Hertz rolling stiffness equals its overlap-dependent Mindlin sliding stiffness; the source selects the independent constant `kr=10000 N/m`.

| Input | First overlap (m) | Source trial (N) | Native trial (N) | Common cap (N) | Cap state |
|---|---:|---:|---:|---:|---|
| dem07-hertz-rolling-3d-si.yaml | 6.06999999995e-05 | 2.5e-05 | 8.76932654122e-09 | 8.8451528572e-06 | Source capped; native elastic |
| dem07-rolling-mu02-3d-si.yaml | 6.06999999995e-05 | 2.5e-05 | 8.76932654122e-09 | 3.53806114288e-05 | Both elastic, unequal |

Both inputs access different rolling forces on their very first nonzero history update. Their stiffness mismatch cannot be hidden by common Coulomb saturation. They remain excluded from faithful native verification. This result does not assume that every coefficient mismatch is observable; it provides a concrete observable counterexample.

## Active tangential damping: conditional probes

The three candidates in the probe queue (retained locally; not included) are not promoted or counted as verified. An independent replay must advance the official Mindlin history on every actual native contact row, use the source `ct=0.5*cn`, verify the source/native normal-force friction bounds and slip signs, and compare the resulting total tangential force with the saved native force at every step. A final analytic PASS alone is insufficient.
