# LAMMPS Hertz oblique: proved saturated-trajectory force equivalence

> Historical source audit / diagnostic (28–29 October 2026), not current-core acceptance. Capability statements apply only to the recorded snapshot. Numerical evidence and outcomes are unchanged.

The official `dem04-hertz-oblique` input selects a direct Hertz normal coefficient, an independent constant tangential stiffness `kt=1.5e10 N/m`, tangential damping zero, and friction0.092. The native Hertz interface does not independently select that tangential stiffness. A faithful **observable-force comparison for this specific supplied trajectory** is nevertheless possible: the original source history is independently replayed, and both total tangential forces stay on the same Coulomb cap in every contact row.

This is not a statement that the complete constitutive models or their internal histories are identical. It does not validate unsaturated stiffness/damping or tangential spring energy. Source parameters are preserved in the replay; no rebound, friction or force is fitted. The source is tagged `granular` and is **not tagged `unvetted`**.

| Native run | Contact rows replayed | Maximum Ft difference (N) | Maximum Fn law difference (N) | Unsaturated source/native rows |
|---|---:|---:|---:|---|
| cpu | 1397 | 5.3290705182e-15 | 5.68434188608e-14 | 0 / 0 |
| gpu | 1397 | 5.3290705182e-15 | 5.68434188608e-14 | 0 / 0 |
| cpu-dt-half | 2793 | 7.1054273576e-15 | 5.68434188608e-14 | 0 / 0 |
| gpu-dt-half | 2793 | 5.3290705182e-15 | 5.68434188608e-14 | 0 / 0 |

All four modes proved: **True**. The complete per-step force/history components, signs, cap values and exact replay-source snapshots are linked in the JSON evidence (retained locally; not included).

The wall torque convention remains different: source R and native R-overlap. Its measured arm/slip difference, actual native torque assembly and angular-impulse balance are independently reported. Passing the official rebound/spin tolerance does not establish equality of this lever arm or eliminate the resulting geometric bias.
