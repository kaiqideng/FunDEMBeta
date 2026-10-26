# `src/material/`

Bulk material properties and contact-pair coefficients are separate. The same
`material` can be assigned to spheres and level-set particles; particle geometry,
not a material marker, selects the contact family.

## Bulk properties

`material(youngModulus, poissonRatio, restitutionCoefficient, density)` stores:

| Property | Unit | Accepted values |
| --- | --- | --- |
| Young's modulus | Pa | Finite and positive |
| Poisson ratio | Dimensionless | Finite, strictly between -1 and 0.5 |
| Restitution coefficient | Dimensionless | Finite, within [0, 1] |
| Density | kg/m^3 | Finite and positive |

A default material is incomplete until its modulus and density are assigned.
Setters validate before changing data; `setProperties()` validates the entire
definition before replacing any value. `materialContainer` uses host AoS/device
SoA storage with exactly these four fields. There is no separate `LSMaterial`
class, material-kind marker, stiffness, or friction in this class.

## Contact-pair properties

`materialPair(firstMaterialIndex, secondMaterialIndex)` represents an unordered
pair, including two particles using the same material. It canonicalizes its
indices so that the smaller index is first.

| Setter | Parameters | Unit |
| --- | --- | --- |
| `setLevelSetStiffness()` | Normal and sliding stiffness per nodal area | N/m^3 |
| `setLinearStiffness()` | Normal, sliding, rolling, torsional stiffness | N/m |
| `setFrictionCoefficients()` | Sliding, rolling, torsional friction | Dimensionless |
| `setDampingCoefficients()` | Normal, sliding, rolling, torsional dashpots | N s/m |

Each pair has exactly two Boolean switches, both initially false:

- `linearStiffnessEnabled()`: selects all four ordinary contact stiffnesses
  together, only for sphere-sphere and sphere-level-set contacts. When false,
  Hertz normal and Mindlin sliding response are used, and the length-scaled
  rolling and torsional stiffnesses each equal that sliding stiffness.
- `dampingEnabled()`: selects all four explicit dashpots together. When false,
  damping is calculated from the material restitution coefficients. When true,
  the stored coefficients are used, including explicitly configured zeros.

Use `setLinearStiffnessEnabled(bool)` and `setDampingEnabled(bool)` to switch
groups without erasing their values. Group value setters and individual ordinary
stiffness or dashpot value setters automatically enable the corresponding whole
group; they do not create per-direction switches. Friction has no switch.

LS normal and sliding stiffnesses per area are independent of both switches and
must both be configured for LS-LS contacts. Their initial `-1` values mean
unset; `hasLevelSetStiffness()` is true when both are nonnegative. Explicit zero
is valid. All other coefficients initially equal zero, and public coefficient
setters accept only finite nonnegative values. The ordinary linear-stiffness
switch has no effect on LS-LS response. With explicit damping enabled, LS nodes
receive their area-weighted shares of the normal/sliding dashpots; rotational
dashpots are ignored for LS-LS contacts.

Rolling and torsional histories use length-scaled displacements, so their
coefficients are not direct angular stiffnesses or angular dashpots.
Using the default sliding stiffness for both rotational springs is FunDEM's
prescribed extension, not a classical Hertz-Mindlin rolling or torsional law.
There is no independent rotational-stiffness override while the linear group
is off; enabling it selects all four stored stiffnesses together.
`materialPairCoefficients` is a plain aggregate of 13 scalar coefficients and
the two Boolean switches, with no bitmask or directional activation flags.
`materialPair::coefficients()` creates this aggregate without pair identity.
Device storage contains 17 fields: these 15 values plus the two material IDs.

## Stable pair storage

`materialPairContainer(materialCount)` owns all `n * (n + 1) / 2` unordered pairs.
`at(i, j)` and `at(j, i)` access the same record.
`resizeForMaterials(newCount)` appends records when materials are added, retaining
existing overrides and pair indices; shrinking is rejected. The stable slot is:

```cpp
max(i, j) * (max(i, j) + 1) / 2 + min(i, j)
```

Host configuration uses checked `at()` access. Kernels use the host/device
`index()` function with validated material IDs. Uploads, capacity management, and
stream synchronization follow the standard `hostAoSDeviceSoA` interface.
Adding materials or changing coefficients requires refreshing the device table
before evaluating contacts. Pair values are effective interaction coefficients;
they are not combined a second time as separate particle-side stiffnesses.
