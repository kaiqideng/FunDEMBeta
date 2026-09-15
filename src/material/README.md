# `src/material/`

> **Purpose:** Define sphere and level-set contact materials with a unified `DeviceLayout`.

## Core Contents

| Type | Responsibility |
| --- | --- |
| `material` | Store normal, sliding, rolling, and torsional contact parameters together with density |
| `LSMaterial` | Reuse `material` storage while exposing level-set stiffness-per-area interfaces |
| `materialContainer` | Store materials in Host AoS/Device SoA form |

## Key Rules

- Stiffnesses and friction coefficients must be finite and nonnegative.
- Restitution coefficients must lie in `[0, 1]`.
- Density must be finite and strictly positive; a material becomes valid only after density is set.
- Infinite mass is represented by fixed level-set geometry or rigid-body mass state, not by zero density.
- `materialType::standard` is the default; every `LSMaterial` is permanently marked as `materialType::levelSet`.
- `LSMaterial` publicly exposes sliding, rolling, and torsional friction coefficient getters/setters.
- LS–LS node contacts scale `LSMaterial` stiffness by nodal area.
- Sphere–LS contacts use all four stiffnesses directly from the sphere material, while friction and restitution still combine both materials.

## Extension Guide

`LSMaterial` must add no data members and must preserve `sizeof(LSMaterial) == sizeof(material)`, allowing one `materialContainer` to serve both CPU and GPU paths.

## Parameter Meaning

| Property | Ordinary `material` | `LSMaterial` |
| --- | --- | --- |
| Normal stiffness | Force per overlap | Stiffness per nodal contact area |
| Sliding stiffness | Tangential force per displacement | Shear stiffness per nodal contact area |
| Rolling/torsional stiffness | Explicit stored values | Hidden and initialized to zero |
| Friction | Sliding, rolling, and torsional values | All three coefficients remain configurable |
| Restitution | Dimensionless value in `[0, 1]` | Same storage and validation |
| Density | Used to construct particle mass | Used with integrated LS volume |

The level-set subclass changes names and access restrictions, not layout. This allows a single `materialContainer` and one device view to serve all contact kernels while the `materialType` marker selects area-scaled LS contact behavior.

## Contact Combination Rules

- Sphere-sphere response uses the ordinary material properties of both spheres.
- Sphere-LS stiffness comes directly from the sphere material; friction and restitution combine both materials.
- LS-LS stiffness is area-scaled using the level-set material marker and the sampled surface-node area.
- `effectiveRadius` remains a physical lever arm and never selects the model.

Material insertion copies the value into solver-owned storage and binds its container identity. A particle must reference the material container owned by the solver that receives it. Do not copy a particle configured against one solver into another solver without assigning its material again.

## Validation and Diagnostics

The default-constructed base material is invalid until a positive density is set. Setters validate before changing storage, so a rejected value does not partially update the material. If `addMaterial()` or particle insertion fails, check type compatibility (`material` for spheres and `LSMaterial` for LS particles), density, coefficient finiteness, and restitution bounds.
