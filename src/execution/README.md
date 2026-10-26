# `src/execution/`

> **Purpose:** Provide shared numerical formulas and backend execution helpers, with solver flow owned by `solver/`.

## Core Contents

| File | Responsibility |
| --- | --- |
| `motionIntegration.h` | Integrate velocity, angular velocity, position, and unit quaternion |
| `contactDetection.h` | Compute effective mass and shared sphere/level-set contact geometry |
| `contactFunctions.h` | Evaluate normal, sliding, rolling, and torsional response with Coulomb limits |
| `contactMaterialFunctions.h` | Resolve material-pair overrides, Hertz-Mindlin defaults, nodal damping and normal elastic energy identically on CPU/CUDA |
| `bondFunctions.h` | Evaluate bond force/torque, stiffness, and BK mixed-mode fracture energy |
| `sphFunctions.h` | Provide Wendland kernels, reference-impedance acoustic Riemann states, density/pressure updates, and fluid/wall interaction formulas |
| `sphCouplingFunctions.h` | Sample prescribed boundary velocity differences without duplicating centripetal acceleration |
| `sphJetFunctions.h` | Test whether a particle lies inside an upstream SPH jet pipe |
| `sphSourceDiscretization.h` | Construct the regular cylindrical-source lattice |
| `cpu/` | Apply generic particle integration and interaction-force assembly over host containers |
| `cuda/` | Launch CUDA kernels that call these formulas |

## Key Rules

- Place mutable output parameters before read-only inputs.
- Shared numerical formulas use the project's host/device inline macro and contain no allocation, exceptions, or global state. CPU container loops and CUDA launchers are separate backend helpers.
- Use `math::defaultTolerance` and helpers such as `math::isFinite()` for degeneracy checks.
- `inverseMass == 0` means infinite mass. Motion integration leaves its velocity, position, and orientation unchanged, with inertia and inverse inertia defaulting to zero matrices.
- Renormalize quaternions after integration to prevent long-term orientation drift.

Tangential damping is limited by projection onto the Coulomb disk: compute the trial tangential force, then constrain the combined elastic and damping force to `mu * |normalForce|`. Rolling and torsional components use their corresponding radius, stiffness, and friction limits.

## Extension Guide

- Implement a formula needed by both CPU and CUDA only once in this directory, then call it from object methods or kernels.
- Public headers must compile with both an ordinary C++ compiler and NVCC.
- Do not introduce kernel launches, device allocation, or container lifecycle operations into shared numerical functions.

## Layer Boundary

The shared numerical formulas at the root of `execution/` consume mathematical values and do not own particle containers, histories, solver modes, or launch configuration. [`cpu/`](cpu/README.md) contains host container loops and temporary assembly storage; [`cuda/`](cuda/README.md) contains kernel launchers. Stateful SPH neighborhoods and boundary-coupling accumulators belong to `interaction/`. A typical call chain is:

```text
solver stage
    -> CPU container loop, interaction component, or CUDA kernel
        -> execution formula
            -> math primitives
```

This boundary keeps the physical equation identical on CPU and GPU. Object methods may organize cached inputs and outputs, but shared equations belong in the root numerical headers. Backend helpers perform the stage requested by the solver; step ordering, simulation time, initialization, and snapshots remain in `solver/`.

## Contact Response Sequence

1. Contact detection supplies point, unit normal, overlap, area, effective mass, and effective radius.
2. Material IDs select one symmetric pair record. Two Boolean switches select the ordinary four-stiffness group and the four-dashpot group; both backends call `resolveContactResponse()`.
3. For positive overlap, normal response retains the signed `k_n*overlap - c_n*v_n` value; its former zero-force truncation is commented out in the shared host/device function. `normalForceMagnitude` records its absolute value for friction, rotational resistance, and output. Nonpositive overlap still returns zero force and clears sliding history.
4. Sliding damping is projected onto the Coulomb disk together with the elastic spring force.
5. Rolling and torsional responses use their own spring state, stiffness, friction coefficient, and lever arm.
6. Force, torque, elastic energies, and updated history values return to the caller.

Degenerate contacts return zero response rather than creating non-finite values. Normals should already be normalized by contact construction; guarded normalization remains appropriate when data crosses a public boundary.

The bulk material stores only Young's modulus, Poisson ratio, restitution and
density. Hertz normal force uses both elastic compliances and `sqrt(R*overlap)`;
Mindlin sliding stiffness uses both shear compliances. The normal dashpot uses
the Hertz **tangent** stiffness, while the elastic normal potential is
`(2/5)*Fn_elastic*overlap`, not the linear `(1/2)*k*overlap^2` expression.
These definitions follow the material-based Hertz/Mindlin and coefficient-of-restitution
forms in the [LAMMPS granular model documentation](https://docs.lammps.org/pair_granular.html).
Tangential history remains the existing displacement-history Mindlin no-slip
approximation, not a full Mindlin-Deresiewicz partial-slip model.

LS nodal detection already stores `m_node = m_pair*A_node/A_patch`.
`nodalDampingWeight()` recovers the same area fraction without a second reduction
or another device field; explicitly supplied pair damping is multiplied by that
fraction. Restitution-derived nodal damping continues using nodal mass and nodal
stiffness. LS-LS always requires both per-area stiffness values; it does not use
the ordinary linear-stiffness switch, and ignores rotational dashpots. For sphere
contacts, disabling the linear-stiffness group selects Hertz-Mindlin, with both
length-scaled rotational stiffnesses equal to the default sliding stiffness.
This equality is FunDEM's prescribed rotational extension, not a classical
Hertz-Mindlin rolling law. Enabling the damping group applies all four
stored dashpots, including zero values. Dashpots can act without their associated
elastic spring, subject to the existing Coulomb limit.

## Bond Response Sequence

Bond formulas consume the current master/slave rigid-body copies, reference endpoint bases, equivalent length, four stiffnesses, and irreversible damage state. They calculate one force plus separate master/slave torques and four elastic energy components. BK damage is monotonic: once the damage factor reaches complete failure, later calls cannot heal the bond.

## WCSPH Sequence

The shared SPH formulas are split intentionally:

- Wendland kernel value and gradient;
- free-surface position-divergence contribution;
- density reinitialization;
- fluid and wall Riemann density-rate terms;
- pressure acceleration;
- viscous acceleration;
- wall reaction from each fluid acceleration contribution;
- acoustic/advection time-step limits.

CPU loops and CUDA kernels compose these functions in the same order. The physical support remains `2h`, even when the neighbor search radius is enlarged with a displacement buffer. `interaction/SPHNeighborhood` tracks actual displacement so the solver can reuse the search structures across acoustic substeps while the buffer remains valid, rebuilding them earlier when needed.

Finite-duration jets reuse the host/device pipe-membership predicate from `sphJetFunctions.h`. The solver owns a bounded upstream buffer, scheduled layer births, stable emitted-particle indices, and source completion; the execution formula owns no source clock. New layers are appended only at common acoustic starts, with timing quantized upward by at most one DEM step and insertion positions adjusted for that delay. Device-tail initialization preserves existing fluid fields and runs on the solver stream. The inlet imposes velocity and reference density before optional advection preparation and again after velocity integration. Outlet release advances a permanent released-layer prefix, so returning particles are not constrained again. Once all sources end, the completion fast path stops source work. The four-argument `addSPHJet` overload creates only a static downstream column and does not enter this sequence.

`isInsideSPHJetOpening` defines a separate virtual-wall mask: the finite upstream buffer cylinder expanded by half a particle spacing. Samples there use zero volume and cleared wall loads; samples that move out regain nominal volume. This is not an LS-grid or triangle edit. Changed quadrature requires a fresh viscous-prior-force preparation after the old owner hydraulic load has been removed. The same predicate drives CPU and CUDA, while observation applies it only to a reversible host view.

## Review Checklist for a New Formula

- Inputs are finite or safely guarded.
- Mutable output references appear before constant inputs.
- No container, allocation, exception, or global mutable state is introduced.
- CPU and device compilation take equivalent branches.
- Units and sign conventions match the calling object.
- A focused invariant test covers zero, symmetric, and limiting cases.
