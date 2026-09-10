# `execution/`

> **Purpose:** Provide stateless numerical formulas shared by CPU and GPU code without owning containers or solver flow.

## Core Contents

| File | Responsibility |
| --- | --- |
| `motionIntegration.h` | Integrate velocity, angular velocity, position, and unit quaternion |
| `contactDetection.h` | Compute effective mass and shared sphere/level-set contact geometry |
| `contactFunctions.h` | Evaluate normal, sliding, rolling, and torsional response with Coulomb limits |
| `bondFunctions.h` | Evaluate bond force/torque, stiffness, and BK mixed-mode fracture energy |
| `sphFunctions.h` | Provide Wendland kernels, reference-impedance acoustic Riemann states, density/pressure updates, and fluid/wall interaction formulas |
| `sphJetFunctions.h` | Test whether a particle lies inside an upstream SPH jet pipe |
| `sphSourceDiscretization.h` | Construct the regular cylindrical-source lattice |
| `cuda/` | Launch CUDA kernels that call these formulas |

## Key Rules

- Place mutable output parameters before read-only inputs.
- Shared functions use the project's host/device inline macro and contain no allocation, exceptions, or global state.
- Use `math::defaultTolerance` and helpers such as `math::isFinite()` for degeneracy checks.
- `inverseMass == 0` means infinite mass. Motion integration leaves its velocity, position, and orientation unchanged, with inertia and inverse inertia defaulting to zero matrices.
- Renormalize quaternions after integration to prevent long-term orientation drift.

Tangential damping is limited by projection onto the Coulomb disk: compute the trial tangential force, then constrain the combined elastic and damping force to `mu * |normalForce|`. Rolling and torsional components use their corresponding radius, stiffness, and friction limits.

## Extension Guide

- Implement a formula needed by both CPU and CUDA only once in this directory, then call it from object methods or kernels.
- Public headers must compile with both an ordinary C++ compiler and NVCC.
- Do not introduce kernel launches, device allocation, or container lifecycle operations into shared numerical functions.

## Layer Boundary

The `execution` layer knows mathematical values but not solver modes, containers, allocation, histories, or launch configuration. A typical call chain is:

```text
solver stage
    -> CPU container loop or CUDA kernel
        -> execution formula
            -> math primitives
```

This boundary keeps the physical equation identical on CPU and GPU. Object methods may organize cached inputs and outputs, but the equation itself belongs here when both backends need it.

## Contact Response Sequence

1. Contact detection supplies point, unit normal, overlap, area, effective mass, and effective radius.
2. Material coefficients produce normal, sliding, rolling, and torsional stiffness/damping values.
3. Normal response is evaluated first and records `normalForceMagnitude`.
4. Sliding damping is projected onto the Coulomb disk together with the elastic spring force.
5. Rolling and torsional responses use their own spring state, stiffness, friction coefficient, and lever arm.
6. Force, torque, elastic energies, and updated history values return to the caller.

Degenerate contacts return zero response rather than creating non-finite values. Normals should already be normalized by contact construction; guarded normalization remains appropriate when data crosses a public boundary.

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

CPU loops and CUDA kernels compose these functions in the same order. The physical support remains `2h`, even when the neighbor search radius is enlarged with a displacement buffer so one neighborhood can be reused across an advection step.

Finite-duration jets reuse one host/device pipe-membership predicate from `sphJetFunctions.h`. The solver performs the one-time upstream particle allocation and owns each source's stable range and end time; the execution formula owns no source state. For every active acoustic substep, the solver applies the prescribed velocity before optional advection preparation and again after velocity integration. Once the final end time is reached, the solver's `jetsCompleted` fast path stops calling the range checks entirely. The legacy four-argument `addSPHJet` overload creates only a static downstream column and does not enter this sequence.

## Review Checklist for a New Formula

- Inputs are finite or safely guarded.
- Mutable output references appear before constant inputs.
- No container, allocation, exception, or global mutable state is introduced.
- CPU and device compilation take equivalent branches.
- Units and sign conventions match the calling object.
- A focused invariant test covers zero, symmetric, and limiting cases.
