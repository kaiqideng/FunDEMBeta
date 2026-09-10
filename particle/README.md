# `particle/`

> **Purpose:** Define moving objects, particle containers, and the multilevel spatial grid without owning contact search or solver control flow.

## Core Contents

```text
pointMass
├── rigidBody
│   └── particle
│       └── LSParticle
├── SPHParticle
└── virtualParticle
```

| Type | Responsibility |
| --- | --- |
| `pointMass` | Position, velocity, force, and inverse mass |
| `rigidBody` | Orientation, angular velocity, torque, and inertia tensors |
| `particle` | Radius, material index, and DEM state for an ordinary sphere |
| `LSParticle` | Level-set geometry reference, grid, surface, and mass properties |
| `SPHParticle` | Density, pressure, density rate, and free-surface state |
| `virtualParticle` | Actual volume, normal, and interpolated state of an SPH wall sample |
| `spatialGrid` | Particle background grid with at most four levels |

## Key Rules

### Containers and Device Data

- User-insertable objects use `isValid()` and container ownership to prevent invalid or duplicate insertion; internal helper objects do not.
- Every particle container uses Host AoS/Device SoA storage, and only fields declared in `DeviceLayout` enter the GPU.
- Host references to material/geometry containers and CPU calculation copies never become device pointers.
- Uniform SPH properties are owned by `SPHDEM`, not repeated in every `SPHParticle`.
- `virtualParticle` objects are generated internally by the solver and cannot be inserted directly by users.

### Level-Set Geometry

Each host-side `LSParticle` copies the referenced geometry range to reduce indirection during CPU search. GPU kernels retain only the geometry index and read a separate `LSGeometryContainer` device view. Geometry or particle changes therefore require device data to be reinitialized by the solver.

Fixed level-set geometry has zero volume and zero unit-density inertia. An `LSParticle` using it keeps zero inverse mass and zero inertia tensors, representing an infinite-mass rigid body.

### Spatial Grid

The spatial grid consists of two storage groups:

- `spatialGridParticle`: stores `hash_` and the original particle `index_` for each sorted slot; the device fields are `particleHash_` and `particleIndex_`.
- `spatialGridCell`: stores `sortedParticleIndexBegin_` and `sortedParticleIndexEnd_` for each cell.
- `spatialGrid`: stores boundaries, cell size, dimensions, and the two containers above.

`particleIndex_[sortedIndex]` maps a sorted slot back to the original particle container. The particle container itself is never reordered, so contact and history indices remain stable after rebuilding the grid.

The multilevel grid supports at most four levels. Finite-mass particles are assigned to the earlier levels by size, while infinite-mass particles occupy a dedicated final level so they become the slave in cross-level pairs.

## Extension Guide

- Add a new state field to `DeviceLayout` only when GPU calculations require it.
- Put formulas shared by CPU and GPU in `execution/` instead of duplicating them in particle classes.
- Add ownership and `isValid()` only to types that users can insert directly; keep internal temporary types simple.
- Use `math::defaultTolerance` for degeneracy checks instead of hard-coded tiny values.

## State Inheritance and Device Layout

Each derived type concatenates the base `DeviceLayout` with only its additional fields:

| Layer | Added dynamic state |
| --- | --- |
| `pointMass` | Position, velocity, force, inverse mass |
| `rigidBody` | Orientation, angular velocity, torque, inertia, inverse inertia |
| `particle` | Radius and material index |
| `LSParticle` | Geometry index and LS grid/surface descriptors needed by device contact detection |
| `SPHParticle` | Density, pressure, density rate, free-surface flag, viscous prior force |
| `virtualParticle` | Local/world boundary state, volume, owner index, acceleration, prior force |

The host object may contain additional references or vectors. Their absence from `DeviceLayout` is intentional. Kernel code receives lightweight device views assembled from the active field arrays.

## Mass and Motion Semantics

`inverseMass == 0` always means infinite mass. `pointMass::setMass()` converts a positive finite mass into its inverse. `rigidBody::setMassAndInertia()` accepts the mass only when the inertia tensor is invertible; otherwise it resets mass, inertia, and inverse inertia to the infinite-mass state.

The four shared integration operations are:

1. Linear velocity from force, inverse mass, gravity, and time step.
2. Angular velocity from torque and current orientation/inertia.
3. Position from velocity.
4. Unit-quaternion orientation from angular velocity.

Infinite-mass rigid bodies skip both translational and rotational integration. Prescribed motion should therefore be assigned explicitly between completed solve stages rather than expected to arise from force.

## Level-Set Particle Representation

An `LSParticle` has two intentionally different geometry access paths:

- Host CPU search reads the copied grid and surface ranges stored with the particle for contiguous access.
- CUDA kernels read the shared `LSGeometryContainer` through `geometryIndex_` to avoid repeating geometry arrays on the device.

`setGeometry()` obtains volume, bounding radius, and unit-density inertia from the descriptor and combines them with material density. Fixed geometry leaves these values zero and produces infinite mass. Surface-node indices used by mapping and LS-LS contact refer to the central geometry surface-node container, not a particle-local expanded list.

## SPH and Virtual Particles

SPH material properties are uniform solver configuration. A user-created `SPHParticle` supplies position and velocity; insertion assigns mass and reference density. CPU and GPU interaction stages update density, pressure, density rate, force, and free-surface state.

Virtual particles are sampled on the uniform SPH lattice anchored by the first fluid particle. The LS grid is interpolated at each candidate point, and solid-side samples within a two-smoothing-length boundary band are retained with volume equal to one SPH lattice cell. This common lattice prevents boundary samples from overlapping fluid particles by half a spacing. Local position and normal follow the LS owner; world position, velocity, acceleration, force, and normal are refreshed by the coupling layer. Users never append virtual particles directly.

## Spatial-Grid Ownership

Grid sorting reorders only `(hash, particleIndex)` records. It never reorders the particle container. Cell ranges refer to sorted slots, and `particleIndex_[sortedIndex]` maps back to stable user insertion indices. When adding a grid field, keep this distinction explicit in both host names and device views.

For rigid multilevel grids, size bands reduce unnecessary cell traversal for large diameter ratios. Infinite-mass particles occupy the final level to enforce slave ordering in cross-level pairs. SPH grids are uniform because smoothing length is shared by the whole fluid.

## Common Failure Modes

- A sphere is invalid: set a positive radius and assign an ordinary material from the destination solver.
- An LS particle is invalid: assign an `LSMaterial` and a valid geometry from the destination solver.
- Mass remains infinite: check geometry volume, material density, and inertia invertibility.
- Device state appears stale: reinitialize or perform the appropriate host/device copy after changing host configuration.
- Contact indices change unexpectedly: do not reorder or erase solver-owned host objects.
