# `src/interaction/`

> **Purpose:** Store contacts, bonds, neighbors, and coupling histories, and provide CPU searches and SPH interaction stages.

## Core Contents

| File | Responsibility |
| --- | --- |
| `contact.*` | Contact geometry, effective parameters, material coefficients, particle copies, force, and torque |
| `bond.*` | Bond endpoints, stiffnesses, BK damage state, force, and master/slave torques |
| `neighborRange.h` | Neighbor counts, inclusive prefix sums, and range history |
| `surfaceNodeMapping.h` | Map expanded surface nodes to geometry nodes and owning LS particles |
| `interactionContainer.*` | Manage container lifecycle, capacity, copies, and prefix sums |
| `contactSearch.*` | Construct sphere–sphere, sphere–LS, and LS–LS contacts on the CPU |
| `sphInteraction.*` | Build compact CPU SPH neighbor lists and evaluate WCSPH stages, integration, and reductions |
| `SPHNeighborhood.*` | Track reference positions and actual displacement while CPU or CUDA SPH searches are reused |
| `virtualParticleCoupling.*` | Accumulate boundary kinematics and resultant SPH forces and torques by LS owner |

## Key Rules

### Contact Conventions

- Sphere–sphere: order master/slave by mass rules; if only one particle has infinite mass, it is always the slave.
- Sphere–LS: the `LSParticle` is always the slave, and the complete spherical model includes sliding, rolling, and torsional torque.
- LS–LS: use one-way surface-node sampling. The smaller or finite-mass side becomes the master, and its surface nodes query the slave signed-distance field. Nodal area weights effective mass and force.
- A level-set master material selects the nodal area-scaled model. Standard master material selects the complete spherical model for both sphere–sphere and sphere–LS contacts.
- Sphere–LS normal, sliding, rolling, and torsional stiffnesses come directly from the sphere material; the LS material still contributes all three friction coefficients and restitution.
- `effectiveRadius_` is only a physical radius/lever arm; it no longer selects the contact model.
- `geometrySurfaceNodeIndex_` is a global index into the geometry surface-node container; `ownerParticleIndex_` indexes the LS-particle container.

CPU search may parallelize candidate generation and per-contact calculations with OpenMP, but accumulation into shared particles must remain race-free. CPU and GPU sphere–LS paths use identical particle ordering and torque conventions.

### Bond Setup

Configure a bond in this order:

1. Call `setEquivalentLength()`.
2. Set connection geometry from two particles or from a contact.
3. Call `setCoefficientB()` with the four stiffnesses.
4. Insert the bond into the solver container matching its particle types.

Changing equivalent length invalidates dependent connection geometry and stiffnesses, which must then be configured again. Each calculation cycle refreshes the master/slave particle copies and clears the previous force and torques before evaluating the bond response.

### Container Lifecycle

`interactionContainer` sizes both neighbor-range arrays and histories from the particle count. Growing the particle container only expands the required arrays and preserves existing contacts and contact history. Shrinking it resets all interaction state to prevent stale indices. The final prefix-sum value gives the active contact count, and contact storage grows only when needed.

After a GPU solve, only valid contacts and bonds are copied to the host. Neighbor ranges, temporary lists, and neighbor history are rebuildable data. Contact history survives normal repeated initialization, while neighbor history is reset each time.

## Extension Guide

- Configure contact geometry and material coefficients separately, without leaving a partially valid contact after failure.
- Keep host particle copies in contact and bond objects out of `DeviceLayout`.
- Put new shared CPU/GPU formulas in `execution/`; object methods should only organize inputs and retain results.
- For every new contact type, explicitly define master/slave ordering, history keys, and torque distribution.

## Interaction Data Model

| Type | Persistent role | Device role |
| --- | --- | --- |
| `contact` | Current geometry, coefficients, history state, loads, and energies | One active detected contact |
| `contactHistory` | Slave key and three spring deformations | Restores history after contacts are rebuilt |
| `bond` | Fixed particle pairing, endpoint bases, stiffness, damage, and loads | Updated every force stage |
| `neighborRange` | Candidate count and inclusive prefix sum | Allocates/writes compact result ranges |
| `particleNeighbor` | Candidate particle index and LS patch area | Intermediate rigid-neighbor list |
| `neighborRangeHistory` | Previous prefix boundary | Locates prior contacts for one master/range |
| `surfaceNodeMapping` | Geometry surface-node index and owning LS particle | Expands LS ranges without particle offsets |

Host-only master/slave particle copies let `contact` and `bond` call the same calculation methods on CPU. They are deliberately excluded from `DeviceLayout`; CUDA kernels read particle containers directly.

## CPU Search Flow

```text
read positions and radii
    -> create candidate particle pairs
    -> enforce master/slave convention
    -> sphere-sphere / sphere-LS / LS-LS geometry construction
    -> append valid contacts
    -> restore spring history
```

Sphere-LS always stores the sphere as master. LS-LS samples master surface nodes in the slave signed-distance field. Infinite-mass ordering is enforced before contact construction so a single fixed member becomes slave. Candidate generation may be parallel, but the final contact order and history key must remain deterministic enough for restoration.

## GPU Storage Lifecycle

1. `initializeDevice()` sizes range arrays from particles and, for LS particles, generates surface-node mappings.
2. Existing bonds are copied to device storage.
3. Existing host contacts are converted into compact contact/range history.
4. Detection kernels write neighbor counts.
5. Prefix sums produce total candidate or contact counts.
6. Contact or neighbor storage grows only when the required count exceeds capacity.
7. Detection writes active contacts and force kernels consume them.
8. `saveCurrentStepToHistory()` stores the current keys and spring states for the next rebuild.

Growing particle counts preserves current contacts and histories while expanding required ranges. Shrinking particle or surface-node counts resets interaction storage because retained indices could be invalid. Neighbor-range history is rebuildable and resets during normal initialization; contact history is preserved when possible.

## Force and Torque Assembly

Every contact produces one force and torque record before loads are accumulated into particles. For sphere-LS and LS-LS contact, torque signs follow the stored master/slave convention. CPU assembly uses race-free reductions or staged accumulation; CUDA uses atomics only where multiple interactions write one particle.

Bond calculations first refresh the two rigid-body copies and clear the previous force/torques. A bond is calculable only after equivalent length, connection geometry, and all four stiffnesses have been accepted. Damage state persists across calls and across normal device reinitialization.

## SPH Neighborhoods and Coupling

`sphInteraction.h` is private implementation and is not installed as a user API. Its CPU engine owns two compact host neighbor lists: fluid candidates per fluid particle and virtual-wall candidates per fluid particle. Both are constructed from sorted uniform-grid hashes at advection preparation and rebuilt earlier when displacement invalidates reuse.

`SPHNeighborhood` stores the fluid and boundary positions at the most recent search build separately from physical particle state. Before a density or pressure stage, `SPHDEM` checks that the maximum actual displacement is no greater than half the search skin. If that bound is exceeded or references are invalidated, the solver rebuilds the CPU lists or CUDA grids and captures new references. This rebuild does not reset density, refresh viscous prior force, or change the advection phase; physical kernel support remains `2h`.

`SPHStateStatistics` and `SPHKinematicsStatistics` are transient reduction results used for time-step and search-buffer calculations. `SPHInteractionParameters` holds immutable values for one CPU stage. Uniform physical configuration remains in `SPHDEM`, which supplies current values whenever a stage runs. Shared equations live in `execution/sphFunctions.h` and `execution/sphCouplingFunctions.h`.

`virtualParticleCoupling` retains stable LS-owner mappings, sampled boundary kinematics, held forces and torques, and new-minus-held load increments. It averages or reduces these values when requested and supports capturing/restoring its host accumulators for temporary observations. `SPHDEM` controls the sampling times, acoustic updates, impulse-correction stage, and complete snapshot lifecycle. CUDA copies use the stream supplied by the solver.

## History Debugging Checklist

- Verify master/slave indices and the slave history key.
- For LS-LS, verify that the mapping index identifies the geometry container's surface node.
- Ensure the active prefix sum's final value matches device contact size.
- Do not retain a neighbor list across a particle-count decrease.
- If tangential response resets unexpectedly, confirm the current contact was matched to a prior history entry rather than treated as new.
- If duplicate contacts appear, inspect pair ownership between grid levels and self-search ranges.
