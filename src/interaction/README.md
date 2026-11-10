# `src/interaction/`

> **Purpose:** Store contacts, bonds, neighbors, and coupling histories, and provide CPU searches and SPH interaction stages.

## Core Contents

| File | Responsibility |
| --- | --- |
| `contact.*` | Contact geometry, effective mass/radius, particle copies, history, force, torque, and energy; pair coefficients are looked up during evaluation |
| `bond.*` | Bond endpoints, stiffnesses, BK damage state, stored force/bond-point torque, and computed master/slave torque accessors |
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
- Particle-container types select the contact family before type erasure; materials carry no geometry marker.
- All contact stiffness, friction and explicit damping come from one symmetric material-pair record. Its two switches select all four ordinary stiffnesses or all four dashpots. With the linear group off, sphere-sphere and sphere-LS use Hertz-Mindlin and set both length-scaled rotational stiffnesses equal to the default sliding stiffness, as a prescribed FunDEM extension. Only restitution-based default damping combines material restitution values.
- LS-LS requires explicitly configured nonnegative per-area normal/sliding stiffness, independent of the ordinary linear-stiffness switch. When the damping group is on, nodal normal/sliding dashpots receive their share of the pair coefficient, weighted by active nodal area. No separate rolling/torsional spring or dashpot is added to LS-LS nodal contacts.
- `effectiveRadius_` is only a physical radius/lever arm; it no longer selects the contact model.
- `geometrySurfaceNodeIndex_` is a global index into the geometry surface-node container; `ownerParticleIndex_` indexes the LS-particle container.

CPU search may parallelize candidate generation and per-contact calculations with OpenMP, but accumulation into shared particles must remain race-free. CPU and GPU sphere–LS paths use identical particle ordering and torque conventions.

### Bond Setup

Configure a bond in this order:

1. Call `setEquivalentLength()`.
2. Call `setConnection()` with particle container(s) and endpoint indices, or with a contact.
3. Call `setStiffness(normal, shear, bending, torsional)` with the four physical stiffnesses.
4. Insert the bond into the solver container matching its particle types.

Changing equivalent length invalidates dependent connection geometry and stiffnesses, which must then be configured again. Each calculation cycle refreshes the connected particle state and clears the previous stored force and bond-point torque before evaluating the bond response. Master/slave particle torques are computed on demand, not retained as separate Bond members or device fields.

Equivalent length is user supplied, independent of endpoint center distance
and the connection reference point. The core accepts four physical
stiffnesses through `setStiffness()`; applications own any conversion from
Beam material and section inputs.

`setStiffness()` converts the physical inputs to the four stored coefficients:
`B1 = kn`, `B2 = ks * L^2`, `B4 = kt`, and
`B3 = kb - B2/4 - B4/2`, where `L` is the equivalent length. Normal and shear
stiffness have units N/m; bending and torsional stiffness have units N·m.
The derived `B3` can be negative even though all four physical inputs must be
finite and nonnegative. Check the Boolean result of every setup call before
inserting the bond; a failed setter does not make a default bond calculable.

The index-based overloads use `setConnection(particles, masterIndex, slaveIndex)` for sphere–sphere or LS–LS bonds, and `setConnection(spheres, sphereIndex, LSParticles, LSParticleIndex)` for mixed bonds. No normal argument is accepted. They reuse sphere-contact geometry with sphere radii or LS bounding radii: `n = (master.position() - slave.position()) / distance`, `overlap = masterRadius + slaveRadius - distance`, and `point = slave.position() + (slaveRadius - overlap / 2) * n`. The point is the center midpoint only when the two radii are equal. Separated bounding spheres remain valid bond endpoints; finite, distinct centers are required, and coincident centers fail because no normal can be defined. Contact-based overloads use the supplied contact's point and normal instead.

`setCrossSectionArea()` assigns the nominal cross-sectional area in square
metres. The matching VTU field is `bondVTUField::crossSectionArea`. Zero area
disables fracture without disabling the elastic response. Set
`fractureMode::unbreakable`, `brittle` or `damage` explicitly: unbreakable
skips fracture, brittle fails at strength onset regardless of stored
toughness, and damage uses the law below. The core default is `damage`,
with zero area preserving the legacy fracture-disabled default. Positive
area requires positive tensile strength and cohesion only when fracture is
enabled. The current `isValid()` retains a positive normal/shear stiffness
guard for damage mode with either positive critical energy; these stiffnesses
are not denominators in the current damage helper. Bending and torsional
stiffness may be zero. Setting the mode does not reset irreversible history.
Software configuration flags are mapped to these core modes by the application.

### Energy-based EVM Damage

The original B1--B4 force, endpoint-couple, and B3/B4 torsion/bending-moment
expressions are unchanged. The undamaged B3/B4 moment is used directly
for the selected bond-point strength convention; it is not replaced by a
two-end envelope or a transported section moment.

The shared CPU/device evaluator retains normal energy and evaluates shear,
bending and torsion using squared vectors equivalent to the original formulas
for unit endpoint frames before the final torsional clamp. Bending uses the
stored original `kb` directly. The original dot-product formulas remain
commented as a reference. The current stored torsional energy is
`max(Ut_raw,0)`, where `Ut_raw` is the tangent-norm difference, not an
unconditionally nonnegative quadratic form.

For the circular section, the undamaged trial force and bond-point moment
give `sigma_a=-dot(F,n)/A`, `sigma_b=|M_perp|*R/I` and
`tau=|F_perp|/A+|dot(M,n)|*R/J`, with `R=sqrt(A/pi)`, `I=A*R²/4`
and `J=2*I`. Strength initiation uses `sigma_criterion=sigma_a+sigma_b`:
the squared tensile ratio plus the squared shear ratio with capacity
`cohesion-friction*sigma_criterion`. Nonpositive shear capacity with
nonzero shear initiates failure directly.

Damage instead uses the distinct opening stress and trial energies:

```text
U_open = (sigma_a>0 ? max(Un,0) : 0) + max(Ub,0)
U_shear = max(Us,0) + max(Ut,0)
sigma_open = max(sigma_a,0) + sigma_b
dn = sigma_open>0 ? 2*U_open/(A*sigma_open) : 0
ds = tau>0 ? 2*U_shear/(A*tau) : 0
d = hypot(dn,ds)
```

Energies are in J, stresses in Pa, area in m² and equivalent lengths in m.
The current helper takes no normal/shear stiffness arguments. Compression
can reduce the strength-criterion stress but cannot cancel the opening
bending term. These lengths are energy-based approximations, not geometric
gaps or a direct implementation of Suarez et al. (2025).

At the first sampled strength crossing, BK `Gc` is computed from the
`U_open`/`U_shear` ratio and frozen. Maximum-scaled weights avoid an absolute
small-energy cutoff and overflow in their unscaled sum. The onset history
stores `sigma_e0=hypot(sigma_open,tau)` and `d0=d`; the failure length is
`df=2*Gc/sigma_e0`. This is not the old radial traction
`2*(U_open+U_shear)/(A*d0)`. `getOnsetEffectiveStress()` returns `sigma_e0`.
Evolution uses `dmax=max(previousMaximum,d)` and
`D=max(previousD,(1-d0/dmax)/(1-d0/df))`, with failure at `dmax>=df`.
There is no within-step onset interpolation or healing.

Brittle mode and damage mode with both critical energies zero check
strength directly. Otherwise zero onset `Gc`, zero onset `d`, or `df<=d0`
causes immediate failure without extending the configured toughness.
The last condition is not an independent `A*Gc<=U0` test.

For `D<1`, pure axial compression and its stored energy are retained; the
remaining force/couple and energy contributions are multiplied by `(1-D)`.
At full failure the entire Bond load and energy state is cleared, leaving any
active native contacts. An axial-compression energy release at that instant
also releases stored compression energy. Removal is not an energy-conserving
transfer to native contacts. The current stress-based equivalent length and
failure length do not establish general mechanical-work equality with
`A*Gc`; work acceptance requires an independent reference for the loading
path, not the former radial-traction budget argument.

All bending energy is assigned to Mode I, including compressed bending;
this does not resolve tensile/compressive fibres or guarantee crack closure.
Torsion shares Mode-II toughness, and EVM has higher-order modal coupling.
See [the Bond guide](../../docs/UserGuide.md#94-strength-initiation-and-bk-softening)
for full definitions. Host `bond::calculateForce()` and CUDA
`bondForceKernel.cu` both use the shared elastic, damage and response-scaling
helpers before assembling force lever-arm moments. This reuse is not a
claim that the current law has passed a parity or work campaign. The
[test guide](../tests/README.md) distinguishes current implementation from
historical constitutive-test expectations; failed campaigns are not rewritten
as passes.

### Container Lifecycle

`interactionContainer` sizes particle ranges from the particle count and, where
needed, surface-node ranges from the expanded node count. Growing a population
preserves accepted contacts and their spring histories. A count decrease drops
rebuildable search/range allocations, not the accepted contact and Bond objects.
The final prefix-sum value gives the active contact count, and contact storage
grows only when needed.

Runtime removal must use the solver's `removeSpheres()` or
`removeLSParticles()` operation, not direct host-vector erasure. The returned
old-to-new map contains `-1` for removed particles and dense indices for
survivors. `remapParticles()` removes incident contacts, copies each surviving
contact with its spring state, and remaps surviving endpoint indices. A Bond
incident to a removed particle is rejected; surviving Bonds retain their
reference geometry, stiffness, energies and irreversible damage history. This
restriction is why applications must not attach Bonds to temporary particles
that will later be removed.

The subsequent runtime initialization refreshes forces without advancing
contact springs or Bond damage, rebuilds disposable search storage, and
uploads the compacted particle and interaction state when CUDA is active.
This is a same-solver continuation, not a restart from initial conditions.

After a GPU solve, only valid contacts and bonds are copied to the host. Neighbor ranges, temporary lists, and neighbor history are rebuildable data. Contact history survives normal repeated initialization, while neighbor history is reset each time.

## Extension Guide

- Configure contact geometry separately from the symmetric material-pair table; pass that table and the contact family to force evaluation.
- Keep host particle copies in contact and bond objects out of `DeviceLayout`.
- Put new shared CPU/GPU formulas in `execution/`; object methods should only organize inputs and retain results.
- For every new contact type, explicitly define master/slave ordering, history keys, and torque distribution.

## Interaction Data Model

| Type | Persistent role | Device role |
| --- | --- | --- |
| `contact` | Current geometry, particle copies, history state, loads, and energies; no persistent material-pair coefficients | One active detected contact |
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

For CPU LS–LS detection, the parallel threshold is 4096 candidate master-surface-node queries summed over all pairs, not 4096 particles or detected contacts. Ordered pair/node ranges are split into blocks of at most 1024 nodes and scheduled at a single OpenMP level. Each pair's rotations and relative position are calculated once per search; its blocks use local-coordinate SDF queries and reusable contact buffers. The CPU-only transform cache preserves the shared CPU/CUDA world-space query interface.

The contact broad phase uses the shared [`cpu::pointGrid`](../execution/cpu/README.md#fixed-radius-point-queries)
occupied-cell index. Its default cell size, z/y/x query order, per-cell original
indices, bounding-sphere inclusion, and direct-search option are unchanged.
The same component provides independent fixed-radius center-distance pairs for
Bond preparation without narrowing them to touching particles or sampling SDFs.

Blocks are merged in pair and surface-node order regardless of worker scheduling. Effective mass is normalized over the **whole active contact area of the pair**, never separately per block: `nodal effective mass = pair effective mass × nodal area / total active pair area`. The area sum follows surface-node order, and spring restoration retains the same `(master index, slave index, master surface-node index)` key. Changing the work partition therefore does not renumber history entries or treat each block as a separate physical contact pair.

## GPU Storage Lifecycle

1. `initializeDevice()` sizes range arrays from particles and, for LS particles, generates surface-node mappings.
2. Existing bonds are copied to device storage.
3. Existing host contacts are converted into compact contact/range history, and the active contact objects are uploaded as well.
4. Detection kernels write neighbor counts.
5. Prefix sums produce total candidate or contact counts.
6. Contact or neighbor storage grows only when the required count exceeds capacity.
7. Detection writes active contacts and force kernels consume them.
8. `saveCurrentStepToHistory()` stores the current keys and spring states for the next rebuild.

Growing particle counts preserves current contacts and histories while
expanding required ranges. Decreasing particle or surface-node counts resets
search storage; an explicit runtime index map first updates accepted
interactions as described above. Compact GPU history is reconstructed from
those surviving host contacts. Neighbor-range history is rebuildable and
resets during normal initialization; it is not the authoritative spring state.

## Force and Torque Assembly

Every contact produces one force and torque record before loads are accumulated into particles. For sphere-LS and LS-LS contact, torque signs follow the stored master/slave convention. CPU assembly uses race-free reductions or staged accumulation; CUDA uses atomics only where multiple interactions write one particle.

Bond calculations first refresh the connected particle state and clear the previous stored force and bond-point torque. A bond is calculable only after equivalent length, connection geometry, and all four stiffnesses have been accepted. Damage state persists across calls and across normal device reinitialization.

`masterTorque()` and `slaveTorque()` return vectors by value, deriving the
original EVM endpoint/shear couples and force lever-arm moments from the
current configuration. Their connected particle containers must remain alive;
the accessors read the current owning host particles, not the cached particle
copies used by CPU force evaluation. GPU callers must first synchronize both
the host particle and Bond state. CUDA
keeps the two endpoint torques as local kernel values and scatters them
directly into particle torque arrays, without persistent Bond torque arrays.
Optional master/slave torque VTU fields remain available as derived values;
their field names and output schema do not change.

## SPH Neighborhoods and Coupling

`sphInteraction.h` is private implementation and is not installed as a user API. Its CPU engine owns two compact host neighbor lists: fluid candidates per fluid particle and virtual-wall candidates per fluid particle. Both are constructed from sorted occupied-cell uniform-grid hashes at advection preparation and rebuilt earlier when displacement invalidates reuse. CUDA instead queries its fluid and virtual-wall grids directly. Neither representation changes particle-container ordering.

`SPHNeighborhood` stores the fluid and boundary positions at the most recent search build separately from physical particle state. Before a density or pressure stage, `SPHDEM` checks that the maximum actual displacement is no greater than half the search skin. If that bound is exceeded or references are invalidated, the solver rebuilds the CPU lists or CUDA grids and captures new references. This rebuild does not reset density, refresh viscous prior force, or change the advection phase; physical kernel support remains `2h`.

[`SPHInteractionTypes.h`](SPHInteractionTypes.h) defines backend-independent `SPHStateStatistics`, `SPHKinematicsStatistics` and `SPHInteractionParameters` in `fundem`. Existing `cpu::` and `cuda::` statistics names remain aliases, not separate copies. Reductions return fluid or boundary extrema; the solver retains the latest fluid statistics for time-step and search-buffer calculation. `SPHInteractionParameters` holds immutable values for one stage, constructed from current uniform properties, gravity and the optional periodic domain. Physical configuration remains in `SPHDEM` rather than being duplicated in cached stage parameters. Shared equations live in `execution/sphFunctions.h` and `execution/sphCouplingFunctions.h`.

With an enabled SPH periodic domain, CPU and CUDA neighbor stages use
minimum-image separations for both fluid and virtual-wall samples, including
free-surface classification, density reset/rate, pressure and viscosity.
Periodic grid queries visit each cell at most once even when the search skin
covers an entire axis; every physical particle/sample keeps its original index
and mass/volume. Wall reaction is accumulated to the original LS owner with
its physical local torque arm, not a period-length image displacement.
Default-disabled domains preserve ordinary nonperiodic searches. Device grid
and interaction launchers receive the caller's `cudaStream_t` explicitly.

`virtualParticleCoupling` retains stable LS-owner mappings, sampled boundary kinematics, and held world forces and torques. At each common acoustic boundary it refreshes geometry from the current LS pose and prepares wall velocity/acceleration from the completed solid interval's samples, or from instantaneous kinematics at startup. `SPHDEM` advances SPH first, collects the new reaction, and installs it before the first DEM half kick. Both halves of every solid substep then use the same held hydraulic load. Replacing the hydraulic part of the cached total does not reevaluate contact/bond history or call external hooks a second time. There are no new-minus-held increment caches or endpoint impulse-correction stage. CUDA copies and launchers use the stream supplied by the solver.

Current-time observations are owned by `SPHDEM`: they expose dense-output fluid state and current solid-boundary geometry without changing these kinematic accumulators, the held reaction, search references, or device state. The interaction component does not select output times or maintain a second observation-time integration path.

## History Debugging Checklist

- Verify master/slave indices and the slave history key.
- For LS-LS, verify that the mapping index identifies the geometry container's surface node.
- Ensure the active prefix sum's final value matches device contact size.
- Do not retain a neighbor list across a particle-count decrease.
- If tangential response resets unexpectedly, confirm the current contact was matched to a prior history entry rather than treated as new.
- If duplicate contacts appear, inspect pair ownership between grid levels and self-search ranges.
