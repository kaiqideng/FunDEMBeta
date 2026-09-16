# `src/geometry/`

> **Purpose:** Store compact level-set geometry that can be shared by multiple `LSParticle` objects.

## Core Contents

| Type | Responsibility |
| --- | --- |
| `LSGeometryDescriptor` | Store grid origin, inverse spacing, `int3` size, three offsets, bounding radius, volume, and unit-density inertia |
| `LSGridNode` | Store one signed-distance value |
| `LSSurfaceNode` | Store a local position and nodal area |
| `LSSurfaceTriangle` | Store local `int3` connectivity |
| `LSGeometryContainer` | Combine the four Host AoS/Device SoA containers above |

## Key Rules

- `LSInfo::buildLSGrid(..., isFixed=false)` integrates volume, centroid, and unit-density inertia, shifts both the grid origin and surface nodes into the centroid frame, and caches the resulting bounding radius and mass properties.
- `LSInfo::buildLSGrid(..., isFixed=true)` preserves the input local frame and stores zero volume and inertia; an `LSParticle` using this geometry has infinite mass.
- Without an explicit `isFixed` argument, particle shapes default to movable geometry and wall types default to fixed geometry, including through `LSInfo&`. `BoxParticle` defaults to movable despite inheriting from `BoxWall`; passing `false` explicitly allows movable wall geometry.
- `add(const levelset::LSInfo&)` validates and copies the built geometry and its cached metadata, then computes nodal areas. It does not integrate mass properties or shift coordinates.
- With triangles, one third of each triangle area is accumulated into each of its three nodes.
- Without triangles, the bounding-sphere area is distributed uniformly as a fallback nodal area.
- Device geometry is read-only; host geometry changes require another host-to-device copy.
- Geometry indices and offsets depend on insertion order. Do not reorder or erase geometry after particles reference it.

## Extension Guide

`data/myLSObject.*` generates `levelset::LSInfo`, which owns the geometry arrays, local frame, bounding radius, volume, and unit-density inertia. The solver's `addGeometry(const levelset::LSInfo&)` passes that built input to `LSGeometryContainer::add()` for compact storage. Neither insertion API takes a fixed flag; `buildLSGrid()` selects the shape's default mode or an explicit override.

## Four-Array Layout

```text
descriptor[i]
    signedDistanceOffset_  -> gridNodes[]
    surfaceNodeOffset_     -> surfaceNodes[]
    surfaceTriangleOffset_ -> surfaceTriangles[]
```

The next descriptor's offset, or the end of the corresponding array for the final geometry, defines each range end. Counts are deliberately not stored. This keeps descriptors compact and gives CPU and GPU code the same range convention.

`gridNodeSize_` is an `int3`; signed-distance storage is x-major and uses `linearIndex(x, y, z, sizeX, sizeY)`. Grid spacing is stored as its inverse because contact detection performs interpolation repeatedly. Surface triangles store local `int3` node connectivity, while each surface node stores a local position and its integration area.

## Insertion Pipeline

1. Validate grid dimensions, spacing, signed distances, surface nodes, and connectivity.
2. Read the already-built local frame and cached `boundingRadius()`, `volume()`, and `unitDensityInertiaTensor()` from `LSInfo`.
3. Accumulate triangle area to nodes, or use the documented fallback when triangles are absent.
4. Append one descriptor followed by the three variable-length ranges, preserving the input coordinates and mass properties.

The returned geometry index is stable. `LSParticle::setGeometry()` copies the host ranges required by CPU contact detection and stores the central geometry index used by device kernels.

## Resolution Guidance

The signed-distance grid must contain at least two nodes per axis. Finer spacing improves interpolation, volume, inertia, and contact accuracy but increases host copies, device memory, and LS-LS search work. For SPH walls, the boundary grid spacing must not exceed the configured smoothing length. Inspect a generated grid with `LSInfo::outputGridVTI()` before running an expensive case when sign convention or padding is uncertain.

Fixed geometry deliberately keeps zero volume and inertia. Use the `isFixed` flag in `LSInfo::buildLSGrid()` only for prescribed or infinite-mass boundaries; using it for a movable body prevents mass construction. For example, use `wall.buildLSGrid(50, true)` for resolution-based construction or `wall.buildLSGrid(spacing, 2, true)` for explicit spacing and padding, then call `simulation.addGeometry(wall)`.
