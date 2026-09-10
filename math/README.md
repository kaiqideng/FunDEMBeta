# `math/`

> **Purpose:** Provide foundational mathematics shared by CPU and GPU code without depending on containers or solvers.

## Core Contents

| File | Responsibility |
| --- | --- |
| `Constants.h` | Constants, default tolerance, finite checks, and approximate comparison |
| `Vector3.h` | `Vec3`, vector algebra, norms, and safe normalization |
| `Matrix3.h` | Row-major `Mat3`, matrix operations, and conditioned inversion |
| `Quaternion.h` | Unit quaternions, rotation, inverse rotation, and normalization |
| `Indexing.h` | Conversion from three-dimensional grid coordinates to x-major linear indices |

## Key Rules

- `linearIndex(x, y, z, sizeX, sizeY)` evaluates `(z * sizeY + y) * sizeX + x`.
- `Mat3` is row-major, and rigid-body inertia tensors are stored in local coordinates.
- `rotateUnit()` and `inverseRotateUnit()` are fast paths that require an already normalized quaternion.
- Prefer `tryNormalize()`, `normalizedOrZero()`, `normalizedOrIdentity()`, and `tryInverse()` when inputs are uncertain.
- Use `defaultTolerance` for degeneracy checks instead of scattering hard-coded tiny values through the algorithms.

## Extension Guide

- Keep public functions free of dynamic allocation, exceptions, and host-only state.
- Shared CPU/GPU functions must compile with both an ordinary C++ compiler and NVCC.
- New types placed in the device SoA should remain standard-layout and trivially copyable.

## Numerical and Coordinate Conventions

- `Real` is `double` on both host and device.
- Vectors use Cartesian components and matrices are row-major.
- Quaternions represent active rotations from a particle-local vector into the world frame.
- Rigid-body inertia tensors are stored in the local frame; world-frame angular updates rotate through the current unit quaternion.
- Angles passed to axis-angle construction are radians. Use `degreesToRadians()` at input boundaries.

`nearlyEqual()` combines absolute and relative tolerances using a scale of at least one. It is appropriate for validation and tests, but it must not replace a physically meaningful convergence tolerance. `defaultTolerance` is intended for numerical degeneracy such as zero-length directions, singular matrices, and normalization guards.

## Safe Operation Guide

| Operation | Preferred interface | Failure behavior |
| --- | --- | --- |
| Normalize uncertain vector | `tryNormalize()` or `normalizedOrZero()` | Returns failure or zero instead of NaN |
| Normalize uncertain quaternion | `normalizedOrIdentity()` | Returns identity for a degenerate input |
| Invert uncertain matrix | `tryInverse()` | Reports an ill-conditioned matrix |
| Compare floating values | `nearlyEqual()` | Uses absolute and relative tolerance |
| Check user input | `isFinite()` | Rejects NaN and infinity explicitly |

Fast methods whose names contain `Unit` assume the caller already supplied a normalized quantity. They avoid redundant normalization inside hot contact and integration loops; use the guarded alternatives at user-input boundaries.

## Device Compatibility

Functions used by CUDA carry `FUNDEM_MATH_HD`, remain allocation-free, and do not throw. When adding a function, use the wrappers in `Constants.h` for math operations that differ between host `<cmath>` and CUDA device compilation. Validate the same function with the ordinary C++ compiler and NVCC whenever it enters a kernel path.
