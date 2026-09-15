# `src/`

This directory contains FunDEMBeta's eight library modules and the test suite. Each module README documents its types, ownership rules, and extension points. Module-prefixed source paths in those READMEs, such as `execution/sphFunctions.h`, are relative to `src/`; bare filenames name files within the module.

| Directory | Responsibility |
| --- | --- |
| [`data/`](data/README.md) | Host Array-of-Structures and device Structure-of-Arrays containers, level-set input objects, and scientific output |
| [`execution/`](execution/README.md) | Shared numerical formulas, [CPU container loops](execution/cpu/README.md), and [CUDA kernels and launchers](execution/cuda/README.md) |
| [`geometry/`](geometry/README.md) | Reusable level-set descriptors, signed-distance grids, and surface storage |
| [`interaction/`](interaction/README.md) | Contacts, bonds, searches, SPH neighborhoods, and virtual-boundary coupling state |
| [`material/`](material/README.md) | Ordinary and level-set material definitions |
| [`math/`](math/README.md) | CPU/CUDA scalar, vector, matrix, quaternion, tolerance, and indexing utilities |
| [`particle/`](particle/README.md) | Rigid, level-set, SPH, and virtual particles, inlet descriptors, and spatial grids |
| [`solver/`](solver/README.md) | Simulation construction, lifecycle, execution modes, time integration, and synchronization |
| [`tests/`](tests/README.md) | CPU/CUDA regression tests and physical integration checks |

## Build and include paths

The project entry point is the repository-root [`CMakeLists.txt`](../CMakeLists.txt). Configure from `FunDEMBeta/`, for example:

```bash
cmake -S . -B ../build -DFUNDEM_ENABLE_CUDA=OFF -DBUILD_TESTING=ON
cmake --build ../build --parallel
ctest --test-dir ../build --output-on-failure
```

CMake exposes this `src/` directory as the build include root. Public include directives keep their module-relative names:

```cpp
#include "solver/solvers.h"
```

Applications should link `FunDEM::FunDEM` to receive the include paths and transitive dependencies. Installation keeps public headers under `include/FunDEM/<module>/`. A single-configuration build places test executables under `<build>/src/tests/`; CTest is run from the build root.

## Repository-level resources

The [user documentation](../docs/UserGuide.md), [CMake package support](../cmake/README.md), [tutorials](../tutorial/README.md), [examples](../example/README.md), and [validation cases](../validation/README.md) remain outside `src/` at the repository root. Start with the [project README](../README.md) for setup and the first simulation.
