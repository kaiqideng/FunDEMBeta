# FunDEMBeta User Guide

[Back to the project overview](../README.md)

FunDEMBeta is a C++17 discrete-element and SPH–DEM framework with one data model for CPU and CUDA execution. User-facing objects are assembled in host Array-of-Structures (AoS) containers, while the same declarations generate Structure-of-Arrays (SoA) fields for CUDA. Geometry, integration, contact, bond, and SPH formulas can therefore be shared without maintaining separate CPU and GPU physics models.

The current solver hierarchy is:

```text
solver
└── LSDEM
    ├── SphereDEM
    └── SPHDEM
```

- `LSDEM` solves level-set particles and level-set walls.
- `SphereDEM` adds ordinary spheres and sphere–level-set coupling.
- `SPHDEM` adds CPU/GPU SPH and SPH–level-set wall coupling; ordinary spheres are intentionally disabled.

## Contents

**Getting started**

- [1. Requirements](#1-requirements)
- [2. Workspace and build](#2-workspace-and-build)
- [3. Tutorials and tests](#3-tutorials-and-tests)
- [Write your own C++ case](ProgrammingGuide.md)
- [4. Development environment and external projects](#4-development-environment-and-external-projects)

**DEM modeling**

- [5. First simulation: an LSParticle falling onto a plane](#5-first-simulation-an-lsparticle-falling-onto-a-plane)
- [6. Model construction](#6-model-construction)
- [7. Solver and execution modes](#7-solver-and-execution-modes)
- [8. Contact models and search](#8-contact-models-and-search)
- [9. Bonds](#9-bonds)
- [10. User-defined external force and torque](#10-user-defined-external-force-and-torque)

**Execution and data**

- [11. Time stepping and solve lifecycle](#11-time-stepping-and-solve-lifecycle)
- [12. Output, visualization, and energy](#12-output-visualization-and-energy)
- [13. Container and ownership model](#13-container-and-ownership-model)

**SPH, performance, and reference**

- [14. SPH–DEM](#14-sphdem)
- [15. Performance and scaling](#15-performance-and-scaling)
- [16. Source tree and API documentation](#16-source-tree-and-api-documentation)
- [17. Troubleshooting](#17-troubleshooting)

## 1. Requirements

### Required

| Component | Requirement |
| --- | --- |
| CMake | 3.24 or newer |
| C++ compiler | A compiler with complete C++17 support |
| Build system | Make, Ninja, Visual Studio, or another CMake-supported generator |

### Optional

| Component | Purpose |
| --- | --- |
| NVIDIA CUDA toolkit | CUDA DEM and SPH execution |
| OpenMP | Parallel CPU loops |
| ParaView | Viewing `.vtu` and `.vti` files |
| VS Code C/C++ extension | IntelliSense and declaration/definition navigation |
| VS Code CMake Tools extension | Automatic CMake configuration; convenient but not required |

On Ubuntu or WSL, a basic CPU toolchain can usually be installed with:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build
```

Verify the tools that will actually be used:

```bash
cmake --version
g++ --version
```

For CUDA builds, also verify:

```bash
nvcc --version
nvidia-smi
```

The NVIDIA driver and CUDA toolkit must be mutually compatible. FunDEM does not download either dependency during configuration.

## 2. Workspace and build

### 2.1 Recommended workspace layout

Keep generated files outside the source directory. The repository is designed for a sibling build directory:

```text
workplace/             # terminal working directory throughout this guide
├── FunDEMBeta/          # source only
├── build/              # generated CPU or CUDA build
├── build-cpu/          # optional second configuration
├── install/            # optional installation prefix
└── .vscode/            # editor configuration for the whole workspace
```

Open the terminal in the parent `workplace/` directory and keep it there for every command in this guide, including configuration, compilation, execution, testing, installation, and result inspection. `-S FunDEMBeta` explicitly selects the repository root containing `CMakeLists.txt`; `-B build` creates the sibling build directory. The core modules are under `FunDEMBeta/src/`.

### 2.2 CUDA-enabled Release build

CUDA is enabled by default, but an explicit option makes the intended configuration clear:

```bash
cmake -S FunDEMBeta -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=ON \
    -DFUNDEM_ENABLE_OPENMP=ON \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DBUILD_TESTING=ON

cmake --build build --parallel
```

If CMake cannot find a CUDA compiler, it reports that CUDA is unavailable and
still produces the CPU modules; CUDA tests are then skipped, and GPU/Hybrid
execution is unavailable in that build. The native C++ template in `example/`
builds as `fundem_example` when `FUNDEM_BUILD_EXAMPLES=ON`.

If the toolkit is installed but `nvcc` is not on `PATH`, configure a new build directory with the compiler path explicitly. For example, for a CUDA 13.2 installation in `/usr/local`:

```bash
cmake -S FunDEMBeta -B build \
    -DFUNDEM_ENABLE_CUDA=ON \
    -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.2/bin/nvcc
```

To select a GPU architecture explicitly, pass its CMake CUDA architecture number:

```bash
cmake -S FunDEMBeta -B build \
    -DFUNDEM_ENABLE_CUDA=ON \
    -DFUNDEM_CUDA_ARCHITECTURES=86
```

For more than one architecture, quote the semicolon-separated list, for example `-DFUNDEM_CUDA_ARCHITECTURES="80;86"`. An empty value uses the CMake/compiler default.

### 2.3 CPU-only Release build

This configuration never enables the CUDA language. OpenMP is used only when its compiler/runtime passes the required minimum/maximum reduction check; otherwise host loops remain serial.

```bash
cmake -S FunDEMBeta -B build-cpu \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DFUNDEM_ENABLE_OPENMP=ON \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DBUILD_TESTING=ON

cmake --build build-cpu --parallel
```

To use Ninja explicitly, add `-G Ninja` to the first configure command. Do not change the generator inside an existing build directory; create a new build directory instead.

### 2.4 Debug builds

CPU debugging:

```bash
cmake -S FunDEMBeta -B build-debug \
    -DCMAKE_BUILD_TYPE=Debug \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DBUILD_TESTING=ON

cmake --build build-debug --parallel
```

CUDA device debugging:

```bash
cmake -S FunDEMBeta -B build-cuda-debug \
    -DCMAKE_BUILD_TYPE=Debug \
    -DFUNDEM_ENABLE_CUDA=ON \
    -DFUNDEM_ENABLE_CUDA_DEVICE_DEBUG=ON \
    -DFUNDEM_ENABLE_CUDA_LINEINFO=OFF \
    -DBUILD_TESTING=ON

cmake --build build-cuda-debug --parallel
```

`FUNDEM_ENABLE_CUDA_DEVICE_DEBUG=ON` enables expensive device-debug flags and should not be used for performance measurements. Release CUDA builds include line information by default so profilers can map kernels to source lines.

### 2.5 CMake options

| Option | Default | Meaning |
| --- | ---: | --- |
| `FUNDEM_ENABLE_CUDA` | `ON` | Enable CUDA when a CUDA compiler is available |
| `FUNDEM_ENABLE_CUDA_DEVICE_DEBUG` | `OFF` | Compile CUDA tests and kernels with device-debug information |
| `FUNDEM_ENABLE_CUDA_LINEINFO` | `ON` | Add CUDA source line information when device debug is off |
| `FUNDEM_ENABLE_OPENMP` | `ON` | Use OpenMP when the compiler/runtime supports the required reductions; otherwise use serial loops |
| `FUNDEM_INSTALL` | `ON` | Generate installation and exported CMake package rules |
| `FUNDEM_BUILD_TUTORIALS` | `ON` | Build the four short CPU-default tutorials |
| `FUNDEM_BUILD_EXAMPLES` | `ON` (top-level); `OFF` (embedded) | Build the reusable native C++ template as `fundem_example` |
| `FUNDEM_BUILD_VALIDATIONS` | `ON` | Build dedicated reference cases and ONDEM checks; full runs remain explicit commands |
| `FUNDEM_CUDA_ARCHITECTURES` | empty | Override the CUDA architectures selected by CMake/compiler |
| `BUILD_TESTING` | `ON` | Build tests registered with CTest |
| `CMAKE_BUILD_TYPE` | generator-dependent | Common values are `Release`, `RelWithDebInfo`, and `Debug` |

Useful CMake operations after the first configuration:

```bash
# Rebuild everything
cmake --build build --parallel

# Build one target
cmake --build build --target tutorial0 --parallel

# Show cached configuration values
cmake -LA -N build
```

## 3. Tutorials and tests

### 3.1 Run the tutorials

```bash
build/tutorial/tutorial0
build/tutorial/tutorial1
build/tutorial/tutorial2
build/tutorial/tutorial3
```

The four small tutorials demonstrate sphere rebound, asymmetric LS contact, a
five-sphere bonded beam and a 960-particle SPH block. They run 0.1–0.3 seconds of
physical time and write 21 or 31 frames. See the
[tutorial guide](../tutorial/README.md) for parameters and output.

Each tutorial includes `solver/solvers.h` and keeps the model and run settings
in one annotated `main()`; tutorial3 also includes `<cmath>` for `std::sqrt`.
Edit the source, rebuild its target and run without arguments. Change
`solve(...)` for the step count, `setTimeStep()` for the integration step,
`setOutputStepInterval()` for the output cadence, and `setVTUOutputFormat()`
for binary or ASCII output. Choose a fresh directory with
`setOutputDirectory()` before rebuilding if existing results must be retained;
the first `solve()` clears old VTU/DAT files inside that directory.

All bundled tutorials select CPU and CUDA device 0, even in a CUDA-enabled
build. To use a compatible CUDA device, change the constructor's mode/device
arguments; tutorial1 instead sets its mode with `setExecutionMode()`. Pure
LSDEM supports CPU/GPU. The other tutorial solvers also support Hybrid. See
[execution modes](#7-solver-and-execution-modes) for the solver APIs.

Tutorial3 is a small teaching case, not the previous experimental bore-in-a-box
benchmark. It has no reaction-force recorder or post-run filter. The original
large tutorials have been replaced.

### 3.2 Write your own native C++ case

Start with the [hands-on programming guide](ProgrammingGuide.md). It walks
through building the minimal template, changing its parameters, registering a
separate CMake target, adding particles and a wall, selecting a backend, and
reading synchronized state in your own code.

[`example/main.cpp`](../example/main.cpp) contains one complete `main()`;
the two-line CMake file links `FunDEM::FunDEM`. Run from the parent directory
containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DCMAKE_BUILD_TYPE=Release -DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build --target fundem_example --parallel
OMP_NUM_THREADS=2 build/example/fundem_example
```

The default is a two-sphere Hertz collision on CPU. Edit source parameters
and rebuild; the example and tutorials take no command-line options.
For a separate case, copy `example/`, give its executable a unique name, and
register the copied directory inside the core's examples CMake block, as
shown in the programming guide. Keep the fresh-output-folder code so each
run preserves earlier results.

On Windows, set `$env:OMP_NUM_THREADS="2"` rather than using the shell prefix.
With a multi-configuration generator, build with `--config Release` and run
`build/example/Release/fundem_example.exe`. The
[example README](../example/README.md) remains the short quick reference.

Each LAMMPS/PFC comparison has its own C++ source and executable target.
See the [validation guide](../validation/README.md) for target names and commands.

### 3.3 Run the tests

```bash
ctest --test-dir build --output-on-failure
```

List tests without running them:

```bash
ctest --test-dir build -N
```

CPU unit tests are available whenever `BUILD_TESTING=ON`. CUDA builds additionally register host/device container and device state-consistency tests. Their sources are under [`src/tests/`](../src/tests/README.md), and a single-configuration build places their executables under `build/src/tests/`.

The Bond regressions use the current shared-helper signatures and stress-based
equivalent lengths. Enabling or listing tests is not evidence that the complete
suite builds or passes: retain the actual build and CTest logs for the source
being delivered. The [test guide](../src/tests/README.md#current-bond-law-and-regression-scope)
separates current constitutive checks from historical campaign outcomes.
For a first tutorial-only build, use `BUILD_TESTING=OFF` as in the root quick
start; do not change the mechanics merely to satisfy a regression expectation.

## 4. Development environment and external projects

### 4.1 VS Code configuration

FunDEM does not require clangd. The Microsoft C/C++ extension (`cpptools`) can provide completion and declaration/definition navigation as long as it reads the `compile_commands.json` generated by the same CMake configuration that builds the code.

#### Open the parent workspace

For the recommended layout, open `workplace/` in VS Code. Place this in `workplace/.vscode/settings.json`:

```json
{
    "terminal.integrated.cwd": "${workspaceFolder}",
    "cmake.sourceDirectory": "${workspaceFolder}/FunDEMBeta",
    "cmake.buildDirectory": "${workspaceFolder}/build",
    "cmake.configureOnOpen": true,
    "cmake.configureSettings": {
        "FUNDEM_ENABLE_CUDA": true,
        "BUILD_TESTING": false
    },
    "C_Cpp.default.compileCommands": "${workspaceFolder}/build/compile_commands.json",
    "C_Cpp.intelliSenseEngine": "default",
    "C_Cpp.workspaceParsingPriority": "high",
    "files.associations": {
        "*.h": "cpp",
        "*.cu": "cuda-cpp",
        "*.cuh": "cuda-cpp"
    }
}
```

If only cpptools is installed, run the CMake configure command manually before opening source files. CMake Tools is not required to consume `compile_commands.json`.

A minimal `workplace/.vscode/c_cpp_properties.json` is:

```json
{
    "configurations": [
        {
            "name": "FunDEMBeta",
            "compilerPath": "/usr/bin/g++",
            "compileCommands": "${workspaceFolder}/build/compile_commands.json",
            "browse": {
                "path": [
                    "${workspaceFolder}/FunDEMBeta/src",
                    "${workspaceFolder}/FunDEMBeta/tutorial",
                    "${workspaceFolder}/FunDEMBeta/validation"
                ],
                "limitSymbolsToIncludedHeaders": false,
                "databaseFilename": "${workspaceFolder}/build/browse.vc.db"
            },
            "cppStandard": "c++17",
            "intelliSenseMode": "linux-gcc-x64"
        }
    ],
    "version": 4
}
```

#### If the source directory itself is opened

When `FunDEMBeta` is the VS Code workspace root, use these settings so new terminals still start in the parent `workplace/` and the compile database resolves to the sibling build:

```json
{
    "terminal.integrated.cwd": "${workspaceFolder}/..",
    "C_Cpp.default.compileCommands": "${workspaceFolder}/../build/compile_commands.json"
}
```

#### Restore declaration/definition navigation

Check these in order:

1. Confirm that `build/compile_commands.json` exists after CMake configuration.
2. Confirm that the paths inside it belong to the current WSL/Linux filesystem and the current source checkout.
3. Run **C/C++: Log Diagnostics** and verify that cpptools reports the expected compile database.
4. Run **C/C++: Reset IntelliSense Database** after changing the build path or compiler.
5. Reload the VS Code window.
6. Reconfigure CMake if a new `.cpp`, `.cu`, or `.cuh` file was added to a `CMakeLists.txt`.

An `includePath` alone is not an adequate replacement for the compile database because it does not reproduce target-specific definitions, CUDA compilation mode, language standard, or generated include paths.

### 4.2 Install

```bash
cmake -S FunDEMBeta -B build-install \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PWD/install" \
    -DFUNDEM_ENABLE_CUDA=ON \
    -DFUNDEM_INSTALL=ON \
    -DFUNDEM_BUILD_TUTORIALS=OFF \
    -DFUNDEM_BUILD_EXAMPLES=OFF \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF

cmake --build build-install --parallel
cmake --install build-install
```

The installation exports one public target, `FunDEM::FunDEM`, together with the module libraries, headers, and CMake package files. Public headers are installed under `include/FunDEM/<module>/`; source organization under `src/` does not add a `src` component to installed paths or to include directives such as `#include "solver/solvers.h"`.

### 4.3 External application

Place the external application's `main.cpp` and the following `CMakeLists.txt` in `workplace/my_fundem_case/`:

```cmake
cmake_minimum_required(VERSION 3.24)
project(my_fundem_case LANGUAGES CXX)

find_package(FunDEM CONFIG REQUIRED)

add_executable(my_fundem_case main.cpp)
target_link_libraries(my_fundem_case PRIVATE FunDEM::FunDEM)
target_compile_features(my_fundem_case PRIVATE cxx_std_17)
```

From `workplace/`, configure the application against `workplace/install/` and place its generated files in `workplace/build-case/`:

```bash
cmake -S my_fundem_case -B build-case \
    -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build build-case --parallel
```

When an application defines its own CUDA kernels or compiles `.cu` sources, declare `LANGUAGES CXX CUDA` in the application project. A package installed with CUDA also asks CMake to find `CUDAToolkit`; a CPU-only package does not. If the CUDA toolkit is outside `PATH` and CMake's default search locations, pass its root when configuring the consuming project:

```bash
cmake -S my_fundem_case -B build-case \
    -DCMAKE_PREFIX_PATH="$PWD/install" \
    -DCUDAToolkit_ROOT=/usr/local/cuda-13.2
```

## 5. First simulation: an LSParticle falling onto a plane

This complete LSDEM example shows direct use of the public API. It is a
separate model from the compact SphereDEM template in `example/` and runs on
CPU by default:

```cpp
#include "data/myLSObject.h"
#include "material/materialPair.h"
#include "particle/LSParticle.h"
#include "solver/LSDEM.h"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        using namespace fundem;

        constexpr math::Real timeStep = 2.0e-5;
        LSDEM simulation;

        const int particleMaterialIndex =
            simulation.addMaterial(material{1.0e9, 0.25, 0.55, 1200.0});
        const int wallMaterialIndex =
            simulation.addMaterial(material{1.0e9, 0.25, 0.55, 1.0});
        materialPair contactProperties{particleMaterialIndex, wallMaterialIndex};
        contactProperties.setLevelSetStiffness(2.4e9, 7.1e8);
        contactProperties.setFrictionCoefficients(0.5);
        simulation.setMaterialPair(contactProperties);

        levelset::Superellipsoid shape{0.055, 0.040, 0.030, 0.45, 0.65};
        shape.buildSurfaceNode(6);
        shape.buildLSGrid(100);
        const int shapeIndex = simulation.addGeometry(shape);

        LSParticle fallingParticle;
        fallingParticle.setPosition({0.0, 0.0, 0.25});
        fallingParticle.setOrientation(math::Quaternion::fromUnitAxisAngle(
            math::Vec3::unitX(), 29.998657884086096 * math::pi / 180.0));
        fallingParticle.setMaterial(simulation.materials(), particleMaterialIndex);
        fallingParticle.setGeometry(simulation.geometries(), shapeIndex);
        simulation.addLSParticle(fallingParticle);

        levelset::PlaneWall plane{{0.0, 0.0, 1.0}, 1.0};
        plane.buildLSGrid(0.1, 2, true);
        const int planeIndex = simulation.addGeometry(plane);

        LSParticle fixedPlane;
        fixedPlane.setMaterial(simulation.materials(), wallMaterialIndex);
        fixedPlane.setGeometry(simulation.geometries(), planeIndex);
        simulation.addLSParticle(fixedPlane);

        simulation.setBoundary({-0.5, -0.5, 0.0}, {0.5, 0.5, 0.45});
        simulation.setGravity({0.0, 0.0, -9.81});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory("falling_particle_files");
        simulation.setOutputStepInterval(2500);
        simulation.setVTUOutputFormat(vtuFormat::binary);

        simulation.solve(250000);
    }
    catch (const std::exception& error)
    {
        std::cerr << "FunDEM case failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
```

For GPU execution, configure a CUDA-enabled core and call
`simulation.setExecutionMode(executionMode::GPU)` before solving. Pure
`LSDEM` does not support Hybrid.

The important ordering is:

1. Add bulk materials, retain their indices, and configure the required material pairs.
2. Fully build the level-set object before adding it as geometry.
3. Assign material and geometry from this solver's containers to the particle.
4. Add the valid particle to the solver.
5. Configure domain, time step, output, and optional GPU execution.
6. Call `solve()`.

## 6. Model construction

### 6.1 Materials

All solid particles use the same bulk material constructor:

```cpp
material solidMaterial{
    youngModulus,
    poissonRatio,
    restitution,
    density};
```

Contact coefficients belong to an unordered material pair, not either material:

```cpp
materialPair pair{firstMaterialIndex, secondMaterialIndex};
pair.setLevelSetStiffness(normalStiffnessPerUnitArea, slidingStiffnessPerUnitArea);
pair.setLinearStiffness(normalStiffness, slidingStiffness, rollingStiffness, torsionalStiffness);
pair.setFrictionCoefficients(slidingFriction, rollingFriction, torsionalFriction);
pair.setDampingCoefficients(normalDamping, slidingDamping, rollingDamping, torsionalDamping);
simulation.setMaterialPair(pair);
```

There are exactly two Boolean switches. `linearStiffnessEnabled` selects all four
ordinary stiffnesses together, only for sphere-sphere and sphere-LS contacts;
`dampingEnabled` selects all four dashpots together. Both default to false.
Group and individual value setters enable their corresponding whole group.
Disabling a group preserves its stored values. For example:

```cpp
pair.setLinearStiffnessEnabled(false); // Hertz/Mindlin; rolling and torsional stiffness each equal default sliding stiffness.
pair.setDampingCoefficients(0.0, 0.0, 0.0, 0.0); // Enables explicit zero damping in every direction.
pair.setDampingEnabled(false); // Restore restitution-derived damping; retain the zeros.
simulation.setMaterialPair(pair);
```

Rules:

- Young's modulus and density must be finite and strictly positive.
- Poisson ratio must lie strictly between -1 and 0.5.
- Stiffness and friction values must be finite and non-negative.
- Restitution must be within `[0, 1]`.
- `material` has no Sphere/LS marker; the same definition may serve either family.
- The solver automatically creates `m*(m+1)/2` pairs, including same-material contacts.
  `simulation.materialPairs().at(i,j)` and `.at(j,i)` refer to the same entry.
- Pair coefficients are final effective values, with no hidden series-stiffness or friction mixing.
- Sphere-sphere and sphere-LS contacts use Hertz-Mindlin when the linear-stiffness group is disabled. The length-scaled rolling and torsional stiffnesses each equal the default sliding stiffness, as a prescribed FunDEM extension rather than a classical Hertz-Mindlin rotational law.
  Sphere-LS uses the detected sphere radius as a local-plane contact approximation,
  not the level-set body's bounding radius. Both materials' elastic compliance is used;
  infinite mass constrains motion, not Young's modulus.
- LS-LS contacts require explicitly configured normal and sliding stiffness per
  unit area, independently of both switches. Initial `-1` values mean unset;
  nonnegative values, including zero, are configured. Initialization reports
  missing material combinations before force assembly.
- Disabled damping falls back to the harmonic combination of material restitution values.
  Enabled zero damping is an explicit override. Explicit LS dashpots are distributed
  over active nodes by `node area / total active pair area`; LS-LS ignores
  rotational dashpots.
- Change material-pair values only before initialization or between solves.
  See [`src/material/README.md`](../src/material/README.md) for complete activation and storage rules.

### 6.2 Level-set geometry

Geometry is built through `fundem::levelset::LSInfo` implementations:

| Type | Main construction interface |
| --- | --- |
| `Sphere` | `Sphere(radius)` then `buildSurfaceNode()` |
| `Superellipsoid` | Five shape parameters then `buildSurfaceNode()` |
| `HexagonalPrism`, `TriangularPrism` | Polygon circumradius and full Z height; refine with `buildSurfaceNode(surfaceSpacing)` |
| `RegularPolygonPrism` | Circumradius, full Z height and at least three sides |
| `Tetrahedron`, `Dodecahedron`, `Icosahedron` | Positive circumscribed-sphere radius; spacing-refined planar faces |
| `Octahedron` | Positive X/Y/Z center-to-vertex distances |
| `Pyramid` | Full rectangular base size and full height |
| `CutSphere` | Sphere radius and native Z cut height; retained solid lies below the cap |
| `Torus` | Major radius and smaller positive tube radius |
| `Link` | Bend radius, smaller positive tube radius and nonnegative straight length |
| `Capsule` | Radius and straight cylindrical length; total Z length includes both caps |
| `TriangleMesh` | Vertex/triangle arrays or `loadOBJ()` |
| `PlaneWall` | Outward normal and finite surface size |
| `BoxWall` | Box size; use `reverseSDFSign()` for an enclosing tank |
| `BoxParticle` | Box size then `buildSurfaceNode(surfaceSpacing)` |
| `CylinderWall` | Bottom center, top center, and radius |
| `ConeWall` | Bottom/top centers and bottom/top radii |

The faceted and ring/capsule shapes use positive physical `surfaceSpacing`,
not an icosphere subdivision level. Surface refinement and Cartesian SDF grid
spacing are independent. `CutSphere(R, 0)` is a lower hemisphere in its native
construction frame; movable grid construction recenters its surface and grid
at the computed mass centroid. A closed `Torus` or `Link` uses a welded ring
mesh, not a star-shaped projection that would fill its hole. `TriangleMesh`
also supplies irregular shapes generated from a superellipsoid base. See the
[input-geometry catalogue](../src/data/README.md) for exact parameter meanings,
limits, mesh construction and random-shape controls.

Build the surface representation first, then build a signed-distance grid before passing the geometry to the solver. Grid construction accepts a resolution across the largest bounding-box extent (default `50`) or explicit spacing with optional padding (default `2`):

```cpp
shape.buildLSGrid(50);          // cells across the largest bounding-box extent
shape.buildLSGrid(0.002, 3);    // explicit spacing and padding
```

Ordinary particle shapes default to movable geometry. `PlaneWall`, `BoxWall`, `CylinderWall`, and `ConeWall` default to fixed geometry, including when called through an `LSInfo&`. `BoxParticle` defaults to movable geometry even though it inherits from `BoxWall`. An explicit `isFixed` argument overrides the shape's default: it is the second argument in the resolution form and the third in the spacing/padding form. Walls can explicitly request movable geometry with `false`.

For movable geometry, `LSInfo` integrates volume, centroid, and unit-density inertia from the grid. It shifts both the surface nodes and grid origin into the centroid frame, then stores the resulting bounding radius and mass properties. These are available before solver insertion through `boundingRadius()`, `volume()`, and `unitDensityInertiaTensor()`.

Useful geometry operations include:

```cpp
shape.reverseSDFSign();
shape.outputGridVTI("shape.vti");
```

Add the built movable geometry:

```cpp
const int geometryIndex = simulation.addGeometry(shape);
```

For an infinite-mass wall, use its fixed default or pass `true` explicitly when building its grid:

```cpp
wall.buildLSGrid(50, true);
const int wallGeometryIndex = simulation.addGeometry(wall);
```

For fixed geometry, `LSInfo` preserves the input local frame and stores zero volume and unit-density inertia. An `LSParticle` using that descriptor consequently has infinite mass and is written to the `fixedLSParticle` output.

`LSDEM::addGeometry(const LSInfo&)` and `LSGeometryContainer::add(const LSInfo&)` accept only the built geometry. The container validates and copies its arrays and cached metadata, and computes nodal areas; insertion does not integrate mass properties or apply another centroid shift.

For an already prepared geometry, use `LSInfo::fromPreparedData()` with its
accepted surface positions/triangles, SDF origin/spacing/size/samples, bounding
radius, volume and unit-density inertia. This data-only restoration preserves
the supplied local frame verbatim: it does not infer an original shape or
repeat centroid correction. Its `signedDistance(localPoint)` interpolates the
saved grid, including boundary nodes, and rejects queries outside it. Such an
object cannot rebuild an analytic surface/grid; retain the original shape
parameters separately when regeneration is needed.

`surfaceAABBMin()` and `surfaceAABBMax()` compute the axis-aligned bounding box
of the current surface points, without padding. They work in the returned local
frame for both native and prepared geometry. Because surface points are already
centroid-corrected when applicable, these queries never subtract the offset
again. Invalid geometry or an empty surface raises `std::logic_error`.

Grid generation uses separate internal sampling bounds in the original shape
frame. In particular, a `PlaneWall` needs a finite three-dimensional sampling
region even when its surface AABB has zero thickness. Before allocating a native
shape's grid, call `estimateGridNodeSize(spacing, padding = 2)`. It shares
`buildLSGrid`'s dimension calculation without resampling, recentering or changing
mass properties. Prepared data has no native shape to regenerate, so estimation
and rebuilding are unavailable.

For an existing grid, read the **entire stored grid extent** directly from
`gridNodeOrigin()` to `origin + spacing * (size - 1)` componentwise. This includes
any saved padding and requires no additional centroid correction. Surface bounds
and grid extent are deliberately distinct; neither changes the existing SDF or
mass-integration algorithm.

### 6.3 Level-set particles

An `LSParticle` becomes valid only after a material and geometry have both been assigned:

```cpp
LSParticle body;
body.setPosition({0.0, 0.0, 1.0});
body.setOrientation(math::Quaternion::identity());
body.setVelocity({0.0, 0.0, 0.0});
body.setAngularVelocity({0.0, 0.0, 0.0});
body.setMaterial(simulation.materials(), materialIndex);
body.setGeometry(simulation.geometries(), geometryIndex);
const int particleIndex = simulation.addLSParticle(body);
```

The material and geometry must belong to the same solver that receives the
particle. The returned index identifies its host-container slot and remains
unchanged by appending particles or rebuilding neighbor grids. Explicit
particle removal compacts slots; use its old-to-new map as described in §11.4.

### 6.4 Spheres

Use `SphereDEM` and the unified `material`:

```cpp
SphereDEM simulation;

const int materialIndex = simulation.addMaterial(material{
    1.0e7, 0.25, 0.5, 2500.0});
materialPair pair{materialIndex, materialIndex};
pair.setFrictionCoefficients(0.3, 0.02, 0.01);
// Leave the linear-stiffness group disabled to use Hertz-Mindlin.
simulation.setMaterialPair(pair);

particle sphere;
sphere.setPosition({0.0, 0.0, 0.2});
sphere.setRadius(0.025);
sphere.setMaterial(simulation.materials(), materialIndex);
const int sphereIndex = simulation.addSphere(sphere);
```

`particle::isValid()` requires a positive radius, a valid assigned material, and finite position, orientation, velocity, and angular velocity. Calling `setInfiniteMass()` after assigning radius and material creates a fixed sphere.

### 6.5 Domain, gravity, and device

```cpp
SphereDEM simulation(executionMode::CPU);
simulation.setBoundary({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
simulation.setGravity({0.0, 0.0, -9.81});
simulation.setTimeStep(1.0e-6);
```

The minimum boundary must be strictly smaller than the maximum boundary in
all three directions. These are search bounds, not physical walls. For GPU,
use `SphereDEM simulation(executionMode::GPU, 0)`. The device index is
non-negative and follows any `CUDA_VISIBLE_DEVICES` remapping. Passing only
`SphereDEM(0)` still selects the default CPU mode.

## 7. Solver and execution modes

### 7.1 Solver choice

| Required particles | Solver |
| --- | --- |
| Only level-set rigid particles/walls | `LSDEM` |
| Spheres, level-set particles, or both | `SphereDEM` |
| SPH fluid and level-set particles/walls | `SPHDEM` |

All solvers use the common `executionMode` enum. `LSDEM`, `SphereDEM`, and `SPHDEM` default to the pure `CPU` mode.

### 7.2 SphereDEM execution modes

`SphereDEM` accepts exactly three modes, selected before `initialize()`, `step()`, or `solve()`:

```cpp
SphereDEM simulation(0);
simulation.setExecutionMode(executionMode::GPU);
```

| Mode | Sphere backend | LSParticle backend | Sphere-sphere | LS-LS | Sphere-LS |
| --- | --- | --- | --- | --- | --- |
| `CPU` | CPU | CPU | CPU | CPU | CPU |
| `GPU` | GPU | GPU | GPU | GPU | GPU |
| `Hybrid` | GPU | CPU | GPU | CPU | GPU |

The removed Sphere CPU-LSParticle GPU combination is not supported. In the mixed mode, the solver copies current LSParticle state to the device for sphere-LS interaction, then returns the resulting LSParticle force and torque to the host and combines them with CPU LS-LS contributions. User code must not duplicate this exchange.

All model objects are appended on the host. Initialization creates the
required device SoA storage and spatial grids. Use `observeCurrentState()`
for a coherent host snapshot across backends (see §11.3). `synchronize()`
only waits for queued device work; it does not copy state to host.
`copyDeviceToHost()` copies during an initialized GPU/Hybrid phase, whereas
the observer also handles resumable post-solve state. Backend changes are
allowed before initialization or after a completed `solve()`, not during an
initialized stepping phase.

### 7.3 SPHDEM combinations

`SPHDEM` accepts the same three execution modes:

```cpp
SPHDEM simulation(0); // Pure CPU by default.
simulation.setExecutionMode(executionMode::GPU);  // CPU, GPU, or Hybrid.
```

| Mode | SPH backend | LSParticle backend | SPH–LS coupling |
| --- | --- | --- | --- |
| `CPU` | CPU | CPU | CPU |
| `GPU` | GPU | GPU | GPU |
| `Hybrid` | GPU | CPU | GPU with host-side resultant transfer |

The CPU path uses a compact neighbor list built from a uniform background grid, parallelized when supported OpenMP is enabled. The GPU path queries the device grid directly. In `Hybrid`, virtual-particle forces are transferred only when SPH advances, reduced to one resultant force and torque per LSParticle, and applied on the CPU at each DEM step; the full LSParticle container is not copied to the device every DEM step.

## 8. Contact models and search

### 8.1 Spatial grid

The CUDA rigid-particle search supports at most four grid levels:

- Finite-mass diameter ratio `<= 4`: one finite-mass level.
- Ratio `> 4` and `<= 16`: two finite-mass levels.
- Ratio `> 16`: three finite-mass levels.
- If infinite-mass particles exist, one level is reserved as the final level for them.

Particles are not reordered in the user container. The grid maintains a separate sorted particle-index array, so existing particle indices, contacts, and bonds remain stable. Infinite-mass particles are placed in the final grid level and are selected as the slave side of generated pairs. CPU contact search applies the same infinite-mass slave convention.

CPU contact search and independent center-distance queries reuse the occupied
cell index in `execution/cpu/pointGrid.h`. `cpu::findPointPairs(points, radius)`
returns exact self pairs with `i < j`; its two-array overload returns local
indices in the respective input arrays. Results are deterministically sorted
and do not reorder either input. Radius zero selects distinct coincident
points. These CPU utilities know nothing about particle radii, masses,
materials, contacts, bonds or periodic boundaries. Use `pointGrid` directly
only when the caller supplies its own exact candidate test; see the
[CPU execution guide](../src/execution/cpu/README.md).

### 8.2 Contact types

- Sphere–sphere uses the complete spherical contact model.
- Sphere–LS uses the complete spherical model, including force, rolling torque, and torsional torque; the level-set particle is the slave.
- LS–LS uses surface-node area weighting.

Contact type is selected from the particle-container types, not a material marker or an inference from `effectiveRadius`.

For sphere–LS contact:

- Both particle families use the same bulk `material`; Young's modulus and Poisson ratio supply both elastic compliances.
- The disabled linear-stiffness group selects Hertz-Mindlin. Sphere-LS uses the sphere radius against a local plane, not an estimated LS curvature.
- The enabled linear group supplies all four stored stiffnesses; the enabled damping group supplies all four stored dashpots. Friction always comes directly from the symmetric `materialPair`, without further mixing.
- Disabled damping uses the harmonic mean of the two restitution coefficients. Explicit zero damping remains zero.
- With the linear group off, rolling and torsional spring stiffness each equal the default sliding stiffness; they cannot be overridden independently. With it on, all four stored stiffnesses apply. The separate damping group supplies explicit rotational dashpots when enabled, subject to the friction limits.

This keeps sphere–sphere and sphere–LS force evaluation on one complete contact path while preserving the area-based LS–LS model.

For positive overlap, the shared force function retains the signed normal
spring-dashpot response, including a possible negative unloading tail; the
former zero-force clamp is disabled. Nonpositive overlap still ends the
contact and clears its sliding history. `normalForceMagnitude()` reports the
absolute normal response for friction/resistance limits and magnitude output,
not its signed projection. A clipped/no-tension model is therefore not
equivalent to the current implementation.

## 9. Bonds

`bond` supports `sphereSphere`, `sphereLSParticle`, and `LSParticleLSParticle` pairs. Create particles first, construct the bond from solver-owned particle containers, then add it to the solver.

Index-based `setConnection()` takes the particle container or containers and the two endpoint indices; callers do not supply a normal. It uses sphere-contact geometry, taking each sphere's physical radius and each LS particle's bounding radius. For center distance `d`, the slave-to-master normal is `n = (master.position() - slave.position()) / d`, overlap is `masterRadius + slaveRadius - d`, and the reference point is `slave.position() + (slaveRadius - overlap / 2) * n`. Equivalently, the point is the center midpoint plus `(slaveRadius - masterRadius) / 2 * n`; it equals the midpoint only for equal radii.

Separated particles or bounding spheres may still be bonded: negative overlap does not reject the connection. The centers must be finite and distinct; coincident centers fail because they do not define a normal. The overloads that take an existing `contact` instead retain that contact's point and normal.

### 9.1 Sphere–sphere bond

```cpp
bond connection{0.01};

if (!connection.setConnection(
        simulation.spheres(),
        masterSphereIndex,
        slaveSphereIndex) ||
    !connection.setStiffness(
        normalStiffness,
        shearStiffness,
        bendingStiffness,
        torsionalStiffness))
{
    throw std::runtime_error("Invalid sphere bond");
}

simulation.addBond(connection);
```

### 9.2 LSParticle–LSParticle bond

```cpp
bond connection{equivalentLength};
connection.setConnection(
    simulation.LSParticles(),
    masterIndex,
    slaveIndex);
connection.setStiffness(kn, ks, kb, kt);
simulation.addBond(connection);
```

### 9.3 Sphere–LSParticle bond

The sphere is the master and the level-set particle is the slave:

```cpp
bond connection{equivalentLength};
connection.setConnection(
    simulation.spheres(),
    sphereIndex,
    simulation.LSParticles(),
    LSParticleIndex);
connection.setStiffness(kn, ks, kb, kt);
simulation.addBond(connection);
```

The equivalent length is a user-supplied mechanical length, independent of
the endpoint center distance and connection reference point. It must be
finite and positive. All four stiffness inputs must be finite and
nonnegative. `setStiffness()` requires that reference length and accepts the
four physical stiffnesses for every particle family. The core has no Beam
input source; applications convert material and section inputs into those
stiffnesses before calling `setStiffness()`.

`setFractureMode()` accepts `fractureMode::unbreakable`, `brittle` or `damage`;
`fractureMode()` reports it. Unbreakable mode retains the elastic response
and positive area for section stress without enabling fracture. Brittle mode
uses strength onset without softening, regardless of retained GIc/GIIc.
Damage mode uses the separation law below. The core default is `damage`
for existing callers; applications should set the intended mode explicitly.

### 9.4 Strength initiation and BK softening

Bonds retain one elastic evaluator and the strength criterion below. The
explicit fracture mode selects unbreakable, strength-only brittle or progressive
behavior; it does not replace the original force/energy formulas. For
progressive damage, configure the mode, section, energies and strengths before
inserting the bond:

```cpp
if (!connection.setFractureMode(fractureMode::damage) ||
    !connection.setCrossSectionArea(area) ||                 // m2
    !connection.setTensileStrength(tensileStrength) ||         // Pa
    !connection.setCohesion(cohesion) ||                      // Pa
    !connection.setInternalFrictionCoefficient(friction) ||   // dimensionless
    !connection.setModeICriticalEnergy(modeI) ||              // J/m2
    !connection.setModeIICriticalEnergy(modeII) ||             // J/m2
    !connection.setModeMixityExponent(exponent))              // positive BK exponent
{
    throw std::runtime_error("Invalid bond fracture parameters");
}
```

All setters reject nonfinite inputs. Strength and friction setters accept
nonnegative values; strengths default to zero (unconfigured). A connected bond
with enabled fracture and positive area requires positive tensile strength
and cohesion; unbreakable mode does not require strengths. The current
`bond::isValid()` also requires positive normal and shear stiffness when
damage mode has either positive critical energy. This is a retained validity
guard, not a stiffness denominator in the current displacement formulas below.
Bending and torsional stiffnesses may be zero. Zero area retains the legacy
fracture-disabled behavior; select `unbreakable` explicitly when a positive
area must remain available without fracture. The software's configuration
flags are translated to this core enum and are not additional core modes.
Changing the enum does not clear irreversible history. Damage is bounded to
[0,1] and never heals during evolution.

This is an **energy-based EVM extension** of the strength-initiation and
irreversible-softening ideas of Suarez et al. (2025), DOI
`10.1016/j.mineng.2025.109488`. It retains FunDEM's original endpoint frames and
B1--B4 elastic response. It does not reproduce the paper's geometric interface
displacement, fixed mixed-mode expression, or multipoint Voronoi discretization.
Native contact and damping laws are unchanged.

`calculateBondForceAndTorque()` records its existing bending/torsion torque
term through its `torsionBendingTorque` output argument, stored as `torque`
in the bond alongside the resultant force:

```text
torque = B3*cross(slaveNormal,masterNormal)
- 0.5*B4*(cross(slaveTangent1,masterTangent1)
+ cross(slaveTangent2,masterTangent2))
```

The endpoint couples and particle-centre torques retain their existing assembly
formulas. `masterTorque()` and `slaveTorque()` return computed vectors by value;
they are not stored members or device arrays. They derive the original EVM
endpoint/shear couples and force lever-arm moments from the current particle
configuration. The connected particle containers must remain alive, and GPU
callers must synchronize the host particle and Bond state before using these
host accessors. The CUDA kernel computes both endpoint torques as local values
and assembles them directly into particle torque arrays. `torque()` exposes the
recorded bond-point term after the current `(1-D)`
reduction; optional `bondVTUField::torque` exports it. A failed bond records
zero. This is an automatic vector state, not a new material parameter.

`calculateBondStrengthCriterion()` uses the UNDAMAGED trial force and this one
`torque`, evaluated once. Its circular equivalent section has
`R=sqrt(A/pi)`, `I=A*R*R/4`, and `J=2*I`. With `M=torque`, it returns
the signed axial stress, nonnegative bending contribution and shear stress:

```text
sigma_a = -dot(F,n)/A
sigma_b = R*norm(M-dot(M,n)*n)/I
tau = norm(F-dot(F,n)*n)/A + R*abs(dot(M,n))/J
sigma_criterion = sigma_a + sigma_b
Q = (max(sigma_criterion,0)/tensileStrength)^2
    + (tau/(cohesion-friction*sigma_criterion))^2
```

The shear term is evaluated only when `tau>0`. A nonpositive shear capacity
with nonzero shear triggers initiation explicitly. The signed
`sigma_criterion` controls both tensile onset and the friction-dependent
shear capacity; it must not be confused with the opening stress used below.
By the selected bond-point convention, the recorded B3/B4 term itself supplies
the moment in this strength calculation. It is not replaced by a transported
endpoint moment or a maximum over two endpoint sections. The combined peak
normal/shear stress approximation is not an integration over section fibres.
For example, particular stiffness combinations can cancel the B3/B4 moment
while endpoint couples remain nonzero; that chosen strength convention remains
in effect.

`setStiffness(kn,ks,kb,kt)` preserves the four original, undamaged inputs as
`normalStiffness()`, `shearStiffness()`, `bendingStiffness()` and
`torsionalStiffness()`, alongside the force/torque coefficients B1--B4.
They are not recovered by coefficient inversion. All inputs and coefficients
are committed together after validation; changing reference length invalidates
both. CPU/device layouts and optional VTU fields retain these original values.
`calculateBondForceAndTorque()` retains the original B1--B4 force and
endpoint-couple expressions. Its recorded torque output uses the original
internal name `torsionBendingTorque`. The four modal energies are:

```text
Un = 0.5*kn*delta^2
Us = 0.5*B2*norm(0.5*(nm-ns)-e)^2
Ub = 0.5*kb*norm(nm+ns)^2
Ut_raw = 0.25*B4*(norm(tm1-ts1)^2+norm(tm2-ts2)^2-norm(nm+ns)^2)
Ut = max(Ut_raw,0)
```

Here `nm/ns` and `tm1/tm2/ts1/ts2` are the endpoint normals and tangents,
`e` is the current master-to-slave unit direction, and `delta` is axial
displacement. Bending uses the stored original `kb` directly. The original
dot-product formulas remain commented in the evaluator as a reference.
The squared-vector forms are equivalent to those reference formulas for unit
endpoint frames before the torsional clamp. The clamp is an explicit part of
the current evaluator: a negative torsional raw difference is stored as zero.
The other three forms are nonnegative for valid nonnegative stiffnesses.

`updateBondDamage()` groups the UNDAMAGED trial energies, in joules, into
opening and shearing contributions. It uses stresses in pascals and area in
square metres to obtain equivalent lengths in metres:

```text
U_open = (sigma_a > 0 ? max(Un,0) : 0) + max(Ub,0)
U_shear = max(Us,0) + max(Ut,0)
sigma_open = max(sigma_a,0) + sigma_b
dn = sigma_open > 0 ? 2*U_open/(A*sigma_open) : 0
ds = tau > 0 ? 2*U_shear/(A*tau) : 0
d = hypot(dn,ds)
```

These are energy-equivalent lengths, not actual bond-point opening/sliding or
the paper's Eq. (10). Normal and shear stiffness are not arguments of this
damage helper and are not divided into the energies. There is no
radius-times-angle term. Axial compression can reduce `sigma_criterion`, but
it does not cancel the bending contribution in `sigma_open` or `U_open`.

At the first sampled `Q>=1`, freeze the BK critical energy release rate:

```text
Gc = GIc + (GIIc-GIc)*(U_shear/(U_open+U_shear))^eta
sigma_e0 = sqrt(sigma_open^2+tau^2)
d0 = d
df = 2*Gc/sigma_e0
dmax = max(dmax,d)
D = max(previousD, (1-d0/dmax)/(1-d0/df))
```

The BK implementation first divides both modal weights by their maximum
before forming the ratio. Their common area scale cancels; exactly zero
weights return zero, and no absolute small-energy cutoff is applied.
`sigma_e0` is the norm of the opening and shear stresses at onset, not
`sqrt(max(sigma_criterion,0)^2+tau^2)` and not the former radial traction
`2*(U_open+U_shear)/(A*d0)`. The stored toughness and onset quantities remain
frozen after initiation. The maximum length and damage cannot decrease;
`dmax>=df` gives complete failure. There is no within-step onset interpolation.

Unbreakable mode skips this update. Brittle mode checks the strength criterion
directly, irrespective of critical energies. Damage mode also uses immediate
strength-only failure when both critical energies are zero. Otherwise, zero
onset `Gc` or zero onset `d` causes immediate failure; a computed `df<=d0`
also leaves no positive softening interval and fails immediately. These
branches set the onset flag and `D=1`; they do not increase the supplied
toughness. The current condition is `df<=d0`, not an independent `A*Gc<=U0`
test.

Before complete failure, the pure axial-compression force and its energy `Uc`
are not reduced. All other force/couple and energy contributions use `(1-D)`.
Force lever-arm moments follow the actual resultant force. With
`Ua=U_open+U_shear`, this corresponds to the recoverable-energy split
`Uc+(1-D)*Ua`; it is not uniform
degradation of the compression branch. At `D=1` the entire Bond is removed:
force, torque, normal and stored energies are cleared, and ordinary native
contacts provide any remaining geometrically active contact response. If
`Uc>0` at removal, this also releases retained compression energy. Native
contact is not an automatic energy-conserving transfer of the removed Bond
spring. The current stress-based equivalent-length and `df` definitions do
not establish a general equality between mechanical fracture work and `A*Gc`;
such an equality requires independent analysis and verification for the
particular loading path. A restarted or initially loaded state also requires
its initial stored energy to be accounted for separately.

The read-only history getters are `getOnsetEffectiveStress()`,
`getOnsetEffectiveDisplacement()`, `getMaximumEffectiveDisplacement()`,
`getFailureEffectiveDisplacement()`, `getEffectiveFractureEnergy()` and
`getDamageStarted()`. The displacement getters describe the ordinary
Euclidean energy-equivalent length `d`; `getOnsetEffectiveStress()` returns
the frozen `sigma_e0` defined above. Neither its product with `d` nor the
former radial-traction work relation should be assumed to be mechanical work.
Optional Bond VTU fields expose these states.

The modal grouping is an approximation: all bending energy is assigned to
Mode I, including bending under axial compression; this is not a fibre-level
crack-closure model. Torsion is grouped with Mode II and shares `GIIc`, rather
than introducing an independent Mode-III toughness. EVM also contains
higher-order coupling between nominal loading modes. Freezing onset mixity is
not general nonproportional mixed-mode fracture validation.

CPU/GPU consistency, loading-increment sensitivity, and independent native
force/torque work checks belong to `src/tests/BondFracture`. Both
`bond::calculateForce()` and the CUDA Bond force kernel call these same shared
host/device helpers and apply damage before assembling endpoint lever-arm
moments. Shared code is not proof of numerical parity. Historical campaigns,
their FAIL records and earlier radial-traction expectations refer to their
own source snapshots, not a pass for this revised law. Current regression
references must use the section-stress equivalent lengths above rather than
the retired radial-traction budget; see the [test guide](../src/tests/README.md).
Prescribed constitutive
loading does not verify a free dynamic integrator or calibrated specimen
behaviour.

## 10. User-defined external force and torque

Contact, bond and coupling forces are assembled before the solver extension
hooks. Gravity is applied separately during velocity integration, so a
particle's stored force does not include its weight. A derived solver can
add fields without replacing the standard force assembly.

CPU sphere example:

```cpp
class DragSphereDEM final : public fundem::SphereDEM
{
public:
    using SphereDEM::SphereDEM;

protected:
    void addSphereExternalForceAndTorque(
        fundem::particleContainer::host_container_type& spheres) override
    {
        for (fundem::particle& sphere : spheres)
        {
            sphere.addForce((-dragCoefficient_) * sphere.velocity());
        }
    }

private:
    fundem::math::Real dragCoefficient_{0.1};
};
```

This override affects host sphere storage only. Use this class on CPU until
you have implemented its device hook; selecting GPU/Hybrid does not translate
the CPU loop into a CUDA kernel.

The level-set CPU hook is:

```cpp
void addLSParticleExternalForceAndTorque(
    fundem::LSParticleContainer::host_container_type&) override;
```

GPU hooks receive device-field views:

```cpp
void addSphereExternalForceAndTorque(
    fundem::particle::device_type, cudaStream_t stream) override;
void addLSParticleExternalForceAndTorque(
    fundem::LSParticle::device_type, cudaStream_t stream) override;
```

Implement device overrides in a CUDA translation unit, launch the user kernel on the supplied stream, and modify the supplied force/torque arrays in place. Override only the overload matching the storage selected for that particle type.

## 11. Time stepping and solve lifecycle

### 11.1 DEM integration

The solver uses shared CPU/GPU motion-integration formulas for velocity, angular velocity, position, and unit-quaternion orientation. Infinite-mass bodies have zero inverse mass and zero inertia/inverse-inertia tensors, so translational and rotational acceleration are skipped.

The user supplies the base DEM time step:

```cpp
simulation.setTimeStep(1.0e-5);
```

Time-step stability remains the case author's responsibility; use the active material-pair law (or elastic moduli for Hertz-Mindlin), smallest mass, contact scale, and intended accuracy when choosing it.

### 11.2 `solve()`, `step()`, and reconfiguration

```cpp
simulation.solve(1000);
```

`solve(numberOfSteps)`:

1. Validates configuration.
2. Initializes CPU and/or GPU storage as required.
3. Writes a due initial frame when periodic output is configured.
4. Advances exactly `numberOfSteps` DEM steps.
5. Synchronizes current device state back to host.
6. Ends the initialized phase so valid new materials, geometry, particles, and bonds may be added.

This allows staged simulations:

```cpp
simulation.solve(firstStageSteps);

// Host state is current here. Existing LSParticle motion may be prescribed by its stable index.
simulation.setLSParticlePosition(particleIndex, newPosition);
simulation.setLSParticleVelocity(particleIndex, newVelocity);
simulation.setLSParticleAngularVelocity(particleIndex, newAngularVelocity);

// New valid objects may also be added using the same solver containers.
simulation.addLSParticle(newParticle);

simulation.solve(secondStageSteps);
```

The three motion setters accept the stable index returned by `addLSParticle()`, validate finite values, and are available before the first `solve()` or after a completed `solve()`. They invalidate the previous integration continuation so CPU and device storage are rebuilt from the modified host state on the next solve. Existing contact history is retained, which supports continuously prescribed moving boundaries.

Time, total step count, output frame count, contact history and existing
particles continue across calls. `solve(1000); solve(500);` advances 1500
additional steps in total. An unchanged next solve resumes cached integration
state; supported model edits invalidate that continuation and rebuild the
affected working storage. Do not edit the model during an initialized
stepping phase.

`step()` initializes on demand and leaves the solver initialized. It advances
one step but does not trigger scheduled output or finish a host snapshot.
A manual stepping loop must call `writeOutput()` at its chosen cadence.
Call `solve(0)` once before such a loop to perform first-solve cleanup and
write any due initial frame. After the loop, `solve(0)` can complete the phase
and synchronize it without advancing time. For ordinary cases, prefer
`solve(numberOfSteps)`.

### 11.3 Inspect state in user code

After adding at least one sphere, use an observer to read synchronized state:

```cpp
math::Vec3 position, velocity;
simulation.observeCurrentState([&](const solver& base)
{
    const auto& state = static_cast<const SphereDEM&>(base);
    const auto& sphere = state.spheres().host().at(0);
    position = sphere.position();
    velocity = sphere.velocity();
});
```

The callback receives a coherent host view at `time()` on CPU, GPU and Hybrid.
Copy the values needed for your output or control logic. Do not retain
references or mutate/advance the solver inside the callback. `synchronize()`
alone does not copy device state. See the
[programming guide](ProgrammingGuide.md#7-read-results-in-your-own-code) for
where to insert this code in the minimal template.

### 11.4 Completed-step resets and particle removal

The following APIs operate between completed advances, never concurrently with
`step()`/`solve()` or from an observation callback:

| API | Effect |
| --- | --- |
| `LSDEM::resetLSParticleVelocities(indices, resetLinear, resetAngular)` | Clear selected LS particles' linear velocity, angular velocity, or both; return only the removed kinetic energy in joules |
| `SphereDEM::resetSphereVelocities(indices, resetLinear, resetAngular)` | The same component-selective reset for spheres |
| `LSDEM::removeLSParticles(indices)` | Compact LS storage and return an old-to-new index map |
| `SphereDEM::removeSpheres(indices)` | Compact sphere storage and remap sphere-sphere and sphere-LS interaction endpoints |
| `SphereDEM::removeLSParticles(indices)` | Remap sphere-LS endpoints as well as the inherited LS-LS storage |

Both reset switches default to `true`; pass `(indices, true, false)` for linear
velocity only or `(indices, false, true)` for angular velocity only. Unselected
components remain unchanged. Indices must be in range and unique. An empty reset
does nothing and returns zero; an empty removal returns the identity map.
Resetting particles whose selected components are already zero also avoids
invalidating continuation. A reset
changes neither pose nor contact springs, Bond damage, forces or torques. It is
a deliberate removal of kinetic energy, not a new damping law or a convergence
test.

Removal preserves survivor order. Each map entry corresponds to an old slot;
`-1` denotes removal, and every other value is its new slot. Drop removed
application references and remap all surviving references. Incident Contacts
are discarded, while surviving Contact states and Bond endpoint indices are
remapped without reconnecting Bonds or resetting their history. Removal of a
particle referenced by a Bond is rejected; it does not silently break that Bond.
The next initialization rebuilds spatial grids, expanded surface mappings and
device search history from the accepted retained state.

These operations use `prepareRuntimeMutation()` to synchronize the represented
completed-step state before a real edit. In SPHDEM, a real edit discards the
private leading acoustic endpoint and restarts its working interval from that
visible state; emitted source ownership and source clocks remain. An unchanged
handoff must not call this function: continuing the same solver preserves its
existing multirate phase. The core does not define desktop Events, temporary
instance ownership, stopping conditions or a separate presentation clock;
those policies belong to the calling application.

## 12. Output, visualization, and energy

### 12.1 Configure output

```cpp
simulation.setOutputDirectory("case_files");
simulation.setOutputStepInterval(1000);
simulation.setVTUOutputFormat(vtuFormat::binary);
```

Binary VTU is the default and is normally preferred for speed and file size.
Use `vtuFormat::ascii` only for readable diagnostics. Set the output directory
before setting its interval. Changing either resets the next due step to the
current step. Scheduled output includes the initial frame but does not force
an off-cadence final frame: at interval 100, `solve(6000)` gives 61 frames,
while `solve(6050)` also gives 61 unless output is explicitly requested.
`writeOutput()` writes a frame immediately and reschedules the next due step
relative to the current one.

Relative output paths are resolved from the executable directory, not from the shell's current working directory. For example, the tutorial's default `tutorial0_files` directory is created beside the tutorial executable.

On the first `solve()` only, FunDEM recursively removes existing `.vtu` and `.dat` files from the configured output directory. It preserves files with other extensions and does not clean again on later `solve()` calls.

Output layout:

```text
case_files/
├── particle/
│   ├── sphere_000000.vtu
│   ├── LSParticle_000000.vtu
│   ├── fixedLSParticle_000000.vtu
│   └── SPHParticle_000000.vtu
├── interaction/
│   ├── sphereContact_000000.vtu
│   ├── sphereLSParticleContact_000000.vtu
│   ├── LSParticleContact_000000.vtu
│   └── ... bond files ...
└── energy.dat
```

Files are written only for particle/interaction groups that exist in the selected solver and current frame.

At every output frame, the console reports calculated step count, current time, frame index, and allocated FunDEM device memory in GB.

### 12.2 Default VTU fields

| Entity | Fields written without user selection |
| --- | --- |
| Sphere | Position, radius, velocity |
| Movable LSParticle | World surface-node positions, triangle topology, surface-node velocity |
| Fixed LSParticle | Same defaults, in a separate `fixedLSParticle` file |
| SPH particle | Position, velocity, mass, density |
| Contact | Point and normal |
| Bond | Point and normal |

`inverseMass` is not a default LSParticle or SPH output field. Movable and infinite-mass level-set particles are separated by file instead.

### 12.3 Optional VTU fields

Select additional fields before solving:

```cpp
simulation.setLSParticleVTUFields({
    LSParticleVTUField::particleIndex,
    LSParticleVTUField::distanceToCenter,
    LSParticleVTUField::force,
    LSParticleVTUField::kineticEnergy});

simulation.setContactVTUFields({
    contactVTUField::force,
    contactVTUField::normalElasticEnergy,
    contactVTUField::slidingElasticEnergy,
    contactVTUField::rollingElasticEnergy,
    contactVTUField::torsionalElasticEnergy});
```

Available fields are:

- Sphere: `orientation`, `angularVelocity`, `force`, `torque`, `inverseMass`, `inertiaTensor`, `materialIndex`, `kineticEnergy`, `gravitationalPotentialEnergy`.
- LSParticle: `particleIndex`, `distanceToCenter`, `orientation`, `angularVelocity`, `force`, `torque`, `inertiaTensor`, `boundingRadius`, `materialIndex`, `volume`, `geometryIndex`, `surfaceNodeArea`, `kineticEnergy`, `gravitationalPotentialEnergy`.
- SPH: `force`, `density`, `densityRate`, `freeSurface`, `kineticEnergy`, `gravitationalPotentialEnergy`.
- Contact: `force`, `torque`, master/slave indices, `overlap`, `area`, `effectiveMass`, `effectiveRadius`, `normalForceMagnitude`, four elastic energies, and three spring deformations.
- Bond: `force`, `torque`, master/slave torque and indices, `equivalentLength`, four original stiffnesses, four elastic energies, global master/slave endpoint normal/tangent bases, `crossSectionArea`, strengths/friction, damage factor and onset flag, mode-I/mode-II critical energy, mode-mixity exponent, onset effective stress/separation, maximum/failure effective separation, and frozen effective fracture energy.

The LSParticle `distanceToCenter` value is the distance from each output surface node to its owning particle center. Bond endpoint orientation fields are transformed to the current global frame before output.
Bond master/slave torque fields are derived from the synchronized host state;
their VTU names and output schema are unchanged by removing the corresponding
persistent Bond members and device arrays.

### 12.4 `energy.dat`

Every output frame appends:

```text
time
kineticEnergy
gravitationalPotentialEnergy
contactNormalElasticEnergy
contactSlidingElasticEnergy
contactRollingElasticEnergy
contactTorsionalElasticEnergy
bondNormalElasticEnergy
bondShearElasticEnergy
bondBendingElasticEnergy
bondTorsionalElasticEnergy
totalEnergy
```

The gravitational value is potential energy, not a per-frame difference. `totalEnergy` is the sum of all columns except time. In `SPHDEM`, these columns intentionally contain only solid-particle, contact, and bond energy; SPH-fluid kinetic, gravitational, and compressibility energy are excluded.

## 13. Container and ownership model

The central storage template is `hostAoSDeviceSoA<T, T::DeviceLayout>`.

```text
Host:   vector<T>                         convenient objects and output
Device: one contiguous array per field   coalesced CUDA access
```

Key semantics:

- `hostSize()` is the number of host objects.
- `deviceSize()` is the number of active device entries.
- Host-to-device copy resizes device capacity when needed and sets `deviceSize()` to `hostSize()`.
- Device-to-host copy resizes the host vector to `deviceSize()` before copying fields.
- Only members declared in `DeviceLayout` are copied to device storage.
- Host-only pointers, validation flags, cached objects, and helper vectors are not device pointers unless explicitly represented in `DeviceLayout`.
- `copyHostToDeviceAsync()` and `copyDeviceToHostAsync()` require an appropriate later stream synchronization.

Solver users should normally operate through `addMaterial()`, `addGeometry()`, `addSphere()`, `addLSParticle()`, `addSPHParticle()`, and `addBond()` rather than mutating internal containers directly. Returned indices refer to the solver-owned host containers.

All public solver container accessors return const references. Use the solver's
`add...`, motion, reset and removal interfaces rather than swapping, reordering,
replacing or erasing internal elements. Insertion indices are stable under
append and neighbor search, not under explicit removal. Particle removal returns
the compaction map and updates native interaction endpoints; application-owned
indices must be remapped by the caller (see §11.4).

Do not retain C++ references or iterators into a host vector across an append: vector reallocation may invalidate them even though the returned numerical indices remain stable.

`LSGeometryContainer` stores reusable geometry in four arrays:

1. Geometry descriptors.
2. Signed-distance grid nodes.
3. Surface nodes.
4. Surface triangles.

Each host-side `LSParticle` also keeps the geometry arrays needed by contiguous CPU contact access. On the GPU, the particle stores a geometry index and kernels read the central device geometry container, avoiding repeated device geometry.

## 14. SPH–DEM

`SPHDEM` uses uniform fluid properties stored by the solver rather than repeated in every SPH particle:

```cpp
SPHDEM simulation(0);
simulation.setSPHProperties(
    spacing,
    smoothingLength,
    referenceDensity,
    dynamicViscosity);
simulation.setSPHMaximumVelocity(maximumExpectedVelocity);
```

`setSPHProperties()` must be called before adding SPH particles. It calculates the shared lattice-kernel normalization and viscous time scale. Particle mass is assigned as `referenceDensity * spacing^3`, and initial density is set to the reference density when a particle enters the solver.

### 14.1 Individual particles

```cpp
SPHParticle fluidParticle{{0.1, 0.1, 0.1}, {0.0, 0.0, 0.0}};
simulation.addSPHParticle(fluidParticle);
```

### 14.2 Rectangular block

`blockMinimum` is the lower corner of the block's lattice domain. Particles are placed at cell centers:

```cpp
simulation.addSPHBlock(
    {0.0, 0.0, 0.0},
    int3{80, 40, 40},
    {0.0, 0.0, 0.0});
```

### 14.3 Cylindrical jets

Use the `SPHJet` value type for a finite-duration inlet:

```cpp
SPHJet jet{
    outletCenter,
    jetDirection,
    jetRadius,
    jetSpeed,
    jetDuration};
const int firstJetParticle = simulation.addSPHJet(jet);
```

The value type normalizes its direction. This is the jet-creation interface: it creates a short upstream inlet buffer and continuously appends actual particles while the finite source runs; it does not preload the whole `jetSpeed * jetDuration` column. For spacing `s` and smoothing length `h`, the total axial layer count is `floor(jetSpeed * jetDuration / s)`, while the initial buffer contains at most `ceil(2h/s)+1` layers. Each layer uses the same circularly clipped regular lattice. Later layers replenish the upstream end without recycling particles that have already left the outlet.

New layers are scheduled every `s / jetSpeed` and appended only at common acoustic starts, before the interval's start snapshot. The acoustic plan is shortened to the first DEM node at or after a birth; insertion is therefore at most one DEM step late, and the new layer position compensates for the delay. No particle is born during a display/output observation. Stable global particle indices are append-only. With several jets or ordinary blocks, one source's particles can occupy several disjoint ranges:

```cpp
const std::vector<int>& emittedIndices = simulation.SPHJetParticleIndices(firstJetParticle);
```

The returned list contains that source's seeds and particles actually appended so far, not future particles. Copy it if it must remain available across further source growth. Do not infer its membership from `firstJetParticle + localIndex` after emission starts.

At each SPH acoustic substep, active buffer particles receive inlet velocity, reference density, and equation-of-state pressure before optional advection preparation and again after velocity integration. They skip their own continuity-equation update while contributing to neighbors' interactions. Crossing the outlet permanently releases a particle; it is not constrained again if the fluid later carries it back into the pipe. The finite layer count prevents an extra trailing reservoir after the requested duration. Once the final jet ends, the completion fast path skips later source work. `SPHJetsCompleted()` reports completion at the solver's visible time, even while the leading internal fluid state represents a later acoustic endpoint. CPU, GPU, and Hybrid use the same source timing and preserve emitted ownership through output, pause, and continuation. GPU capacity growth preserves already integrated device particles and initializes only the new tail.

Each finite source also opens its short inlet cylinder in the **SPH virtual-boundary samples only**. Samples inside the buffer, expanded radially and axially by half a particle spacing, receive zero quadrature volume. The opening covers the outlet plane but ends after the half-spacing downstream clearance, so a downstream target wall is not removed. Several sources form a union of these finite openings. Sample indices, owner mappings, level-set grids, and DEM contacts are unchanged; this is not a geometric hole or a modification of the LS signed-distance field.

Openings remain after a finite source stops emitting. Moving wall samples are tested again at acoustic starts; samples leaving an opening regain their normal volume. A changed mask refreshes the neighborhood and cached viscous loads before fluid integration. CPU, GPU, and Hybrid apply the same rule. Current-time observations rebuild only the visible host mask and do not advance the physical source or coupling state.

### 14.4 Level-set walls

SPH boundaries use fixed or movable `LSParticle` owners. `SPHDEM` samples internal virtual particles on the same uniform lattice and spacing as the fluid, then interpolates the level-set grid to retain solid-side samples within the kernel boundary band. Users do not add virtual particles directly.

An enclosing box reverses the signed-distance convention:

```cpp
levelset::BoxWall tank{{length, width, height}};
tank.buildLSGrid(spacing, 3, true);
tank.reverseSDFSign();
const int tankGeometry = simulation.addGeometry(tank);
```

### 14.5 SPH time-step schedule

The DEM time step is the minimum scheduling unit in `SPHDEM`:

- The SPH acoustic step is an integer multiple of the DEM step.
- It never exceeds the acoustic or advection/viscous stability limit calculated from smoothing length, sound speed, velocity scale, density, and viscosity.
- It is never smaller than one DEM step.
- The SPH advection step is an integer multiple of the SPH acoustic step and does not exceed its advection/viscous limit.
- Minimum smoothing length and the shared viscous time scale are stored once; the solver does not loop over all particles merely to recover uniform parameters.

If one DEM step is already larger than either SPH limit, initialization throws instead of silently violating the limit. The smoothing length must remain between one and two particle spacings, the configured sound speed must be at least ten times the configured velocity scale, and a level-set boundary grid must not be coarser than the SPH smoothing length. Background grids or compact neighbors normally rebuild once per SPH advection step. Actual fluid and wall displacement is checked against the search buffer; exhausting that buffer triggers an earlier search-only rebuild without resetting density or the advection phase. The physical kernel remains truncated at `2 * smoothingLength`.

Current values are available through `SPHTimeStep()`, `SPHTimeStepLimit()`, `SPHAdvectionTimeStep()`, and their interval getters.

SPH uses DKD (position half-step, velocity full-step, position half-step); solids use DEM KDK (velocity half-step, position/orientation full-step, velocity half-step). At a common acoustic boundary the fluid advances first over `H = N * DEM_dt`. Its pressure and viscous reaction force and world torque are installed before the first solid half kick and held through all `N` solid substeps. Other solid forces continue to update every DEM step. No endpoint impulse correction is required. Wall geometry uses the current solid pose; wall velocity/acceleration use the completed solid interval's averages, with instantaneous kinematics at startup. This explicit coupling does not guarantee exact energy conservation or unconditional fluid-solid stability.

For current-time host observations between common acoustic endpoints, use `observeCurrentState(callback)` rather than reading leading integration buffers directly. Interior fluid positions use the DKD quadratic through the stored endpoint positions and velocities; velocity and density use linear interpolation, and pressure follows the equation of state. These interior values are a dense-output approximation, not an independent SPH solve at the observation time. At a common endpoint the committed fluid solution is used. Observation leaves the force evaluation schedule, device state, and solid histories unchanged.

`solve()` returns current-time host data even if it stops inside an acoustic interval. An unchanged subsequent solve on the same solver object continues the same interval, while append-only model changes or configuration changes restart from the synchronized visible state. Loading only a saved visible frame into a new solver is a restart, not restoration of the private leading endpoint and multirate phase; its future trajectory is not guaranteed to be bitwise identical to uninterrupted integration. See [SPH state and coupling corrections](SPH_CORRECTIONS.md) for the full observation contract, held-load convention, and focused regression coverage.

### 14.6 Optional periodic SPH neighbor images

Configure the translational cell independently of the integration coordinates:

```cpp
math::periodicDomain domain;
domain.minimum_ = {0.0, 0.0, 0.0};
domain.length_ = {1.0, 1.0, 1.0};
domain.periodicAxes_[0] = true;
simulation.setSPHPeriodicDomain(domain);
```

Each enabled length must be finite, positive and strictly larger than `4*h`,
the diameter of the physical SPH kernel support. CPU and CUDA searches visit
wrapped cells without duplicate neighbor contributions; physical interaction
distances and inlet openings use the nearest translated image. Integrated
particle positions remain continuous rather than being wrapped in place.
`math::periodicDomain` also exposes checked metadata, `wrapPosition()`,
`minimumImage()` and `imagePosition()` for frontend use.

This does not automatically periodicize rigid-body contacts or particle
integration. A derived solver/frontend remains responsible for solid wrapping,
cross-boundary DEM contacts and admissible interaction ranges. The separate
ON-DEM validation runner has its own CPU periodic contact law and is not a
general implementation of these boundary policies.

## 15. Performance and scaling

### 15.1 Choose the execution mode once

| Solver | CPU | GPU | Hybrid |
| --- | --- | --- | --- |
| `LSDEM` | LS on CPU | LS on GPU | Not accepted |
| `SphereDEM` | Sphere and LS on CPU | Sphere and LS on GPU | Sphere on GPU, LS on CPU |
| `SPHDEM` | SPH and LS on CPU | SPH and LS on GPU | SPH on GPU, LS on CPU |

Select the mode before initialization. Mode dispatch is cached by the base solver, so per-step code follows the selected path without repeatedly rediscovering the execution plan.

### 15.2 Build configuration

- Use `Release` for timing and throughput comparisons.
- Keep `FUNDEM_ENABLE_CUDA_DEVICE_DEBUG=OFF` outside device debugging.
- Request OpenMP for CPU search and interaction loops; check CMake confirms compiler/runtime support for the required reductions.
- Compile only the CUDA architectures that must run the binary to reduce build time and binary size.
- Disable tutorials, examples, validation executables and tests in a production consumer build when they are not needed. Set `FUNDEM_BUILD_EXAMPLES=OFF` to omit `fundem_example`.

### 15.3 Search and geometry resolution

Rigid-particle broad-phase cost depends on populated grid cells and candidate count. Configure boundaries tightly enough to avoid a mostly empty search domain, but leave sufficient room for the entire trajectory. The CUDA grid uses up to three finite-size levels plus a final infinite-mass level, avoiding a single cell size dictated by the largest particle.

Level-set resolution controls both memory and narrow-phase work. Refine the signed-distance grid and surface only as far as contact accuracy requires. Reuse one geometry descriptor for particles with the same shape; host LS particles intentionally retain contiguous geometry data for fast CPU access, while GPU kernels use the centralized geometry container.

SPH neighbor construction normally occurs once per advection interval. Acoustic substeps reuse the current neighborhood while its actual-displacement bound is valid; unusually rapid motion can require an earlier rebuild. A very small smoothing length increases particle count; an unnecessarily large smoothing length increases neighbors per particle. Both can dominate runtime.

### 15.4 Output and synchronization

VTU serialization and device-to-host snapshots can dominate short simulations. Prefer binary VTU, select only fields needed for analysis, and use a physically meaningful output interval. Each frame is a deliberate synchronization point; avoid manual `copyDeviceToHost()` calls inside the time loop unless host-side inspection is required.

The frame log reports allocated FunDEM device memory. This is container capacity, not total process memory reported by the CUDA driver. Capacity can exceed active size after growth so that repeated contact and neighbor allocation does not reallocate every step.

### 15.5 Stability before speed

A larger DEM step is useful only while resolving the stiffest contact or bond mode. For SPH, the solver separately evaluates acoustic and advection limits and quantizes both to DEM-step multiples. Do not raise the design velocity merely to obtain a different sound speed without checking density variation and Mach number. Invalid SPH state statistics indicate a physical or numerical failure and should be diagnosed instead of hidden by reducing output frequency.

## 16. Source tree and API documentation

### 16.1 Source-tree guide

[`src/README.md`](../src/README.md) introduces the eight library modules and the test suite. The root `CMakeLists.txt` configures the project, and `docs/`, `cmake/`, `tutorial/`, `example/`, and `validation/` remain at the repository root.

| Directory | Responsibility |
| --- | --- |
| [`src/math/`](../src/math/README.md) | Shared CPU/GPU scalars, vectors, matrices, quaternions, tolerances, and indexing |
| [`src/data/`](../src/data/README.md) | Host AoS/device SoA storage, VTU/DAT output, energy output, and level-set input objects |
| [`src/material/`](../src/material/README.md) | Unified bulk materials and symmetric contact-pair coefficients |
| [`src/geometry/`](../src/geometry/README.md) | Four-array reusable level-set geometry storage |
| [`src/particle/`](../src/particle/README.md) | Point masses, rigid bodies, spheres, level-set particles, SPH particles, virtual particles, and spatial grids |
| [`src/interaction/`](../src/interaction/README.md) | Contacts, bonds, neighbor ranges, histories, mapping, and CPU contact search |
| [`src/execution/`](../src/execution/README.md) | Shared numerical formulas and backend execution helpers |
| [`src/execution/cpu/`](../src/execution/cpu/README.md) | Host container loops and temporary assembly storage |
| [`src/execution/cuda/`](../src/execution/cuda/README.md) | CUDA kernels and host launch interfaces |
| [`src/solver/`](../src/solver/README.md) | Common lifecycle and `LSDEM`, `SphereDEM`, and `SPHDEM` |
| [`src/tests/`](../src/tests/README.md) | CPU/CUDA unit and state-consistency tests |
| [`tutorial/`](../tutorial/README.md) | Short CPU-default Sphere, LS contact, Bond and SPH tutorials |
| [`example/`](../example/README.md) | Simple native C++ case built with the core |
| [`validation/`](../validation/README.md) | LAMMPS/PFC reference cases and ONDEM benchmarks |
| [`cmake/`](../cmake/README.md) | Installed CMake package configuration |

All public code is in namespace `fundem`; geometry input objects are in `fundem::levelset`, mathematical types in `fundem::math`, shared physics functions in `fundem::execution`, and CUDA launch interfaces in `fundem::cuda`.

### 16.2 Documentation layers and Doxygen

The documentation is intentionally layered:

1. The root [`README`](../README.md) explains the project, configuration, and first build.
2. This guide documents solver construction, lifecycle, output, and backend behavior.
3. The source index and each module README under `src/` describe ownership and invariants.
4. Header Doxygen records class purpose, non-obvious member semantics, constraints, side effects, and synchronization requirements.

Obvious getters and setters are intentionally not annotated one by one when their behavior is completely expressed by the declaration. Operators and mechanical `DeviceLayout` field tags are treated the same way. This keeps generated reference pages focused on information that cannot be recovered from the identifier alone.

To generate reference pages locally, install Doxygen 1.9.8 or newer and Python 3.9 or newer. Python must be available as `python3` on `PATH`; no additional Python packages, CUDA toolkit, or compiled solver are needed. Run from `workplace/`:

```bash
doxygen --version
python3 --version
python3 FunDEMBeta/docs/buildApi.py --output build-api
```

The builder runs Doxygen from the core directory and requires a fresh output directory outside the source tree. Open `build-api/html/index.html` beside the source directory. Linked source files, numerical records, and existing media are copied into that external output; repository files remain unchanged. The command checks local file links and anchors, and fails on documentation errors or source changes during generation. Intentionally undocumented getters and setters remain allowed.

The Doxygen configuration explicitly parses `.cu` and `.cuh` as C++, includes the module READMEs, and uses GitHub-compatible heading IDs. Documentation-only input filters preserve source line numbers while adapting Markdown math and resolving container aliases for the API parser. They do not change compiled code or the original files. Equations use MathJax, which requires network access to its CDN when viewing the pages.

For source navigation in VS Code, use the same CMake build directory for compilation and `compile_commands.json`; see [Section 4.1](#41-vs-code-configuration). Public user construction belongs in `src/solver/*.h`, model state in `src/particle/`, `src/material/`, and `src/interaction/`, shared formulas in `src/execution/`, and CUDA launch contracts in `src/execution/cuda/`.

## 17. Troubleshooting

### CMake finds no CUDA compiler

Symptoms:

```text
CUDA compiler not found: the CPU modules will still be built.
```

Actions:

1. Run `nvcc --version` in the same shell used for CMake.
2. Check that the CUDA toolkit `bin` directory is on `PATH`.
3. Remove or use a new build directory after changing compiler/toolkit discovery.
4. Reconfigure with `-DFUNDEM_ENABLE_CUDA=ON`.
5. If `nvcc` is intentionally outside `PATH`, add `-DCMAKE_CUDA_COMPILER=/absolute/path/to/nvcc`.

### GPU mode selection throws that the backend is unavailable

If `setExecutionMode()` reports this error, the library was configured without CUDA or CUDA detection failed. Reconfigure and rebuild FunDEM with a working CUDA compiler. Linking a CPU-only FunDEM package cannot enable GPU at runtime.

### CUDA architecture build failure

Unset `FUNDEM_CUDA_ARCHITECTURES` to use the compiler default, or set it to an architecture supported by both the installed toolkit and the target GPU.

### VS Code cannot jump between declarations and definitions

The most common cause is a missing or stale compile database. Reconfigure CMake, point cpptools to the resulting `compile_commands.json`, reset its IntelliSense database, and verify that VS Code and the compiler run in the same WSL/Linux environment.

### `Cannot add an invalid ...`

- Material: set positive Young's modulus and density, valid Poisson ratio and restitution.
- Sphere: assign positive radius and a material from this solver.
- LSParticle: assign a unified material and geometry from this solver.
- Missing LS-LS stiffness: assign both nonnegative per-area stiffnesses on each required material pair; there is no LS-stiffness activation switch.
- Geometry: build both its surface representation and signed-distance grid as required by the selected shape.
- SPH: call `setSPHProperties()` before adding particles.
- Bond: set equivalent length, connection, and stiffness; use this solver's particle containers.

### Output appears in an unexpected directory

Relative output directories are based on the executable location. Pass an absolute path to `setOutputDirectory()` when output must go to a specific working directory.

The example reserves a fresh run directory, so earlier runs are retained.

### First solve removed old results

This is intentional. The first `solve()` removes `.vtu` and `.dat` files below the configured output directory. Use a new output directory or move results that must be preserved before starting a new solver instance.

### CPU build is unexpectedly serial

Check the CMake configuration log for:

```text
OpenMP not found: host loops will use serial execution.
```

Install a compiler/runtime supporting the required OpenMP minimum/maximum reductions and configure again with `FUNDEM_ENABLE_OPENMP=ON`. A found but insufficient implementation also falls back to serial execution. `FUNDEM_ENABLE_OPENMP=OFF` explicitly selects serial loops.

### A short tutorial seems to produce only the initial frame

Each tutorial sets its output-step interval in the source. A diagnostic run
shorter than that interval can write only the initial frame. Increase the count
passed to `solve(...)` or reduce `setOutputStepInterval()`, then rebuild.
