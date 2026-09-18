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

If CMake cannot find a CUDA compiler, it reports that CUDA is unavailable and still produces the CPU modules; CUDA tests and GPU examples are then skipped.

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

This configuration never enables the CUDA language. OpenMP is used when CMake finds it; otherwise host loops remain serial.

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
| `FUNDEM_ENABLE_OPENMP` | `ON` | Use OpenMP for supported CPU loops when available |
| `FUNDEM_INSTALL` | `ON` | Generate installation and exported CMake package rules |
| `FUNDEM_BUILD_TUTORIALS` | `ON` | Build the CPU LSDEM tutorials |
| `FUNDEM_BUILD_EXAMPLES` | `ON` | Build the larger validation and application cases |
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
build/tutorial/tutorial0 --steps 1000
build/tutorial/tutorial1
build/tutorial/tutorial2 --steps 1000
build/tutorial/tutorial3 --steps 1000
```

Supported tutorial0 arguments are:

```text
--time-step DT  positive DEM time step
--duration T    positive simulation duration; defaults to 120 s
--steps N       non-negative DEM step count; overrides --duration
--output DIR    output directory
--ascii         write readable ASCII VTU instead of binary VTU
```

All bundled tutorials use CPU mode by default even in a CUDA-enabled build. Tutorial2 and tutorial3 can select another supported backend directly:

```bash
build/tutorial/tutorial2 --mode gpu --device 0
build/tutorial/tutorial3 --mode hybrid --device 0
```

Tutorial3 reproduces the classical three-dimensional bore-in-a-box configuration for `3.0 s`. In addition to particle VTU frames and `energy.dat`, `columnForce.dat` records the fixed square-column reaction uniformly at every `2e-5 s` DEM base step; the most recently computed SPH load is retained between acoustic updates. The example also writes `columnForceFiltered.dat` using a centered `13 ms` triangular moving average for inspecting the force history.

For a compact, literature-backed verification of the LS nodal contact
formulation, see the [two-sphere central-compression benchmark](LSDEMValidation.md).

### 3.2 Run the tests

```bash
ctest --test-dir build --output-on-failure
```

List tests without running them:

```bash
ctest --test-dir build -N
```

CPU unit tests are available whenever `BUILD_TESTING=ON`. CUDA builds additionally register host/device container and device state-consistency tests. Their sources are under [`src/tests/`](../src/tests/README.md), and a single-configuration build places their executables under `build/src/tests/`.

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
                    "${workspaceFolder}/FunDEMBeta/example",
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

This is the complete construction sequence used by the tutorial:

```cpp
#include "data/myLSObject.h"
#include "material/LSMaterial.h"
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
            simulation.addMaterial(LSMaterial{2.4e9, 7.1e8, 0.5, 0.55, 1200.0});
        const int wallMaterialIndex =
            simulation.addMaterial(LSMaterial{0.0, 0.0, 0.5, 0.55, 1.0});

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

        // Omit this line for the default CPU mode.
        simulation.setExecutionMode(executionMode::GPU);

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

The important ordering is:

1. Add the material and retain its returned index.
2. Fully build the level-set object before adding it as geometry.
3. Assign material and geometry from this solver's containers to the particle.
4. Add the valid particle to the solver.
5. Configure domain, time step, output, and optional GPU execution.
6. Call `solve()`.

## 6. Model construction

### 6.1 Materials

The ordinary sphere material constructor is:

```cpp
material sphereMaterial{
    normalStiffness,
    slidingStiffness,
    rollingStiffness,
    torsionalStiffness,
    slidingFriction,
    rollingFriction,
    torsionalFriction,
    restitution,
    density};
```

The level-set material constructor is:

```cpp
LSMaterial levelSetMaterial{
    normalStiffnessPerUnitArea,
    shearStiffnessPerUnitArea,
    slidingFriction,
    restitution,
    density};
```

`LSMaterial` adds no storage to `material`; it changes the public parameter names and records `materialType::levelSet`. It exposes sliding friction only. Rolling/torsional stiffness and friction interfaces are hidden, and their shared-storage values remain zero:

```cpp
LSMaterial wallMaterial{2.4e9, 7.1e8, 0.25, 0.55, 1200.0};
wallMaterial.setSlidingFrictionCoefficient(0.25);
```

Rules:

- Density must be finite and strictly positive.
- Stiffness and friction values must be finite and non-negative.
- Restitution must be within `[0, 1]`.
- `LSParticle::setMaterial()` accepts only a material marked as level-set.
- `LSDEM::addMaterial()` accepts `LSMaterial`; `SphereDEM` additionally accepts ordinary `material`.
- A sphere passed to `SphereDEM::addSphere()` must use ordinary `material`. `LSMaterial` is reserved for `LSParticle`, keeping contact-type selection unambiguous.

### 6.2 Level-set geometry

Geometry is built through `fundem::levelset::LSInfo` implementations:

| Type | Main construction interface |
| --- | --- |
| `Sphere` | `Sphere(radius)` then `buildSurfaceNode()` |
| `Superellipsoid` | Five shape parameters then `buildSurfaceNode()` |
| `TriangleMesh` | Vertex/triangle arrays or `loadOBJ()` |
| `PlaneWall` | Outward normal and finite surface size |
| `BoxWall` | Box size; use `reverseSDFSign()` for an enclosing tank |
| `BoxParticle` | Box size then `buildSurfaceNode(surfaceSpacing)` |
| `CylinderWall` | Bottom center, top center, and radius |
| `ConeWall` | Bottom/top centers and bottom/top radii |

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

### 6.3 Level-set particles

An `LSParticle` becomes valid only after a level-set material and geometry have both been assigned:

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

The material and geometry must belong to the same solver that receives the particle. The returned particle index is the stable host-container index used when creating bonds.

### 6.4 Spheres

Use `SphereDEM` and an ordinary `material`:

```cpp
SphereDEM simulation;

const int materialIndex = simulation.addMaterial(material{
    1.0e5, 5.0e4, 1.0e3, 1.0e3,
    0.3, 0.02, 0.01, 0.5, 2500.0});

particle sphere;
sphere.setPosition({0.0, 0.0, 0.2});
sphere.setRadius(0.025);
sphere.setMaterial(simulation.materials(), materialIndex);
const int sphereIndex = simulation.addSphere(sphere);
```

`particle::isValid()` requires a positive radius, a valid assigned material, and finite position, orientation, velocity, and angular velocity. Calling `setInfiniteMass()` after assigning radius and material creates a fixed sphere.

### 6.5 Domain, gravity, and device

```cpp
SphereDEM simulation(0);  // CUDA device index 0 if any particle type uses CUDA
simulation.setBoundary({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0});
simulation.setGravity({0.0, 0.0, -9.81});
simulation.setTimeStep(1.0e-6);
```

The minimum boundary must be strictly smaller than the maximum boundary in all three directions. The constructor device index is non-negative and is interpreted after any `CUDA_VISIBLE_DEVICES` remapping performed by the process environment.

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

All model objects are still appended on the host. Initialization creates the required device SoA storage and spatial grids. `solve()`, `writeOutput()`, `copyDeviceToHost()`, and `synchronize()` are the public synchronization boundaries for GPU execution.

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

The CPU path uses an OpenMP-parallel compact neighbor list built from a uniform background grid. The GPU path queries the device grid directly. In `Hybrid`, virtual-particle forces are transferred only when SPH advances, reduced to one resultant force and torque per LSParticle, and applied on the CPU at each DEM step; the full LSParticle container is not copied to the device every DEM step.

## 8. Contact models and search

### 8.1 Spatial grid

The CUDA rigid-particle search supports at most four grid levels:

- Finite-mass diameter ratio `<= 4`: one finite-mass level.
- Ratio `> 4` and `<= 16`: two finite-mass levels.
- Ratio `> 16`: three finite-mass levels.
- If infinite-mass particles exist, one level is reserved as the final level for them.

Particles are not reordered in the user container. The grid maintains a separate sorted particle-index array, so existing particle indices, contacts, and bonds remain stable. Infinite-mass particles are placed in the final grid level and are selected as the slave side of generated pairs. CPU contact search applies the same infinite-mass slave convention.

### 8.2 Contact types

- Sphere–sphere uses the complete spherical contact model.
- Sphere–LS uses the complete spherical model, including force, rolling torque, and torsional torque; the level-set particle is the slave.
- LS–LS uses surface-node area weighting.

Contact type is determined from the material marker rather than inferred from `effectiveRadius`.

For sphere–LS contact:

- Normal, sliding, rolling, and torsional stiffness come directly from the sphere's ordinary material.
- Rolling and torsional friction coefficients come directly from the sphere material; no LS-side coefficients are required.
- Sliding friction combines the two materials using the harmonic mean.
- Restitution combines the two materials.
- A non-level-set material cannot be assigned to an `LSParticle`.

This keeps sphere–sphere and sphere–LS force evaluation on one complete contact path while preserving the area-based LS–LS model.

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

The equivalent length must be finite and positive. All four stiffness inputs must be finite and non-negative. `setStiffness()` requires the equivalent length to have been set. Optional BK damage parameters are configured with:

```cpp
connection.setCrossSectionArea(area);
connection.setModeICriticalEnergy(modeI);
connection.setModeIICriticalEnergy(modeII);
connection.setModeMixityExponent(exponent);
connection.setDamageInitiationRatio(ratio);
connection.setMaximumEnergyReleaseRatio(ratioLimit);
```

`crossSectionArea` is the nominal bond cross-sectional area in square metres. Setting it to zero disables fracture; otherwise it converts elastic energy into energy release per unit area. `damageFactor` and `damageInitiationRatio` are bounded to `[0, 1]` and `(0, 1]`, respectively. Always check boolean setters when values originate outside the program.

## 10. User-defined external force and torque

Built-in gravity, contact, bond, and coupling forces are assembled first. Solver extension hooks run last, allowing a derived solver to add fields without replacing the standard assembly.

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

Time-step stability remains the case author's responsibility; use material stiffness, smallest mass, contact scale, and intended accuracy when choosing it.

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

Time, total step count, output frame count, contact history, and existing particles continue across calls. The next initialization resizes or rebuilds the affected working containers. Do not add objects or reset LSParticle motion while an explicitly initialized/stepping phase is active.

`step()` initializes on demand and leaves the solver initialized. Use it for controlled per-step workflows. `synchronize()` waits for the solver CUDA stream, and `copyDeviceToHost()` exposes the current device state when needed during an initialized GPU workflow.

## 12. Output, visualization, and energy

### 12.1 Configure output

```cpp
simulation.setOutputDirectory("case_files");
simulation.setOutputStepInterval(1000);
simulation.setVTUOutputFormat(vtuFormat::binary);
```

Binary VTU is the default and is normally preferred for speed and file size. Use `vtuFormat::ascii` only for readable diagnostics.

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
- Bond: `force`, master/slave torque and indices, `equivalentLength`, four elastic energies, global master/slave endpoint normal/tangent bases, `crossSectionArea`, damage factor, mode-I/mode-II critical energy, mode-mixity exponent, damage-initiation ratio, and maximum energy-release ratio.

The LSParticle `distanceToCenter` value is the distance from each output surface node to its owning particle center. Bond endpoint orientation fields are transformed to the current global frame before output.

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

All public solver container accessors return const references. Users may only append through the solver's `add...` interfaces; they must not swap, reorder, replace, or erase existing elements. Insertion indices stay stable because particles, contacts, bonds, histories, and spatial-grid index arrays all refer to them.

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

The legacy four-argument overload creates one static downstream particle column. The velocity defines both its axis and its initial particle velocity; `jetLength` controls only the initially populated geometry. This overload does not maintain an inlet velocity after insertion:

```cpp
simulation.addSPHJet(
    inletCenter,
    jetVelocity,
    jetRadius,
    jetLength);
```

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

The value type normalizes its direction, and this overload pre-fills an upstream virtual pipe of length `jetSpeed * jetDuration` once during setup. It does not inject particles or resize storage while solving. The solver records the jet's stable particle range and constrains only particles from that range that are still inside its pipe, so unrelated fluid and particles that have crossed the outlet remain unconstrained.

At each SPH acoustic substep, the inlet velocity is applied at the substep start before any advection-step preparation and again immediately after velocity integration. A jet active at the substep start remains active through that complete substep, so its end time is rounded upward by at most one acoustic substep. After the final jet ends, the completion fast path in `SPHJetState` bypasses all later pipe checks. `SPHJetsCompleted()` reports completion at the solver's visible time, even when internal fluid work is deferred. Output observations restore the internal completion flag along with the saved fluid state and do not end a jet prematurely.

### 14.4 Level-set walls

SPH walls use fixed `LSParticle` geometry. `SPHDEM` samples internal virtual particles on the same uniform lattice and spacing as the fluid, then interpolates the level-set grid to retain solid-side samples within the kernel boundary band. Users do not add virtual particles directly.

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

For current-time host observations between acoustic updates, use `observeCurrentState(callback)` rather than reading deferred particle containers directly. See [SPH state and coupling corrections](SPH_CORRECTIONS.md) for the observation contract, held-load convention, impulse correction, and focused regression coverage.

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
- Enable OpenMP for CPU search and interaction loops; check CMake output to confirm it was found.
- Compile only the CUDA architectures that must run the binary to reduce build time and binary size.
- Disable tutorials, examples, and tests in a production consumer build when they are not needed.

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

[`src/README.md`](../src/README.md) introduces the nine source modules. The root `CMakeLists.txt` configures the project, and `docs/`, `cmake/`, `tutorial/`, `example/`, and `validation/` remain at the repository root.

| Directory | Responsibility |
| --- | --- |
| [`src/math/`](../src/math/README.md) | Shared CPU/GPU scalars, vectors, matrices, quaternions, tolerances, and indexing |
| [`src/data/`](../src/data/README.md) | Host AoS/device SoA storage, VTU/DAT output, energy output, and level-set input objects |
| [`src/material/`](../src/material/README.md) | Ordinary and level-set material definitions |
| [`src/geometry/`](../src/geometry/README.md) | Four-array reusable level-set geometry storage |
| [`src/particle/`](../src/particle/README.md) | Point masses, rigid bodies, spheres, level-set particles, SPH particles, virtual particles, and spatial grids |
| [`src/interaction/`](../src/interaction/README.md) | Contacts, bonds, neighbor ranges, histories, mapping, and CPU contact search |
| [`src/execution/`](../src/execution/README.md) | Shared numerical formulas and backend execution helpers |
| [`src/execution/cpu/`](../src/execution/cpu/README.md) | Host container loops and temporary assembly storage |
| [`src/execution/cuda/`](../src/execution/cuda/README.md) | CUDA kernels and host launch interfaces |
| [`src/solver/`](../src/solver/README.md) | Common lifecycle and `LSDEM`, `SphereDEM`, and `SPHDEM` |
| [`src/tests/`](../src/tests/README.md) | CPU/CUDA unit and state-consistency tests |
| [`tutorial/`](../tutorial/README.md) | Compact CPU level-set tutorials |
| [`example/`](../example/README.md) | Larger application cases |
| [`validation/`](../validation/README.md) | Analytical, convergence, and backend-consistency cases |
| [`cmake/`](../cmake/README.md) | Installed CMake package configuration |

All public code is in namespace `fundem`; geometry input objects are in `fundem::levelset`, mathematical types in `fundem::math`, shared physics functions in `fundem::execution`, and CUDA launch interfaces in `fundem::cuda`.

### 16.2 Documentation layers and Doxygen

The documentation is intentionally layered:

1. The root [`README`](../README.md) explains the project, configuration, and first build.
2. This guide documents solver construction, lifecycle, output, and backend behavior.
3. The source index and each module README under `src/` describe ownership and invariants.
4. Header Doxygen records class purpose, non-obvious member semantics, constraints, side effects, and synchronization requirements.

Obvious getters and setters are intentionally not annotated one by one when their behavior is completely expressed by the declaration. Operators and mechanical `DeviceLayout` field tags are treated the same way. This keeps generated reference pages focused on information that cannot be recovered from the identifier alone.

To generate reference pages locally, install Doxygen 1.9.8 or newer and Python 3. Python must be available as `python3` on `PATH`; no additional Python packages, CUDA toolkit, or compiled solver are needed. Run from `workplace/`:

```bash
doxygen --version
python3 --version
cmake -E chdir FunDEMBeta doxygen Doxyfile
```

`cmake -E chdir` sets the working directory only for the Doxygen process, so the terminal remains in `workplace/`. Open `FunDEMBeta/docs/api/html/index.html`. Generated documentation is ignored by Git. Documentation errors fail the command rather than silently producing a successful build; intentionally undocumented getters and setters remain allowed.

The Doxygen configuration explicitly parses `.cu` and `.cuh` as C++, includes the module READMEs, and uses GitHub-compatible heading IDs. Documentation-only input filters preserve source line numbers while adapting Markdown math/media and resolving container aliases for the API parser. They do not change compiled code or the original files. Referenced images and the tutorial video are copied into the HTML output; equations use MathJax, which requires network access to its CDN when viewing the pages.

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

- Material: set a positive density and valid coefficients.
- Sphere: assign positive radius and a material from this solver.
- LSParticle: assign an `LSMaterial` and geometry from this solver.
- Geometry: build both its surface representation and signed-distance grid as required by the selected shape.
- SPH: call `setSPHProperties()` before adding particles.
- Bond: set equivalent length, connection, and stiffness; use this solver's particle containers.

### Output appears in an unexpected directory

Relative output directories are based on the executable location. Pass an absolute path to `setOutputDirectory()` when output must go to a specific working directory.

### First solve removed old results

This is intentional. The first `solve()` removes `.vtu` and `.dat` files below the configured output directory. Use a new output directory or move results that must be preserved before starting a new solver instance.

### CPU build is unexpectedly serial

Check the CMake configuration log for:

```text
OpenMP not found: host loops will use serial execution.
```

Install an OpenMP-capable compiler/runtime and configure again with `FUNDEM_ENABLE_OPENMP=ON`.

### A short tutorial seems to produce only the initial frame

The tutorial converts its `0.05 s` output interval to a step interval from the selected time step. A shorter smoke run writes frame zero and may finish before frame one. Increase `--steps` or reduce the output interval in the case source.
