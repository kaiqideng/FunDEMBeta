# `example/structureIceLoad/`

> **Purpose:** Keep the application-specific level-ice/structure solver, its physical extensions, and its executable outside the core FunDEM solver library.

## Ownership

`structureIceLoadSolver` derives from the public `SphereDEM` solver. It is intentionally an example-local class: configuring FunDEM with `FUNDEM_BUILD_EXAMPLES=OFF` removes this entire feature without changing the installed library or the `solver/solvers.h` umbrella header.

| File | Responsibility |
| --- | --- |
| `structureIceLoad.cpp` | Builds and runs the example case |
| `structureIceLoadSolver.h/.cpp` | Defines the `SphereDEM` specialization, level-ice generator, water state, and structure-force output |
| `structureIceLoadFunctions.h` | CPU/GPU-shared spherical submergence, buoyancy, translational drag, and rotational drag formulas |
| `structureIceLoadKernel.cuh/.cu` | Streamed CUDA launchers for ice hydrodynamics and LS-structure force accumulation |
| `CMakeLists.txt` | Builds one CPU-only or CUDA-enabled example executable |

## Physical Representation

- Every ordinary sphere is one level-ice element.
- Unique nearest-neighbor edges in the staggered triangular lattice become sphere-sphere bonds.
- Selected boundary ice elements are assigned infinite mass to anchor the sheet.
- Every LSParticle is one structure component.
- The fluid is represented by a prescribed uniform water velocity, density, free-surface elevation, and drag coefficient; this is not an SPH fluid solve.
- Buoyancy uses the exact submerged spherical-cap volume.
- Translational and rotational drag retain the original IceDEM quadratic model.
- `structureIceLoad.dat` reports the resultant internal force already assembled on each LSParticle before any user LS external load is added.

## Execution Modes

| Mode | Ice spheres | LS structures | Load accumulator |
| --- | --- | --- | --- |
| `CPU` | CPU | CPU | Host |
| `GPU` | GPU | GPU | Device |
| `Hybrid` | GPU | CPU | Host |

All CUDA launchers receive the solver stream explicitly. The example-local CUDA source is compiled only when `FUNDEM_HAS_CUDA` is true.

## Build

Open a terminal in `workplace/`, the directory containing `FunDEMBeta/`, and keep it there for all commands:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_EXAMPLES=ON -DFUNDEM_ENABLE_CUDA=ON
cmake --build build --target structureIceLoad -j
```

For a CPU-only build:

```bash
cmake -S FunDEMBeta -B build-cpu -DFUNDEM_BUILD_EXAMPLES=ON -DFUNDEM_ENABLE_CUDA=OFF
cmake --build build-cpu --target structureIceLoad -j
```

## Run

```bash
build/example/structureIceLoad/structureIceLoad \
    --mode hybrid \
    --device 0 \
    --output output/structureIceLoad
```

Use `--mode cpu`, `--mode gpu`, or `--mode hybrid`. `--steps N` limits a diagnostic run to exactly `N` DEM steps; omitting it runs the configured two-second case. Binary VTU is the default, while `--ascii` selects readable XML arrays.

## Construction Flow

1. Configure water with `setWaterCondition()`.
2. Set ice thickness and density through `simulation.ice()`.
3. Call `levelIce::build()` to generate buoyancy-equilibrium sphere centers and unique bond edges.
4. Add all spheres before creating bonds so generated indices remain stable.
5. Convert the physical ice elastic/fracture parameters with the `levelIce` stiffness helpers.
6. Add one or more LS geometries and LSParticles for the structure.
7. Configure the global boundary, gravity, DEM time step, and output schedule.
8. Call `solve()`.

Changing ice thickness or density clears previously generated points and connections. Rebuild the lattice before adding it to a solver.

## Output

The inherited solver writes:

- sphere VTU files below `particle/`;
- movable and fixed LSParticle VTU files below `particle/`;
- contact and bond VTU files below `interaction/`;
- `energy.dat` for the solid DEM system.

The derived solver additionally writes `structureIceLoad.dat`. It contains time followed by `Fx`, `Fy`, and `Fz` for every LS structure component. Forces are averaged over the DEM steps since the preceding output frame and converted to MN. Frame zero contains zero accumulated load because it is written before time integration.

## Adaptation Notes

The example assumes all spheres are ice elements and all LSParticles are structure components. If a future case mixes unrelated spheres or LS objects, introduce explicit case-local membership arrays rather than adding ice-specific flags to the core particle containers.
