# `example/`

> **Purpose:** Store larger validation and application examples built when `FUNDEM_BUILD_EXAMPLES=ON` (the default).

## Core Contents

| Path | Scenario | Backend |
| --- | --- | --- |
| [`structureIceLoad/`](structureIceLoad/README.md) | Bonded level ice interacting with a fixed LS structure | CPU, GPU, or Hybrid |
| `zhang2017ObstacleRWCSPH.cpp` | Zhang 2017 weakly compressible SPH dam break impacting a fixed obstacle | SPH on CUDA, LS boundaries on CPU |
| `canelas2016StackedCubesRWCSPH.cpp` | Canelas 2016 dam break impacting three vertically stacked PVC cubes | SPH on CUDA, movable LS particles on CPU |

These examples demonstrate application-specific solver derivation, bonded level ice, SPH water blocks, LS structures, fluid-solid coupling, binary VTU output, and CUDA device selection.

## Key Rules

- Examples demonstrate interfaces and data flow; they do not replace mesh convergence, parameter calibration, or quantitative literature validation.
- The structure-ice example is available in CPU-only builds; CMake skips only the SPH examples when CUDA is unavailable.
- Case-specific solver classes and kernels remain below `example/` and are not installed as core FunDEM APIs.
- Relative output directories are resolved from the executable directory.
- VTU output is binary by default; `--ascii` selects ASCII.
- The Canelas case uses roughly 645,000 fluid particles at its reference spacing and is intentionally much larger than a tutorial.

## Build

Open a terminal in `workplace/`, the directory containing `FunDEMBeta/`, and keep it there for all commands.

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_ENABLE_CUDA=ON -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build --target structureIceLoad zhang2017ObstacleRWCSPH canelas2016StackedCubesRWCSPH -j
```

A CPU-only structure-ice build is also supported:

```bash
cmake -S FunDEMBeta -B build-cpu -DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build-cpu --target structureIceLoad -j
```

## Run

```bash
build/example/structureIceLoad/structureIceLoad \
    --mode hybrid \
    --steps 1000 \
    --device 0 \
    --output output/structureIceLoad
```

```bash
build/example/zhang2017ObstacleRWCSPH \
    --steps 1000 \
    --device 0 \
    --output output/zhang2017
```

```bash
build/example/canelas2016StackedCubesRWCSPH \
    --steps 1000 \
    --device 0 \
    --output output/canelas2016
```

- `--steps N`: total number of DEM steps; omit it to use the example's full default duration.
- `--mode cpu|gpu|hybrid`: structure-ice execution mode; its default is `cpu`.
- `--device N`: non-negative CUDA device index.
- `--output DIR`: output directory.
- `--ascii`: write ASCII VTU instead of the default binary format.

## Case Construction Patterns

The structure-ice example derives `structureIceLoadSolver` from `SphereDEM`, adds a triangular bonded-sphere ice sheet, and represents each structure component with one LSParticle. Its specialized hydrodynamics and load accumulator are intentionally local to that example.

The two SPH examples follow this sequence:

1. Parse a small set of runtime arguments.
2. Define the physical scale, material properties, spacing, smoothing length, and DEM time step.
3. Create an `SPHDEM` solver and add fixed or movable LS boundaries.
4. Add the SPH water block.
5. Select `executionMode::Hybrid`, keeping LS motion and solid-solid contact on CPU while SPH runs on CUDA.
6. Configure binary VTU output and call `solve()` for the requested duration.

## Scenario Notes

### Structure–ice load

- Keeps the derived solver, level-ice generator, hydrodynamic formulas, CUDA kernels, and executable in one directory.
- Uses exact spherical-cap buoyancy and quadratic translational/rotational water drag.
- Supports `CPU`, `GPU`, and `Hybrid`; Hybrid runs ice spheres on the GPU and LS structures on the CPU.
- Writes one interval-averaged force vector per LS structure component to `structureIceLoad.dat` in MN.
- See the [case-specific README](structureIceLoad/README.md) for its file ownership and parameter flow.

### Zhang 2017 obstacle impact

- Uses fixed LS tank/obstacle geometry.
- Demonstrates the cheaper Hybrid wall-coupling path where only virtual-boundary kinematics and reactions cross the host/device boundary.
- Uses a fixed `1e-5 s` DEM step and writes at `0.05 s` intervals.

### Canelas 2016 stacked cubes

- Adds three movable LS cubes above the downstream floor.
- Derives the DEM time step from the selected contact frequency.
- Exercises SPH pressure/viscosity, fluid-solid reaction, LS rigid-body motion, and cube contact together.
- Uses a much larger fluid count and should be treated as a performance/validation workload rather than a smoke test.

## Reproducible Runs

Record the executable revision, CUDA toolkit, GPU model, selected device, command line, particle spacing, time step, and output format alongside results. `--steps` is useful for verifying initialization and the first few frames, but it does not reproduce the published physical duration unless the full default step count is used.

For performance measurements, use a Release build, binary output, a dedicated output drive when possible, and enough steps to amortize initialization. Measure simulation time separately from VTU writing, and report particle counts with the timing.

## Output Inspection

For the structure-ice example, open sphere and LSParticle VTU sequences together and inspect `structureIceLoad.dat`. For SPH cases, open SPH particle and LSParticle sequences together in ParaView. Check the first frame for sign convention, geometry placement, and unexpected particles outside the domain before starting a full run. `energy.dat` contains only solid and interaction energy; it is not the total fluid-solid energy balance for SPH cases.

## Extension Guide

New examples should follow the same command-line style and declare their backend requirements in CMake. Keep case-specific derived solvers and kernels beside their executable. Short introductory workflows belong in `tutorial/`.
