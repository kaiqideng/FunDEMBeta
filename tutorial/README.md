# Tutorials

Four small, independent programs introduce the public FunDEM API. Read each
`main()` from top to bottom: material, geometry/particles, interactions, solver
settings, then `solve()`. Each source uses one FunDEM header,
`solver/solvers.h`; tutorial3 also includes the standard `<cmath>` header for
`std::sqrt`. Comments explain the main classes, function arguments and setup
order. No OBJ asset, custom solver subclass or tutorial-local CUDA kernel is
needed.

| Program | Lesson | Moving particles | Time step | Steps | Physical duration | Frames |
| --- | --- | --- | --- | --- | --- | --- |
| `tutorial0` | Sphere rebound on a fixed plane | 1 sphere | 1e-4 s | 3,000 | 0.3 s | 31 |
| `tutorial1` | Asymmetric LS pyramid and centroid-corrected contact | 1 LS particle | 1e-4 s | 3,000 | 0.3 s | 31 |
| `tutorial2` | Five-particle beam with a fixed root and four bonds | 4 spheres | 2e-5 s | 5,000 | 0.1 s | 21 |
| `tutorial3` | Small SPH block released inside a closed LS tank | 960 SPH particles | 2e-5 s | 10,000 | 0.2 s | 21 |

These are short teaching transients, not equilibrium solutions, experimental
benchmarks or accuracy certifications. All programs default to CPU execution.
The original large tutorial models have been replaced, not retained as examples.
For a reusable starting point for your own native case, see the
[`example/`](../example/README.md) C++17/CMake template. Its `fundem_example`
target builds through the core with `FUNDEM_BUILD_EXAMPLES=ON`. It uses the public core API and writes VTU and `energy.dat`;
the default model is a two-sphere Hertz collision. Follow the
[Programming Guide](../docs/ProgrammingGuide.md) to edit that example or copy it
and register a unique target in the core CMake build. The template's model and
run settings are edited in source.

## Build and run

Keep the terminal in `workplace/`, containing `FunDEMBeta/`:

```bash
cmake -S FunDEMBeta -B build -DCMAKE_BUILD_TYPE=Release -DFUNDEM_BUILD_TUTORIALS=ON
cmake --build build --target tutorial0 tutorial1 tutorial2 tutorial3 --parallel
build/tutorial/tutorial0
build/tutorial/tutorial1
build/tutorial/tutorial2
build/tutorial/tutorial3
```

With a multi-configuration Windows generator, build with `--config Release`
and use the executables in `build/tutorial/Release/`.

Small CPU jobs may be slower when many OpenMP threads are started repeatedly.
For repeatable timing on Linux/WSL, try `OMP_NUM_THREADS=2`; this is a runtime
setting, not a requirement or a universal best thread count.

## Edit the source

All model and run settings are visible in each `main()`. Change them in the
C++ source, rebuild the corresponding target, then run it without arguments.
There is no command-line parser or shared options helper.

| Setting | Where to edit |
| --- | --- |
| CPU, GPU or Hybrid; CUDA device | Solver constructor; tutorial1 uses `setExecutionMode()` |
| Material, contacts and geometry | The corresponding setup blocks |
| Number of steps | `solve(...)` |
| Time step and output cadence | `setTimeStep()` and `setOutputStepInterval()` |
| Output directory and VTU format | `setOutputDirectory()` and `setVTUOutputFormat()` |

For a ten-step diagnostic run, temporarily change `solve(...)` to `solve(10)`,
choose a fresh output directory in the source, and rebuild. Restore the case's
step count and rebuild for the complete run.

All programs select CPU and CUDA device 0 by default. GPU modes require a
CUDA-enabled build and compatible device. The pure `LSDEM` solver in tutorial1
supports CPU or GPU; tutorial0, tutorial2 and tutorial3 also support Hybrid.

## What to look for

- **tutorial0:** the sphere falls and rebounds; the Hertz-Mindlin response uses
  material elastic properties and restitution, plus the configured sliding
  friction. No linear stiffness override is needed.
- **tutorial1:** a tilted pyramid contacts the plane. `LSInfo` shifts its surface
  and SDF origin into one centroid frame. The pair supplies LS normal/shear
  coefficients in N/m^3, multiplied by nodal area during contact assembly.
- **tutorial2:** one infinite-mass root anchors four movable spheres. Beam
  stiffnesses come from a circular section; ordinary contact stiffness is zero
  to isolate bonds. Zero bond cross-sectional area disables fracture, not the
  explicitly supplied elastic stiffness. The transient is not a beam-convergence
  test.
- **tutorial3:** a regular, initially resting SPH block slumps and spreads.
  Reversing the closed box SDF places boundary solid outside the fluid cavity.
  Acoustic/advection intervals remain automatic integer multiples of the DEM
  step. No column-load recorder or post-run filter is included.

## Output

The default `tutorialN_files/` directory is beside its executable. On the first
`solve()`, the core removes old VTU/DAT files **inside the selected directory**,
so edit `setOutputDirectory()` to choose a fresh dedicated folder before
rebuilding when earlier results must be retained.

Only the basic binary VTU fields are selected. Open the matching time series in
ParaView, click **Apply**, reset the camera and press **Play**:

| Program | Main series under its output directory |
| --- | --- |
| tutorial0 | `particle/sphere_*.vtu`, `particle/fixedLSParticle_*.vtu` |
| tutorial1 | `particle/LSParticle_*.vtu`, `particle/fixedLSParticle_*.vtu` |
| tutorial2 | `particle/sphere_*.vtu`, `interaction/sphereBond_*.vtu` |
| tutorial3 | `particle/SPHParticle_*.vtu`, `particle/fixedLSParticle_*.vtu` |

The core also writes `energy.dat` at the output cadence. It contains the **solid
DEM subsystem** only, not SPH fluid energy. A run shorter than the output
interval may write just the initial frame; changing the step count does not
force an extra final frame.
