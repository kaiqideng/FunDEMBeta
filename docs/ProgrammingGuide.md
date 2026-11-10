# Write your own FunDEM case

[User guide](UserGuide.md) · [Minimal source](../example/main.cpp)

Use the three-file [`example/`](../example/README.md) as your starting point.
Its `main.cpp` contains the whole model; there is no argument parser or hidden
setup function. This guide uses the public C++17 API and SI units.

## 1. Build and run the starting case

Keep the terminal in the directory containing `FunDEMBeta/`. Put builds and
results outside the source tree:

```text
workplace/
├── FunDEMBeta/
│   └── example/       # main.cpp, CMakeLists.txt, README.md
└── build/            # generated files
```

```bash
cmake -S FunDEMBeta -B build -DCMAKE_BUILD_TYPE=Release -DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build --target fundem_example --parallel
OMP_NUM_THREADS=2 build/example/fundem_example
```

On Windows, run CMake in a Developer PowerShell with your C++ toolchain.
For MSVC, add `-DFUNDEM_ENABLE_OPENMP=OFF` to the CPU configuration above;
this uses serial loops because classic MSVC OpenMP lacks the required grid reductions.
For a multi-configuration generator, build with `--config Release`, then run:

```powershell
$env:OMP_NUM_THREADS = "2"
.\build\example\Release\fundem_example.exe
```

For a single-configuration Windows generator, the executable is directly under
`build/example/`. The small thread count avoids excessive overhead for two
spheres; it is not a requirement for larger cases.

The unmodified case collides two spheres with restitution 0.8. It runs 6000
steps at 1e-6 s and saves every 100 steps, producing 61 sphere frames including
the initial state. The console prints the actual output directory.

## 2. Know what to edit

Read `main.cpp` from top to bottom:

| Part | What you control |
| --- | --- |
| Constants | Size, density, elastic properties, restitution, friction and run settings |
| `SphereDEM` constructor | CPU, GPU or Hybrid; optional device index |
| Material and pair setup | Bulk properties and contact law |
| Particle loop | Number, shape, position and velocity of bodies |
| Boundary/gravity/time step | Search region and integration settings |
| Fresh output directory | Separate results for each invocation |
| Output calls and `solve(steps)` | Saved fields, cadence and number of additional steps |

For example, change `restitution` to `0.9` and rebuild the same target. Source
changes take effect only after recompilation. The template and all four
tutorials keep model and run settings in `main()` and take no command-line
options. Each uses `solver/solvers.h` as its single FunDEM include; necessary
standard C++ headers remain explicit.

Material construction is **E, Poisson ratio, restitution, density**:

```cpp
const int materialIndex = simulation.addMaterial(
    material{young, poisson, restitution, density});
```

Lengths are metres, time seconds, density kg/m^3, Young's modulus Pa, velocity
m/s and angular velocity rad/s. Restitution and friction are dimensionless.
Keep density and E positive, restitution in [0,1], and Poisson ratio between
-1 and 0.5. Reducing particle size or increasing stiffness generally requires
reconsidering the time step. Check convergence by halving it while keeping
the physical duration and output times comparable.

## 3. Keep your own case separate

Copy `FunDEMBeta/example/` to `FunDEMBeta/my_case/` using your file manager.
Keep `main.cpp`; change the copied `CMakeLists.txt` to:

```cmake
add_executable(my_case main.cpp)
target_link_libraries(my_case PRIVATE FunDEM::FunDEM)
```

In the core's root `FunDEMBeta/CMakeLists.txt`, extend the existing examples
block; do not add a second copy of that block:

```cmake
if(FUNDEM_BUILD_EXAMPLES)
    add_subdirectory(example)
    add_subdirectory(my_case)
endif()
```

Reconfigure after adding the directory/target, then build only your case:

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build --target my_case --parallel
OMP_NUM_THREADS=2 build/my_case/my_case
```

Use a unique target name for each case. The two-line case CMake file is a
**subdirectory of the core**, so do not configure it with `cmake -S my_case`.
For a separately installed library and standalone application, use the full
[external-project CMake recipe](UserGuide.md#43-external-application) instead.

## 4. Add particles and choose contacts

The snippets below go inside the existing `main()`; they are not independent
programs. Keep the existing output-folder and solve code.

To replace the original two-sphere loop with three resting spheres, use:

```cpp
for (int i = 0; i < 3; ++i)
{
    particle sphere;
    sphere.setRadius(radius);
    sphere.setMaterial(simulation.materials(), materialIndex);
    sphere.setPosition({(i - 1) * 2.2 * radius, 0, 0});
    sphere.setVelocity({0, 0, 0});
    simulation.addSphere(sphere);
}
```

Assign radius and material before insertion so mass and inertia are valid.
For initial non-overlapping packing, center spacing must exceed the sum of
radii. `addSphere()` returns a container index; retain it to address a body later.
Appending and neighbor searches do not change that index. Explicit removal
compacts the population and returns an old-to-new map; update your saved indices
from that map before using them again. Do not retain references into containers
while adding or removing objects.

The original material pair is for contacts between bodies of `materialIndex`.
For a second material, add it and configure the cross-pair explicitly:

```cpp
const int secondMaterial = simulation.addMaterial(material{2e9, 0.25, 0.6, 2700});
materialPair mixed{materialIndex, secondMaterial};
mixed.setFrictionCoefficients(0.3);
simulation.setMaterialPair(mixed);
```

Assign `secondMaterial` to the particles that should use it. Self-pairs and
cross-pairs are created automatically, but each pair's friction must be set
when you need a nonzero value. Pair coefficients are effective contact values.
Calling `setMaterialPair()` again replaces the complete pair, not just its
last changed field.

Leave `setLinearStiffness()` and `setDampingCoefficients()` out for automatic
Hertz-Mindlin stiffness and restitution-derived damping. Material restitution
is combined by harmonic mean. Explicit damping setters override all four
dashpots; passing zeros suppresses damping. Do not set them to zero when you
intend automatic damping. See [material-pair rules](UserGuide.md#61-materials).

## 5. Add a real wall

`setBoundary()` defines the search region; it does **not** create physical walls.
The existing `solver/solvers.h` include also supplies the level-set geometry
and particle types. Insert this block after the sphere loop and before
`setBoundary()`:

```cpp
levelset::PlaneWall plane{math::Vec3::unitZ(), 0.4};
plane.buildLSGrid(0.02, 2);
const int planeIndex = simulation.addGeometry(plane);

LSParticle wall;
wall.setMaterial(simulation.materials(), materialIndex);
wall.setGeometry(simulation.geometries(), planeIndex);
wall.setPosition({0, 0, -0.03});
simulation.addLSParticle(wall);
```

This uses the existing material and pair. The plane's normal points upward;
its reference plane is translated to z=-0.03 m. `PlaneWall` builds as fixed
geometry by default: zero volume gives infinite mass. Its default velocity
and spin are zero. A dense material alone does not make a fixed wall.
Sphere-wall Hertz contact needs no LS per-area stiffness setting.

To see the original two spheres reach and rebound from this floor, change
`steps` to `100000`, `outputEvery` to `1000`, keep dt=1e-6 s, and replace the
existing boundary/gravity calls with:

```cpp
simulation.setBoundary({-0.25, -0.25, -0.1}, {0.25, 0.25, 0.1});
simulation.setGravity({0, 0, -9.81});
```

This runs 0.1 s and schedules 101 frames. For a complete independent
sphere-on-plane example, read [`tutorial0.cpp`](../tutorial/tutorial0.cpp).
For movable non-spherical bodies, build surface nodes and the SDF grid, then
assign geometry and material as in [model construction](UserGuide.md#62-level-set-geometry).
LS-LS contact additionally requires explicit stiffness per unit area.

## 6. Select CPU, GPU or Hybrid

For GPU execution, replace the constructor in your source with:

```cpp
SphereDEM simulation(executionMode::GPU, 0);
```

Use `executionMode::CPU` or `executionMode::Hybrid` in place of `GPU` to
select either of those backends.

GPU/Hybrid also require configuring the core with `FUNDEM_ENABLE_CUDA=ON`,
a detected CUDA compiler, and a compatible device. Follow the
[CUDA build instructions](UserGuide.md#22-cuda-enabled-release-build), rebuild
your target, then run it. There is no silent fallback to CPU. `SphereDEM(0)`
only supplies device index 0 and still defaults to CPU.

Construct all objects on the host. `SphereDEM` supports all three modes;
Hybrid places spheres on GPU and LS bodies on CPU. Pure `LSDEM` supports CPU
or GPU, not Hybrid. User-defined GPU kernels need their own `.cu` source,
but selecting a built-in GPU solver does not require changing `main.cpp` to `.cu`.

## 7. Read results in your own code

Open `particle/sphere_*.vtu` as a time series in ParaView. The same fresh run
directory contains `energy.dat` and any contact/wall series that exist.
Relative paths are based on the executable location. Keep the template's
fresh-folder code: a solver's first `solve()`, including `solve(0)`, removes
existing VTU/DAT files inside its selected directory.

To inspect a body in C++, insert this **after** `simulation.solve(steps)`:

```cpp
math::Vec3 position, velocity;
simulation.observeCurrentState([&](const solver& base)
{
    const auto& state = static_cast<const SphereDEM&>(base);
    const auto& sphere = state.spheres().host().at(0);
    position = sphere.position();
    velocity = sphere.velocity();
});
std::cout << "x=" << position.x << ", vx=" << velocity.x << '\n';
```

The observer gives a coherent host view on CPU/GPU/Hybrid. Copy values out;
do not retain references or mutate/advance the solver inside the callback.
`synchronize()` only waits for device work; it does not copy state to host.
Use the copied values for your CSV, checks or control logic.

For periodic custom sampling, split a run into additional-step segments and
observe after each one. `solve(1000); solve(500);` advances 1500 total steps;
it does not restart time. Make supported model changes only before the first
solve or after a solve returns. See [solve lifecycle](UserGuide.md#112-solve-step-and-reconfiguration).

For an explicit completed-step velocity reset or particle removal, use the
[runtime mutation APIs](UserGuide.md#114-completed-step-resets-and-particle-removal)
instead of editing container storage. They preserve the remaining interaction
history; removal rejects particles that are still referenced by Bonds.

Scheduled output includes step zero but does not force a final frame when
the last step is off cadence. `writeOutput()` writes immediately and resets
the next scheduled step. `step()` alone does not trigger periodic output;
prefer `solve()` until you need an explicit stepping loop.

## 8. Extend the model

| Next task | Working reference |
| --- | --- |
| Non-spherical geometry and fixed walls | [User guide §6](UserGuide.md#6-model-construction), [tutorial1](../tutorial/tutorial1.cpp) |
| Bonds between bodies | [User guide §9](UserGuide.md#9-bonds), [tutorial2](../tutorial/tutorial2.cpp) |
| Add your own force/torque | [User guide §10](UserGuide.md#10-user-defined-external-force-and-torque) |
| Select extra VTU fields | [User guide §12](UserGuide.md#12-output-visualization-and-energy) |
| SPH setup and coupling | [User guide §14](UserGuide.md#14-sphdem), [tutorial3](../tutorial/tutorial3.cpp) |
| Numerical verification | [Validation cases](../validation/README.md) |

Start with a short CPU run, inspect the initial geometry and motion, then
increase duration or change backend. Keep the same physical setup when
comparing CPU/GPU and time-step convergence. A successful build or output
file is execution evidence; use a physical comparison to assess accuracy.
