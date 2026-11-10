# FunDEMBeta

FunDEMBeta is a C++17 discrete-element framework for spherical and level-set particles, bonded and frictional contact, and SPH–DEM coupling. One host-side object definition supplies both CPU data and CUDA structure-of-arrays storage, so the same physical model can be executed by CPU, GPU, or supported hybrid solvers.

This README provides one complete path from a new development environment through the introductory tutorials and their results. Detailed interfaces and advanced workflows are documented in the [User Guide](docs/UserGuide.md). To write your own C++ case, follow the [Programming Guide](docs/ProgrammingGuide.md).

The CMake package version is **1.0.2**. Current interfaces include an expanded
native level-set primitive catalog, restoration of already prepared geometry
without another centroid shift, reusable CPU fixed-radius point queries,
completed-step particle resets/removal, and periodic SPH neighbor images. See
the [model-construction guide](docs/UserGuide.md#6-model-construction),
[runtime-edit contract](docs/UserGuide.md#114-completed-step-resets-and-particle-removal)
and [SPH guide](docs/UserGuide.md#14-sphdem) for their scope and limitations.

Using a Mac? Follow the dedicated [macOS quick start](README_macOS.md) for a CPU-only build on Apple Silicon or Intel.

Create a parent directory named `workplace` and open a terminal there before following this guide. Keep the terminal in `workplace/` for every command below; the source checkout and build directory will be `workplace/FunDEMBeta/` and `workplace/build/`.

## Program framework

```text
Bulk materials + material pairs + geometries + particles
                │
                ▼
        Host AoS containers ─────────► CPU / OpenMP execution
                │
                │ DeviceLayout fields
                ▼
        Device SoA containers ───────► CUDA execution
                │
                ▼
      Solver synchronization and output
                │
                ├── particle/*.vtu
                ├── interaction/*.vtu
                └── energy.dat
```

The solver hierarchy is:

```text
solver
└── LSDEM          level-set DEM
    ├── SphereDEM  spheres coupled to level-set particles
    └── SPHDEM     SPH coupled to level-set particles
```

The core library and its tests live under [`src/`](src/README.md). The repository root keeps the CMake entry point, documentation, build support, tutorials, examples, and validation cases:

| Directory | Responsibility |
| --- | --- |
| [`src/math/`](src/math/README.md) | CPU/CUDA mathematical types and functions |
| [`src/data/`](src/data/README.md) | Host/device containers and VTU/DAT output |
| [`src/material/`](src/material/README.md) | Unified bulk materials and symmetric contact-pair coefficients |
| [`src/geometry/`](src/geometry/README.md) | Reusable level-set geometry storage |
| [`src/particle/`](src/particle/README.md) | DEM, level-set, SPH, virtual particles, SPH source descriptors, and spatial grids |
| [`src/interaction/`](src/interaction/README.md) | Contacts, bonds, searches, SPH neighborhoods, and wall-coupling state |
| [`src/execution/`](src/execution/README.md) | Shared physical formulas, CPU container loops, and CUDA kernels |
| [`src/solver/`](src/solver/README.md) | Simulation lifecycle and execution modes |
| [`src/tests/`](src/tests/README.md) | CPU/CUDA regression and integration tests |
| [`docs/`](docs/UserGuide.md) | User guide, API documentation, and validation notes |
| [`cmake/`](cmake/README.md) | Installed CMake package configuration |
| [`tutorial/`](tutorial/README.md) | Small runnable demonstrations |
| [`example/`](example/README.md) | Simple native C++ case built with the core |
| [`validation/`](validation/README.md) | LAMMPS, PFC and ONDEM reference cases |

From `workplace/`, configure CMake with `-S FunDEMBeta -B build`. Public include names such as `#include "solver/solvers.h"` are unchanged; CMake supplies `src/` as the build include root and installs public headers under `include/FunDEM/<module>/`.

## 1. Configure the environment

The commands below assume a Linux or WSL2 environment. Required tools are:

- Git;
- CMake 3.24 or newer;
- a C++17 compiler such as GCC;
- a build tool such as Ninja;
- an NVIDIA driver and CUDA toolkit for the GPU backend;
- ParaView for viewing VTU results.

On Ubuntu or WSL2, install the basic CPU tools with:

```bash
sudo apt update
sudo apt install -y build-essential git cmake ninja-build
```

The optional energy-plot command at the end also uses:

```bash
sudo apt install -y python3-numpy python3-matplotlib
```

Install the NVIDIA driver and CUDA toolkit using the instructions for the operating system and GPU. The driver must support the installed CUDA toolkit. Confirm the complete toolchain before configuring FunDEM:

```bash
git --version
cmake --version
g++ --version
ninja --version
nvcc --version
nvidia-smi
```

If `cmake --version` reports a version older than 3.24, install a newer CMake before continuing. `nvcc` and `nvidia-smi` must both work when CPU and GPU support are to be compiled together.

## 2. Download the source code

With the terminal still in `workplace/`, clone the source:

```bash
git clone https://github.com/kaiqideng/FunDEMBeta.git
```

The resulting layout will be:

```text
workplace/      terminal working directory for this guide
├── FunDEMBeta/  source code
└── build/       generated in the next step
```

## 3. Build CPU and GPU support

Use one Release configuration requesting CUDA and OpenMP support. Tutorials and the
native C++ example are built, while validation executables and tests are omitted
from this first build:

```bash
cmake -S FunDEMBeta -B build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=ON \
    -DFUNDEM_ENABLE_OPENMP=ON \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DFUNDEM_BUILD_EXAMPLES=ON \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF

cmake --build build --parallel
```

CPU code is always compiled. When CUDA is found, the same configuration also compiles the CUDA libraries and kernels. During configuration, confirm that CMake reports an NVIDIA CUDA compiler. If it prints `CUDA compiler not found`, the build contains CPU support only; make `nvcc` available on `PATH`, then repeat the configure command with `cmake --fresh` to refresh the CMake cache without deleting saved simulation results.

OpenMP is enabled only when the compiler and runtime support the required
minimum/maximum reductions. Otherwise host loops use serial execution. Set
`FUNDEM_ENABLE_OPENMP=OFF` to request a serial build explicitly.

## 4. Run the tutorials

### Tutorial 0: sphere rebound

`tutorial0` drops one sphere onto a fixed level-set plane. It demonstrates material
properties, Hertz-Mindlin contact, gravity and output in 0.3 seconds of physical
time: 3,000 steps at `1e-4 s`, with one frame every `0.01 s`.

Run it from `workplace/`:

```bash
build/tutorial/tutorial0
```

### The other small tutorials

| Program | Main concept | Default size | Physical duration |
| --- | --- | --- | --- |
| `tutorial1` | Asymmetric LS geometry, centroid correction and contact | One tilted pyramid and a plane | 0.3 s |
| `tutorial2` | Elastic bonds with an infinite-mass root | Five spheres and four bonds | 0.1 s |
| `tutorial3` | SPH setup and a closed LS boundary | 960 fluid particles and a tank | 0.2 s |

All four default to CPU. The same sources can select GPU in a CUDA build;
tutorial0, tutorial2 and tutorial3 also support their solver's Hybrid mode.
Each program uses only public APIs and writes 21 or 31 frames. They are teaching
transients, not large demonstration runs or experimental benchmark reproductions.

```bash
build/tutorial/tutorial1
build/tutorial/tutorial2
build/tutorial/tutorial3
```

See the [tutorial guide](tutorial/README.md) for the source walkthrough and
settings. Each tutorial includes `solver/solvers.h` and defines the model and
run controls directly in `main()`. Edit the source and rebuild to change them;
run the executable without arguments. To start your own native C++ case, edit
[`example/main.cpp`](example/main.cpp). The complete setup is in one `main()`:
materials, particles, backend, time step, output and `solve()`.

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build --target fundem_example --parallel
build/example/fundem_example
```

The default is a two-sphere Hertz collision on CPU. Edit the model and run
settings in `main.cpp`, then rebuild. For another case, copy the directory,
give its target a unique name and register it in the core CMake build. Each run
writes native VTU and `energy.dat` in a fresh directory. See the
[example guide](example/README.md).

For the tutorials, the terminal reports the calculated step count, physical time, output-frame index, and allocated device memory at every frame. Their default results are written beside the executable:

```text
build/tutorial/tutorial0_files/
├── particle/
│   ├── sphere_000000.vtu
│   ├── sphere_000001.vtu
│   ├── fixedLSParticle_000000.vtu
│   └── ...
├── interaction/
│   ├── sphereLSParticleContact_*.vtu
│   └── ...
└── energy.dat
```

## 5. View the DEM result

Open ParaView and load these two file series:

```text
build/tutorial/tutorial0_files/particle/sphere_*.vtu
build/tutorial/tutorial0_files/particle/fixedLSParticle_*.vtu
```

ParaView normally recognizes the numbered files as time series. Select both readers in the Pipeline Browser, click **Apply**, press **Reset Camera**, and use **Play** to animate the rebound. Load `interaction/sphereLSParticleContact_*.vtu` as another series when contact points and normals are also required.

VTU is appended binary by default. ParaView reads it directly; no conversion is necessary.

## 6. View the DEM system energy

`energy.dat` is a whitespace-separated text table. Inspect its header and first frames with:

```bash
head -n 6 build/tutorial/tutorial0_files/energy.dat
```

Each row contains simulation time followed by solid-particle kinetic energy, gravitational potential energy, contact elastic energy, bond elastic energy, and `totalEnergy`. SPH-fluid energy is intentionally excluded from this DEM energy record.

To plot the important energy groups with the NumPy and Matplotlib packages installed in Section 1, run:

```bash
python3 - <<'PY'
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

path = Path("build/tutorial/tutorial0_files/energy.dat")
energy = np.genfromtxt(path, names=True)
elastic_names = [name for name in energy.dtype.names if name.endswith("ElasticEnergy")]
elastic = sum((energy[name] for name in elastic_names), np.zeros_like(energy["time"]))

plt.plot(energy["time"], energy["kineticEnergy"], label="kinetic")
plt.plot(energy["time"], energy["gravitationalPotentialEnergy"], label="gravitational")
plt.plot(energy["time"], elastic, label="contact + bond elastic")
plt.plot(energy["time"], energy["totalEnergy"], label="total", linewidth=2)
plt.xlabel("time [s]")
plt.ylabel("energy [J]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
PY
```

`totalEnergy` is the sum of the reported DEM mechanical-energy components. It may decrease because restitution, friction, and damping dissipate energy; dissipated energy is not stored as an additional state component.

## Verification and known limits

Core regression tests live in [src/tests](src/tests/README.md). The
[validation guide](validation/README.md) documents reference comparisons and
their physical-validation limits. The
[LAMMPS](validation/LAMMPSComparison/README.md) and
[PFC](validation/PFCComparison/README.md) categories include English Markdown
reports and dedicated C++ programs for 35 LAMMPS-derived variants and three
PFC-derived cases, including two explicitly failed Mindlin sliding-to-rolling
cases. Each program directly constructs its model using the public core API. Those reports retain their
original source pins, tolerances, dates and build limitations; they are not a
claim that every current core or desktop feature has been physically validated,
and they do not claim a local execution of the proprietary PFC solver.

Reusable periodic-coordinate transforms do not by themselves add periodic
rigid-particle contact assembly to `LSDEM` or `SphereDEM`. The native SPH
periodic API covers fluid/virtual-boundary neighbor images; a frontend still
owns solid wrapping, DEM seam contacts and the boundary policy.
