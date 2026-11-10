# FunDEMBeta: Quick Start for macOS

This guide gets FunDEMBeta building and running on the CPU on both Apple Silicon and Intel Macs. For the complete API and advanced workflows, see the [User Guide](docs/UserGuide.md).

Open a terminal in `workplace/`, the parent directory that will contain `FunDEMBeta/` and its sibling build directories. Keep the terminal in `workplace/` for every command in this guide.

> [!IMPORTANT]
> Modern macOS does not provide a supported NVIDIA CUDA environment. This guide explicitly disables CUDA. The core library, CPU tutorials and native C++ example remain available.

## 1. Install the build tools

Install the Apple Command Line Tools, which provide Apple Clang and Git:

```bash
xcode-select --install
```

If they are already installed, macOS will tell you that no installation is necessary.

Next, install [Homebrew](https://brew.sh/) if it is not already available, then install CMake and Ninja:

```bash
brew install cmake ninja
```

Verify the toolchain:

```bash
git --version
clang++ --version
cmake --version
ninja --version
```

FunDEMBeta requires CMake 3.24 or newer and a compiler with C++17 support.

## 2. Download the source

Keep the source and generated build files in separate directories:

```bash
git clone https://github.com/kaiqideng/FunDEMBeta.git
```

The `workplace/` directory will have this layout after configuration:

```text
workplace/
├── FunDEMBeta/   # Source code
└── build-macos/  # Generated build files
```

If you already have the source, open the terminal in its parent `workplace/` directory.

## 3. Configure and build

For the shortest first build, disable CUDA, OpenMP, validation cases, tests and
the optional C++ example. Build only the CPU tutorials:

```bash
cmake -S FunDEMBeta -B build-macos -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DFUNDEM_ENABLE_OPENMP=OFF \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DFUNDEM_BUILD_EXAMPLES=OFF \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF

cmake --build build-macos --target tutorial0 --parallel
```

After a successful build, the executable is located at:

```text
build-macos/tutorial/tutorial0
```

## 4. Run the first tutorial

Run the small sphere-rebound case:

```bash
build-macos/tutorial/tutorial0
```

The tutorial creates one sphere and a small analytic plane grid. It advances
3,000 steps over 0.3 seconds of physical time and writes 31 frames. Wall-clock
runtime depends on the machine and build. A successful return confirms that
the program initialized, advanced and wrote output; it is not an accuracy test.

The output is written beside the executable:

```text
build-macos/tutorial/tutorial0_files/
├── particle/     # Particle VTU files
└── energy.dat    # Solid DEM energy
```

Run settings are edited in `tutorial0.cpp`, then compiled again. For a shorter
diagnostic run, temporarily change `solve(3000)` to `solve(10)`, choose a fresh
folder in `setOutputDirectory()`, and rebuild. The core clears old VTU/DAT files
inside the selected output directory on the first `solve()`. Ten steps normally
produce only the initial frame; restore 3,000 steps and rebuild for the full run.

## 5. View the result in ParaView

Download ParaView from the [official website](https://www.paraview.org/download/), or install it with Homebrew:

```bash
brew install --cask paraview
```

For the complete run, open these numbered file series in ParaView:

```text
build-macos/tutorial/tutorial0_files/particle/sphere_*.vtu
build-macos/tutorial/tutorial0_files/particle/fixedLSParticle_*.vtu
```

Select the file series in the file dialog, click **Apply**, and then click **Reset Camera**. After the full simulation has finished, use **Play** to animate the result.

FunDEM writes binary VTU files by default. ParaView reads them directly; no conversion is required.

## 6. Optional: enable OpenMP

The minimal configuration disables OpenMP so that the first Apple Clang build needs no additional runtime setup. To enable CPU multithreading, install the Homebrew OpenMP runtime:

```bash
brew install libomp
```

Create a separate OpenMP build directory. `brew --prefix` selects the correct Homebrew path for Apple Silicon or Intel automatically:

```bash
LIBOMP_PREFIX="$(brew --prefix libomp)"

cmake -S FunDEMBeta -B build-macos-openmp -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DFUNDEM_ENABLE_OPENMP=ON \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF \
    -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp" \
    -DOpenMP_CXX_LIB_NAMES=omp \
    -DOpenMP_omp_LIBRARY="$LIBOMP_PREFIX/lib/libomp.dylib" \
    -DOpenMP_CXX_INCLUDE_DIR="$LIBOMP_PREFIX/include"

cmake --build build-macos-openmp --target tutorial0 --parallel
```

CMake enables OpenMP only when the compiler and runtime support the required
minimum/maximum reductions; finding the runtime alone is insufficient. If the
check fails, host loops remain serial. To select the thread count when enabled:

```bash
OMP_NUM_THREADS="$(sysctl -n hw.logicalcpu)" \
    build-macos-openmp/tutorial/tutorial0
```

Small simulations do not always become faster as the thread count increases. For performance measurements, compare several `OMP_NUM_THREADS` values.

## 7. Run the other tutorials

Build all tutorial targets:

```bash
cmake --build build-macos --parallel
```

Then run any of the remaining tutorials:

```bash
build-macos/tutorial/tutorial1
build-macos/tutorial/tutorial2
build-macos/tutorial/tutorial3
```

Keep `executionMode::CPU` in each source on macOS. Edit settings in `main()`
and rebuild; the executables take no options. GPU and Hybrid require CUDA.

See the [tutorial guide](tutorial/README.md) for each model, source setting and
output series.

## 8. Start a native C++ case

[`example/main.cpp`](example/main.cpp) contains a complete two-sphere collision
in one `main()`. Build it through the core:

```bash
cmake -S FunDEMBeta -B build-macos \
    -DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build-macos --target fundem_example --parallel
build-macos/example/fundem_example
```

Edit `main.cpp` to change the model and run settings, then rebuild. Keep CPU
execution on macOS. To start another case, copy the directory, give its target
a unique name and register it in the core CMake build. Follow the
[Programming Guide](docs/ProgrammingGuide.md) for the steps. Each run saves
native VTU and `energy.dat` in a fresh directory.

## 9. Troubleshooting

### `cmake: command not found`

Confirm that Homebrew is available, then run:

```bash
brew install cmake ninja
```

The default Apple Silicon Homebrew prefix is `/opt/homebrew`. If you have just installed Homebrew, run the shell setup command printed by its installer and open a new Terminal window in `workplace/`.

### CMake is older than 3.24

Check which CMake is active:

```bash
which cmake
cmake --version
```

Run `brew upgrade cmake`, then make sure the Homebrew executable directory appears before an older CMake installation in `PATH`.

### CMake still looks for CUDA

Make sure the configuration command contains:

```text
-DFUNDEM_ENABLE_CUDA=OFF
```

After changing the compiler, generator, or another major option, regenerate the build tree from scratch with CMake 3.24 or newer:

```bash
cmake --fresh -S FunDEMBeta -B build-macos -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DFUNDEM_ENABLE_OPENMP=OFF \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF
```

### OpenMP is not found

For the first build, leave `FUNDEM_ENABLE_OPENMP=OFF`. This only makes supported CPU loops run serially; it does not change the model.

To enable multithreading, confirm that Homebrew installed both the header and runtime library, then repeat the complete configuration command from Section 6:

```bash
test -f "$(brew --prefix libomp)/include/omp.h" && echo "omp.h found"
test -f "$(brew --prefix libomp)/lib/libomp.dylib" && echo "libomp found"
```

### The compiler and libraries have incompatible architectures

Check the architecture used by Terminal, Homebrew, and `libomp`:

```bash
uname -m
brew config
file "$(brew --prefix libomp)/lib/libomp.dylib"
```

On Apple Silicon, prefer a native ARM Homebrew installation and avoid mixing it with Intel libraries installed under Rosetta.

## Next steps

- Read the [tutorial guide](tutorial/README.md) for the four introductory cases.
- Read the [User Guide](docs/UserGuide.md) for modeling, contacts, solvers, output, and external-project integration.
- Return to the [project overview](README.md) for the source layout, GPU backend and reusable C++ example.
