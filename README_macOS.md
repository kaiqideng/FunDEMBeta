# FunDEMBeta: Quick Start for macOS

This guide gets FunDEMBeta building and running on the CPU on both Apple Silicon and Intel Macs. For the complete API and advanced workflows, see the [User Guide](docs/UserGuide.md).

> [!IMPORTANT]
> Modern macOS does not provide a supported NVIDIA CUDA environment. This guide explicitly disables CUDA. CPU tutorials and CPU-compatible examples remain fully available.

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
mkdir -p ~/FunDEM-workspace
cd ~/FunDEM-workspace
git clone https://github.com/kaiqideng/FunDEMBeta.git
cd FunDEMBeta
```

The workspace will have this layout after configuration:

```text
FunDEM-workspace/
├── FunDEMBeta/   # Source code
└── build-macos/  # Generated build files
```

If you already have the source, start in its parent directory and adjust the path passed to `-S` in the commands below.

## 3. Configure and build

For the shortest first build, disable CUDA, OpenMP, large examples, validation cases, and tests. Build only the CPU tutorials:

```bash
cd ~/FunDEM-workspace/FunDEMBeta

cmake -S . -B ../build-macos -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DFUNDEM_ENABLE_OPENMP=OFF \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DFUNDEM_BUILD_EXAMPLES=OFF \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF

cmake --build ../build-macos --target tutorial0 --parallel
```

After a successful build, the executable is located at:

```text
../build-macos/tutorial/tutorial0
```

## 4. Run a smoke test

Run ten time steps to confirm that the program can initialize, advance the simulation, and write output:

```bash
../build-macos/tutorial/tutorial0 \
    --steps 10 \
    --output tutorial0_smoke
```

On its first run, the program must generate the Gömböc level-set grid. This can take a minute or longer on some Macs, and the terminal may remain quiet during this stage.

When the solver prints its status and returns to the prompt without an error, the smoke test has passed. The output is written beside the executable:

```text
../build-macos/tutorial/tutorial0_smoke/
├── particle/     # Particle VTU files
└── energy.dat    # System energy
```

A ten-step smoke test normally contains only the initial output frame. To generate the complete Gömböc self-righting simulation, run:

```bash
../build-macos/tutorial/tutorial0
```

The complete case simulates 120 seconds of physical time and takes much longer than the smoke test.

## 5. View the result in ParaView

Download ParaView from the [official website](https://www.paraview.org/download/), or install it with Homebrew:

```bash
brew install --cask paraview
```

For the complete run, open these numbered file series in ParaView:

```text
../build-macos/tutorial/tutorial0_files/particle/LSParticle_*.vtu
../build-macos/tutorial/tutorial0_files/particle/fixedLSParticle_*.vtu
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
cd ~/FunDEM-workspace/FunDEMBeta
LIBOMP_PREFIX="$(brew --prefix libomp)"

cmake -S . -B ../build-macos-openmp -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DFUNDEM_ENABLE_CUDA=OFF \
    -DFUNDEM_ENABLE_OPENMP=ON \
    -DFUNDEM_BUILD_TUTORIALS=ON \
    -DFUNDEM_BUILD_EXAMPLES=OFF \
    -DFUNDEM_BUILD_VALIDATIONS=OFF \
    -DBUILD_TESTING=OFF \
    -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp" \
    -DOpenMP_CXX_LIB_NAMES=omp \
    -DOpenMP_omp_LIBRARY="$LIBOMP_PREFIX/lib/libomp.dylib" \
    -DOpenMP_CXX_INCLUDE_DIR="$LIBOMP_PREFIX/include"

cmake --build ../build-macos-openmp --target tutorial0 --parallel
```

CMake should report that OpenMP was found. To select the thread count for one run:

```bash
OMP_NUM_THREADS="$(sysctl -n hw.logicalcpu)" \
    ../build-macos-openmp/tutorial/tutorial0 --steps 1000
```

Small simulations do not always become faster as the thread count increases. For performance measurements, compare several `OMP_NUM_THREADS` values.

## 7. Run the other tutorials

Build all tutorial targets:

```bash
cmake --build ../build-macos --parallel
```

Then run any of the remaining tutorials:

```bash
../build-macos/tutorial/tutorial1
../build-macos/tutorial/tutorial2 --steps 1000
../build-macos/tutorial/tutorial3 --steps 1000
```

Keep the default CPU execution mode on macOS. Do not pass `--mode gpu` or `--mode hybrid`; those modes require CUDA.

See the [tutorial guide](tutorial/README.md) for the physical model, expected runtime, and command-line options of each case.

## 8. Troubleshooting

### `cmake: command not found`

Confirm that Homebrew is available, then run:

```bash
brew install cmake ninja
```

The default Apple Silicon Homebrew prefix is `/opt/homebrew`. If you have just installed Homebrew, run the shell setup command printed by its installer and open a new Terminal window.

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
    -DFUNDEM_BUILD_EXAMPLES=OFF \
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
- Return to the [project overview](README.md) for the source layout, GPU backend, and larger examples.
