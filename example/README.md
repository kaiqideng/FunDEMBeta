# Minimal case template

`main.cpp` explains the model in eight numbered steps: parameters, solver,
material, contact pair, particles, environment, output and solve. Comments beside
the calls describe their arguments, return values and required order.

The only FunDEM include is `solver/solvers.h`; standard C++ headers for folders,
strings and console output remain explicit. Edit the parameters and select CPU,
GPU or Hybrid in the `SphereDEM` constructor; GPU/Hybrid require a CUDA-enabled core.
For another case, copy this directory, rename its target, and register the copy
with `add_subdirectory()` in the core's root `CMakeLists.txt`.
Follow the [Programming Guide](../docs/ProgrammingGuide.md) for the steps and
public API explanations.

Build through the core (run beside `FunDEMBeta/`):

```bash
cmake -S FunDEMBeta -B build -DFUNDEM_ENABLE_CUDA=OFF -DFUNDEM_BUILD_EXAMPLES=ON
cmake --build build --target fundem_example --parallel
OMP_NUM_THREADS=2 build/example/fundem_example
```

On Windows/MSVC, configure with `-DFUNDEM_ENABLE_OPENMP=OFF` for a serial CPU build
and omit `OMP_NUM_THREADS=2`; multi-config builds use `--config Release`
and `build/example/Release/fundem_example.exe`.

The default is a two-sphere collision in SI units: restitution 0.8, automatic
stiffness/damping, zero gravity/friction, dt=1e-6 s and 6000 steps. Output every
100 steps gives 61 frames. VTU and `energy.dat` go into a fresh
`example_result/run_N` beside the executable. This is a starting template;
numerical benchmarks live in [`validation/`](../validation/README.md).
