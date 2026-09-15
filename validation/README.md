# Validation cases

This directory separates numerical verification and physical validation from tutorials and unit tests. A tutorial demonstrates an API workflow; a validation case defines a reference quantity, an error measure, and an explicit acceptance criterion.

## Build and run

Validation executables are enabled by default. Open a terminal in `workplace/`, the directory containing `FunDEMBeta/`, and keep it there for all commands:

```bash
cmake -S FunDEMBeta -B build
cmake --build build --target freeFall normalCollision frictionLaw bondModes LSCentralCompression backendConsistency --parallel
```

Run the complete suite:

```bash
cmake --build build --target runValidations
```

Or run one case:

```bash
build/validation/freeFall
build/validation/normalCollision
build/validation/frictionLaw
build/validation/bondModes
build/validation/LSCentralCompression
build/validation/backendConsistency
```

Every executable prints `PASS`, `FAIL`, or an explicit backend `SKIP`. By default, numerical histories are written beside the executable in `<case>_files/`. Use `--output DIR` to select another directory.

Set `-DFUNDEM_BUILD_VALIDATIONS=OFF` only when a minimal library build is required.

## Implemented cases

| Executable | Reference and primary checks | Default output |
| --- | --- | --- |
| `freeFall` | Constant-acceleration translation, constant-axis torque, angular velocity, and orientation | `freeFall.dat` |
| `normalCollision` | Unilateral linear-dashpot rebound, analytical contact duration, and pair momentum | `normalCollision.dat` |
| `frictionLaw` | Tangential history, Coulomb-disk return, and history removal when contact opens | `frictionLaw.dat` |
| `bondModes` | Single-bond axial, shear, pure-bending, and pure-torsion responses | `bondModes.dat` |
| `LSCentralCompression` | Two LS spheres, nodal-force integration, surface refinement, and the small-overlap continuum reference | `LSCentralCompression.dat` |
| `backendConsistency` | Identical sphere–LS trajectory in CPU, GPU, and Sphere-GPU/LS-CPU Hybrid modes | `backendConsistency.dat` |

`normalCollision` reports both the configured restitution parameter and the analytical rebound after the noncohesive normal force is truncated at zero. These values differ slightly because damping stops before geometric separation; the acceptance check uses the analytical response of the force law that is actually implemented.

`backendConsistency` always records the CPU result. GPU comparisons are skipped only when CUDA was not built or no CUDA device is present.

## Interpreting the LS compression case

Two equal spheres of radius `R` are compressed by `delta`. For the one-sided nodal penalty integral used here, the leading small-overlap force is

```text
P = (pi / 2) R k_n delta^2,
```

where `k_n` is stiffness per unit area. The DAT file contains the numerical force, this continuum reference, relative error, and active nodal-contact count for four icosphere surface resolutions. It is a discretization check, not a Hertz-material calibration.

## Required next layer

The following cases need traceable laboratory data and must not be treated as analytical unit checks:

1. confined compression for bulk stiffness;
2. direct shear at several normal stresses for peak strength, steady strength, and dilation;
3. triaxial compression at several confining pressures;
4. angle of repose or rotating drum as an independent bulk-flow validation;
5. hopper discharge or granular-column collapse as a predictive flow case;
6. bonded Brazilian splitting and mixed-mode fracture after the elastic bond modes pass.

Each physical case should keep calibration data separate from validation data, preserve the experimental particle-size and shape distributions, report time-step and resolution studies, and compare measured curves rather than a single final number.
