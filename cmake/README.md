# `cmake/`

> **Purpose:** Store the installed package configuration used by `find_package(FunDEM CONFIG)`.

## Core Contents

| File | Responsibility |
| --- | --- |
| `FunDEMConfig.cmake.in` | Load `FunDEMTargets.cmake` and locate CUDA Toolkit and OpenMP dependencies required by the build |

The root `CMakeLists.txt` generates the version and exported-target files during configuration and installation. Core targets are defined under [`src/`](../src/README.md), which supplies the build-tree include root. This directory contains no solver code.

## Key Rules

- Add every new public library target to the root `FunDEMTargets` export set.
- Add new public dependencies here with `find_dependency()`.
- Use `BUILD_INTERFACE` and `INSTALL_INTERFACE` for build-tree and install-tree include paths.
- Install public headers through `FILE_SET public_headers` to preserve `FunDEM/<module>/...` include paths.

## Extension Guide

When adding an installable module, verify target export, public dependencies, and installed include paths, then test `find_package()` from a separate consumer project.

## Configuration and Export Flow

```text
Root CMake options
    -> detect optional CUDA and OpenMP support
    -> create module targets
    -> compose FunDEM::FunDEM
    -> install headers, libraries, and FunDEMTargets.cmake
    -> generate FunDEMConfig.cmake and version metadata
```

The installed configuration records whether the library was built with CUDA and OpenMP. A consuming project therefore resolves only the dependencies required by that installed package: `CUDAToolkit` is requested for CUDA builds, and `OpenMP::OpenMP_CXX` is requested when host parallelism was enabled successfully.

The primary consumer target is:

```cmake
find_package(FunDEM CONFIG REQUIRED)
target_link_libraries(my_case PRIVATE FunDEM::FunDEM)
```

Individual exported module targets are implementation-level building blocks. Applications should normally link the aggregate target so compile definitions such as `FUNDEM_HAS_CUDA`, include paths, the C++17 requirement, and transitive libraries remain consistent.

## Installation Checklist

1. Configure with `FUNDEM_INSTALL=ON` and an explicit `CMAKE_INSTALL_PREFIX`.
2. Build before invoking `cmake --install`.
3. Test the install from a separate source and build directory.
4. Point the consumer at the prefix with `CMAKE_PREFIX_PATH` rather than including FunDEM source directories manually.
5. Reinstall after adding a public header or changing exported dependencies.

If `find_package(FunDEM)` succeeds but a header is missing, confirm that the header belongs to a target `FILE_SET public_headers`. Private implementation headers, such as the host-only SPH interaction engine, intentionally are not installed.
