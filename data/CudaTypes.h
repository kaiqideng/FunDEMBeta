/**
 * @file CudaTypes.h
 * @brief Provides CUDA-compatible stream and integer-vector types in every build.
 */
#pragma once

#ifndef FUNDEM_HAS_CUDA
#define FUNDEM_HAS_CUDA 0
#endif

#if FUNDEM_HAS_CUDA
#include <cuda_runtime.h>
#else
struct CUstream_st;
using cudaStream_t = CUstream_st*;

#if defined(__has_include)
#if __has_include(<vector_functions.h>)
#include <vector_functions.h>
#define FUNDEM_HAS_CUDA_VECTOR_TYPES 1
#endif
#endif

#if !defined(FUNDEM_HAS_CUDA_VECTOR_TYPES)
/** CUDA-compatible three-integer vector used when CUDA headers are unavailable. */
struct int3 {
    int x; ///< x component.
    int y; ///< y component.
    int z; ///< z component.
};

constexpr int3 make_int3(int x, int y, int z) noexcept { return {x, y, z}; }
#endif
#endif

#if defined(FUNDEM_HAS_CUDA_VECTOR_TYPES)
#undef FUNDEM_HAS_CUDA_VECTOR_TYPES
#endif
