/**
 * @file atomicKernel.cuh
 * @brief Defines CUDA atomic helpers for scalar and vector accumulation.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem::cuda::detail
{

/** Atomically adds a `Real`, including the pre-sm_60 double fallback. */
__device__ inline math::Real atomicAddReal(math::Real* address, math::Real increment) noexcept
{
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 600
    return atomicAdd(address, increment);
#else
    auto* addressAsInteger = reinterpret_cast<unsigned long long*>(address);
    unsigned long long oldValue = *addressAsInteger;
    unsigned long long assumedValue;
    do
    {
        assumedValue = oldValue;
        oldValue = atomicCAS(addressAsInteger, assumedValue, __double_as_longlong(increment + __longlong_as_double(assumedValue)));
    } while (oldValue != assumedValue);
    return __longlong_as_double(oldValue);
#endif
}

/** Atomically accumulates all components of @p increment into one vector. */
__device__ inline void atomicAddVec3(math::Vec3* values, int index, const math::Vec3& increment) noexcept
{
    atomicAddReal(&values[index].x, increment.x);
    atomicAddReal(&values[index].y, increment.y);
    atomicAddReal(&values[index].z, increment.z);
}

} // namespace fundem::cuda::detail
