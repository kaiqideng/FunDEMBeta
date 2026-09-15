#include "spatialGridKernel.cuh"

#include "data/HostAoSDeviceSoA.h"
#include "math/Indexing.h"

#include <thrust/device_ptr.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/transform.h>

#include <stdexcept>

namespace fundem::cuda::detail
{
namespace
{

constexpr int blockSize = 256;

struct equalInt {
    __host__ __device__ bool operator()(int first, int second) const noexcept { return first == second; }
};

struct addInt {
    __host__ __device__ int operator()(int first, int second) const noexcept { return first + second; }
};

__device__ int clampSpatialGridCoordinate(int coordinate, int size) noexcept { return coordinate < 0 ? 0 : (coordinate < size ? coordinate : size - 1); }

__device__ int calculateSpatialGridHash(const math::Vec3& position, const math::Vec3& minimumBoundary, const spatialGridLevel& level) noexcept
{
    const math::Vec3 local = math::hadamardProduct(position - minimumBoundary, level.inverseCellSize_);
    const int x = clampSpatialGridCoordinate(static_cast<int>(math::detail::floor(local.x)), level.size_.x);
    const int y = clampSpatialGridCoordinate(static_cast<int>(math::detail::floor(local.y)), level.size_.y);
    const int z = clampSpatialGridCoordinate(static_cast<int>(math::detail::floor(local.z)), level.size_.z);
    return math::linearIndex(x, y, z, level.size_.x, level.size_.y);
}

__global__ void calculateParticleHashKernel(spatialGridContainer::device_type spatialGrid, const math::Vec3* particlePositions, int levelIndex)
{
    const spatialGridLevel level = spatialGrid.levels_[levelIndex];
    const int gridParticleIndex = level.particleBegin_ + static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (gridParticleIndex >= level.particleEnd_)
    {
        return;
    }
    const int particleIndex = spatialGrid.particleIndex_[gridParticleIndex];
    spatialGrid.particleHash_[gridParticleIndex] = calculateSpatialGridHash(particlePositions[particleIndex], spatialGrid.minimumBoundary_, level);
}

__global__ void buildDenseCellRangesKernel(spatialGridContainer::device_type spatialGrid, int levelIndex)
{
    const spatialGridLevel& level = spatialGrid.levels_[levelIndex];
    const int localSortedParticleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int particleCount = level.particleEnd_ - level.particleBegin_;
    if (localSortedParticleIndex >= particleCount)
    {
        return;
    }

    const int sortedParticleIndex = level.particleBegin_ + localSortedParticleIndex;
    const int hash = spatialGrid.particleHash_[sortedParticleIndex];
    const int cellIndex = level.occupiedCellBegin_ + hash;
    if (localSortedParticleIndex == 0 || spatialGrid.particleHash_[sortedParticleIndex - 1] != hash)
    {
        spatialGrid.cellSortedParticleIndexBegin_[cellIndex] = sortedParticleIndex;
    }
    if (localSortedParticleIndex + 1 == particleCount || spatialGrid.particleHash_[sortedParticleIndex + 1] != hash)
    {
        spatialGrid.cellSortedParticleIndexEnd_[cellIndex] = sortedParticleIndex + 1;
    }
}

} // namespace

void launchBuildSpatialGrid(spatialGridContainer& spatialGrid, const math::Vec3* particlePositions, int particleCount, cudaStream_t stream)
{
    if (particleCount != spatialGrid.particleCount())
    {
        throw std::invalid_argument("Spatial-grid particle count does not match the particle container.");
    }

    spatialGridContainer::device_type spatialGridDevice = spatialGrid.mutableDevice();
    for (int levelIndex = 0; levelIndex < spatialGrid.levelCount(); ++levelIndex)
    {
        const spatialGridLevel& level = spatialGrid.level(levelIndex);
        const int levelParticleCount = level.particleEnd_ - level.particleBegin_;
        if (spatialGrid.usesDenseCells())
        {
            const int cellCount = spatialGrid.cellCount(levelIndex);
            const std::size_t cellBytes = static_cast<std::size_t>(cellCount) * sizeof(int);
            host_device_detail::checkCuda(cudaMemsetAsync(spatialGridDevice.cellSortedParticleIndexBegin_ + level.occupiedCellBegin_, 0xFF, cellBytes, stream),
                                          "cudaMemsetAsync dense spatial-grid cell begin");
            host_device_detail::checkCuda(cudaMemsetAsync(spatialGridDevice.cellSortedParticleIndexEnd_ + level.occupiedCellBegin_, 0xFF, cellBytes, stream),
                                          "cudaMemsetAsync dense spatial-grid cell end");
        }
        if (levelParticleCount <= 0)
        {
            if (!spatialGrid.usesDenseCells())
            {
                spatialGrid.setOccupiedCellCount(levelIndex, 0);
            }
            continue;
        }

        const int blockCount = (levelParticleCount + blockSize - 1) / blockSize;
        calculateParticleHashKernel<<<blockCount, blockSize, 0, stream>>>(spatialGridDevice, particlePositions, levelIndex);
        host_device_detail::checkCuda(cudaGetLastError(), "calculateParticleHashKernel launch");

        auto policy = thrust::cuda::par.on(stream);
        auto particleHash = thrust::device_pointer_cast(spatialGridDevice.particleHash_ + level.particleBegin_);
        auto particleIndex = thrust::device_pointer_cast(spatialGridDevice.particleIndex_ + level.particleBegin_);
        thrust::sort_by_key(policy, particleHash, particleHash + levelParticleCount, particleIndex);
        host_device_detail::checkCuda(cudaGetLastError(), "spatial-grid particle sort");

        if (spatialGrid.usesDenseCells())
        {
            buildDenseCellRangesKernel<<<blockCount, blockSize, 0, stream>>>(spatialGridDevice, levelIndex);
            host_device_detail::checkCuda(cudaGetLastError(), "buildDenseCellRangesKernel launch");
            continue;
        }

        auto cellHash = thrust::device_pointer_cast(spatialGridDevice.cellHash_ + level.occupiedCellBegin_);
        auto cellSortedParticleIndexEnd = thrust::device_pointer_cast(spatialGridDevice.cellSortedParticleIndexEnd_ + level.occupiedCellBegin_);
        const auto reducedEnd =
            thrust::reduce_by_key(policy, particleHash, particleHash + levelParticleCount, thrust::make_constant_iterator(1), cellHash, cellSortedParticleIndexEnd, equalInt{}, addInt{});
        const int occupiedCellCount = static_cast<int>(reducedEnd.first - cellHash);

        auto cellSortedParticleIndexBegin = thrust::device_pointer_cast(spatialGridDevice.cellSortedParticleIndexBegin_ + level.occupiedCellBegin_);
        thrust::exclusive_scan(policy, cellSortedParticleIndexEnd, cellSortedParticleIndexEnd + occupiedCellCount, cellSortedParticleIndexBegin, level.particleBegin_);
        thrust::transform(policy, cellSortedParticleIndexBegin, cellSortedParticleIndexBegin + occupiedCellCount, cellSortedParticleIndexEnd, cellSortedParticleIndexEnd, addInt{});
        host_device_detail::checkCuda(cudaGetLastError(), "spatial-grid occupied-cell ranges");
        spatialGrid.setOccupiedCellCount(levelIndex, occupiedCellCount);
    }
}

} // namespace fundem::cuda::detail
