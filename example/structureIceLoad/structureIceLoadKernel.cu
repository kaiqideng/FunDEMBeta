#include "structureIceLoadKernel.cuh"

#include "data/HostAoSDeviceSoA.h"
#include "structureIceLoadFunctions.h"

namespace fundem::cuda
{
namespace
{

constexpr int blockSize = 256;

__global__ void sphereBuoyancyAndDragKernel(particle::device_type particles, math::Vec3 gravity, math::Vec3 waterVelocity, math::Real waterDensity, math::Real waterLevel, math::Real dragCoefficient)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particles.size_)
    {
        return;
    }

    execution::addSphereBuoyancyAndDrag(particles.force_[particleIndex],
                                        particles.torque_[particleIndex],
                                        particles.position_[particleIndex],
                                        particles.velocity_[particleIndex],
                                        particles.angularVelocity_[particleIndex],
                                        particles.radius_[particleIndex],
                                        gravity,
                                        waterVelocity,
                                        waterDensity,
                                        waterLevel,
                                        dragCoefficient);
}

__global__ void structureForceAccumulationKernel(math::Vec3* accumulatedForces, LSParticle::device_type structures)
{
    const int structureIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (structureIndex < structures.size_)
    {
        accumulatedForces[structureIndex] += structures.force_[structureIndex];
    }
}

} // namespace

void launchSphereBuoyancyAndDrag(particle::device_type particles,
                                 const math::Vec3& gravity,
                                 const math::Vec3& waterVelocity,
                                 math::Real waterDensity,
                                 math::Real waterLevel,
                                 math::Real dragCoefficient,
                                 cudaStream_t stream)
{
    if (particles.size_ <= 0)
    {
        return;
    }
    const int gridSize = (particles.size_ + blockSize - 1) / blockSize;
    sphereBuoyancyAndDragKernel<<<gridSize, blockSize, 0, stream>>>(particles, gravity, waterVelocity, waterDensity, waterLevel, dragCoefficient);
    host_device_detail::checkCuda(cudaGetLastError(), "sphereBuoyancyAndDragKernel launch");
}

void launchStructureForceAccumulation(math::Vec3* accumulatedForces, LSParticle::device_type structures, cudaStream_t stream)
{
    if (structures.size_ <= 0)
    {
        return;
    }
    const int gridSize = (structures.size_ + blockSize - 1) / blockSize;
    structureForceAccumulationKernel<<<gridSize, blockSize, 0, stream>>>(accumulatedForces, structures);
    host_device_detail::checkCuda(cudaGetLastError(), "structureForceAccumulationKernel launch");
}

void launchStructureForceClear(math::Vec3* accumulatedForces, int structureCount, cudaStream_t stream)
{
    if (structureCount <= 0)
    {
        return;
    }
    host_device_detail::checkCuda(cudaMemsetAsync(accumulatedForces, 0, static_cast<std::size_t>(structureCount) * sizeof(math::Vec3), stream), "clear structure force accumulator");
}

} // namespace fundem::cuda
