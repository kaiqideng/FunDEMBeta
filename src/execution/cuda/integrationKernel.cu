#include "integrationKernel.cuh"

#include "data/HostAoSDeviceSoA.h"
#include "execution/motionIntegration.h"

namespace fundem::cuda::detail
{
namespace
{

constexpr int blockSize = 256;

__global__ void velocityAndAngularVelocityIntegrationKernel(integrationDeviceView particles, math::Vec3 gravity, math::Real timeStep, int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount || particles.inverseMass_[particleIndex] == 0.0)
    {
        return;
    }

    execution::integrateVelocity(particles.velocity_[particleIndex], particles.force_[particleIndex], particles.inverseMass_[particleIndex], gravity, timeStep);
    execution::integrateAngularVelocity(particles.angularVelocity_[particleIndex],
                                        particles.torque_[particleIndex],
                                        particles.orientation_[particleIndex],
                                        particles.inertiaTensor_[particleIndex],
                                        particles.inverseInertiaTensor_[particleIndex],
                                        timeStep);
}

__global__ void positionAndOrientationIntegrationKernel(integrationDeviceView particles, math::Real timeStep, int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }

    execution::integratePosition(particles.position_[particleIndex], particles.velocity_[particleIndex], timeStep);
    execution::integrateQuaternion(particles.orientation_[particleIndex], particles.angularVelocity_[particleIndex], timeStep);
}

} // namespace

void launchVelocityAndAngularVelocityIntegration(integrationDeviceView particles, const math::Vec3& gravity, math::Real timeStep, int particleCount, cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }
    const int gridSize = (particleCount + blockSize - 1) / blockSize;
    velocityAndAngularVelocityIntegrationKernel<<<gridSize, blockSize, 0, stream>>>(particles, gravity, timeStep, particleCount);
    host_device_detail::checkCuda(cudaGetLastError(), "velocityAndAngularVelocityIntegrationKernel launch");
}

void launchPositionAndOrientationIntegration(integrationDeviceView particles, math::Real timeStep, int particleCount, cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }
    const int gridSize = (particleCount + blockSize - 1) / blockSize;
    positionAndOrientationIntegrationKernel<<<gridSize, blockSize, 0, stream>>>(particles, timeStep, particleCount);
    host_device_detail::checkCuda(cudaGetLastError(), "positionAndOrientationIntegrationKernel launch");
}

} // namespace fundem::cuda::detail
