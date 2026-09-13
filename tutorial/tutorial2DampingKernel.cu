/**
 * @file tutorial2DampingKernel.cu
 * @brief Applies the bonded cloth tutorial's damping on the solver CUDA stream.
 */
#include "tutorial2DampingKernel.cuh"

#include "data/HostAoSDeviceSoA.h"

namespace fundem::tutorial
{
namespace
{

constexpr int blockSize = 256;

/** Accumulates independent translation and rotation damping for one finite sphere. */
__global__ void clothDampingKernel(particle::device_type particles, math::Real dampingRate, math::Real rotationalDragCoefficient)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particles.size_ || !(particles.inverseMass_[particleIndex] > 0.0))
    {
        return;
    }

    const math::Real velocityCoefficient = dampingRate / particles.inverseMass_[particleIndex];
    particles.force_[particleIndex] -= velocityCoefficient * particles.velocity_[particleIndex];
    particles.torque_[particleIndex] -= rotationalDragCoefficient * particles.angularVelocity_[particleIndex];
}

} // namespace

void launchClothDamping(particle::device_type particles, math::Real dampingRate, math::Real rotationalDragCoefficient, cudaStream_t stream)
{
    if (particles.size_ <= 0)
    {
        return;
    }
    const int gridSize = (particles.size_ + blockSize - 1) / blockSize;
    clothDampingKernel<<<gridSize, blockSize, 0, stream>>>(particles, dampingRate, rotationalDragCoefficient);
    host_device_detail::checkCuda(cudaGetLastError(), "clothDampingKernel launch");
}

} // namespace fundem::tutorial
