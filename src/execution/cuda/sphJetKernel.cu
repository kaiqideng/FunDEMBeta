#include "sphJetKernel.cuh"

#include "data/HostAoSDeviceSoA.h"
#include "execution/sphJetFunctions.h"

namespace fundem::cuda
{
namespace
{

constexpr int blockSize = 256;

__global__ void applySPHJetVelocityKernel(const math::Vec3* positions,
                                          math::Vec3* velocities,
                                          int* constrained,
                                          int particleBegin,
                                          int particleCount,
                                          bool active,
                                          math::Vec3 outletCenter,
                                          math::Vec3 unitDirection,
                                          math::Real radius,
                                          math::Real pipeLength,
                                          math::Vec3 velocity)
{
    const int localParticleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (localParticleIndex >= particleCount)
    {
        return;
    }

    const int particleIndex = particleBegin + localParticleIndex;
    const bool particleConstrained = active && execution::isInsideSPHJetPipe(positions[particleIndex], outletCenter, unitDirection, radius, pipeLength);
    constrained[particleIndex] = particleConstrained ? 1 : 0;
    if (particleConstrained)
    {
        velocities[particleIndex] = velocity;
    }
}

} // namespace

void launchApplySPHJetVelocity(SPHParticleContainer& particles,
                               int particleBegin,
                               int particleCount,
                               bool active,
                               const math::Vec3& outletCenter,
                               const math::Vec3& unitDirection,
                               math::Real radius,
                               math::Real pipeLength,
                               const math::Vec3& velocity,
                               cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }

    const int blockCount = (particleCount + blockSize - 1) / blockSize;
    applySPHJetVelocityKernel<<<blockCount, blockSize, 0, stream>>>(particles.device<pointMass::positionField>(),
                                                                    particles.device<pointMass::velocityField>(),
                                                                    particles.device<SPHParticle::constrainedField>(),
                                                                    particleBegin,
                                                                    particleCount,
                                                                    active,
                                                                    outletCenter,
                                                                    unitDirection,
                                                                    radius,
                                                                    pipeLength,
                                                                    velocity);
    host_device_detail::checkCuda(cudaGetLastError(), "applySPHJetVelocityKernel launch");
}

} // namespace fundem::cuda
