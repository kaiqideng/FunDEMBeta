#include "sphJetKernel.cuh"

#include "data/HostAoSDeviceSoA.h"
#include "execution/sphJetFunctions.h"

namespace fundem::cuda
{
namespace
{

constexpr int blockSize = 256;

__global__ void updateSPHJetOpeningsKernel(virtualParticle::device_type virtualParticles,
                                           const execution::SPHJetOpening::geometry* openings,
                                           int openingCount,
                                           math::Real nominalVolume,
                                           math::Real clearance,
                                           int* changed)
{
    const int index = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= virtualParticles.size_)
        return;
    math::Real volume = nominalVolume;
    for (int opening = 0; opening < openingCount; ++opening)
    {
        if (execution::isInsideSPHJetOpening(virtualParticles.position_[index], openings[opening], clearance))
        {
            volume = 0.0;
            break;
        }
    }
    if (virtualParticles.volume_[index] == volume)
        return;
    virtualParticles.volume_[index] = volume;
    virtualParticles.force_[index] = math::Vec3::zero();
    virtualParticles.priorForce_[index] = math::Vec3::zero();
    virtualParticles.couplingTorque_[index] = math::Vec3::zero();
    atomicExch(changed, 1);
}

__global__ void applySPHJetVelocityKernel(const math::Vec3* positions,
                                          math::Vec3* velocities,
                                          int* constrained,
                                          math::Real* densities,
                                          math::Real* pressures,
                                          math::Real* densityRates,
                                          int particleBegin,
                                          int particleCount,
                                          bool active,
                                          math::Vec3 outletCenter,
                                          math::Vec3 unitDirection,
                                          math::Real radius,
                                          math::Real pipeLength,
                                          math::Vec3 velocity,
                                          math::Real referenceDensity)
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
        densities[particleIndex] = referenceDensity;
        pressures[particleIndex] = 0.0;
        densityRates[particleIndex] = 0.0;
    }
}

} // namespace

bool launchUpdateSPHJetOpenings(virtualParticleContainer& virtualParticles,
                                 const execution::SPHJetOpeningContainer& openings,
                                 math::Real nominalVolume,
                                 math::Real clearance,
                                 cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count == 0)
        return false;
    host_device_detail::deviceArray<int> changedFlag;
    changedFlag.allocate(1);
    host_device_detail::checkCuda(cudaMemsetAsync(changedFlag.data(), 0, sizeof(int), stream), "cudaMemsetAsync SPH opening flag");
    updateSPHJetOpeningsKernel<<<(count + blockSize - 1) / blockSize, blockSize, 0, stream>>>(deviceFields(virtualParticles),
                                                                                         openings.device<execution::SPHJetOpening::valueField>(),
                                                                                         static_cast<int>(openings.deviceSize()),
                                                                                         nominalVolume,
                                                                                         clearance,
                                                                                         changedFlag.data());
    host_device_detail::checkCuda(cudaGetLastError(), "updateSPHJetOpeningsKernel launch");
    int changed = 0;
    host_device_detail::checkCuda(cudaMemcpyAsync(&changed, changedFlag.data(), sizeof(int), cudaMemcpyDeviceToHost, stream), "cudaMemcpyAsync SPH opening flag");
    host_device_detail::synchronize(stream);
    return changed != 0;
}

void launchApplySPHJetVelocity(SPHParticleContainer& particles,
                               int particleBegin,
                               int particleCount,
                               bool active,
                               const math::Vec3& outletCenter,
                               const math::Vec3& unitDirection,
                               math::Real radius,
                               math::Real pipeLength,
                               const math::Vec3& velocity,
                               math::Real referenceDensity,
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
                                                                    particles.device<SPHParticle::densityField>(),
                                                                    particles.device<SPHParticle::pressureField>(),
                                                                    particles.device<SPHParticle::densityRateField>(),
                                                                    particleBegin,
                                                                    particleCount,
                                                                    active,
                                                                    outletCenter,
                                                                    unitDirection,
                                                                    radius,
                                                                    pipeLength,
                                                                    velocity,
                                                                    referenceDensity);
    host_device_detail::checkCuda(cudaGetLastError(), "applySPHJetVelocityKernel launch");
}

} // namespace fundem::cuda
