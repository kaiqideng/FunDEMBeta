#include "bondForceKernel.cuh"

#include "atomicKernel.cuh"
#include "data/HostAoSDeviceSoA.h"
#include "execution/bondFunctions.h"

namespace fundem::cuda::detail
{
namespace
{

constexpr int blockSize = 256;

__device__ inline void clearFailedBondState(bondDeviceView bonds, int bondIndex) noexcept
{
    bonds.force_[bondIndex] = math::Vec3::zero();
    bonds.masterTorque_[bondIndex] = math::Vec3::zero();
    bonds.slaveTorque_[bondIndex] = math::Vec3::zero();
    bonds.normal_[bondIndex] = math::Vec3::zero();
    bonds.normalElasticEnergy_[bondIndex] = 0.0;
    bonds.shearElasticEnergy_[bondIndex] = 0.0;
    bonds.bendingElasticEnergy_[bondIndex] = 0.0;
    bonds.torsionalElasticEnergy_[bondIndex] = 0.0;
}

__device__ inline void clearNewlyFailedBondEnergy(bondDeviceView bonds, int bondIndex) noexcept
{
    bonds.normal_[bondIndex] = math::Vec3::zero();
    bonds.normalElasticEnergy_[bondIndex] = 0.0;
    bonds.shearElasticEnergy_[bondIndex] = 0.0;
    bonds.bendingElasticEnergy_[bondIndex] = 0.0;
    bonds.torsionalElasticEnergy_[bondIndex] = 0.0;
}

__global__ void bondForceKernel(bondDeviceView bonds, bondParticleDeviceView masterParticles, bondParticleDeviceView slaveParticles, int bondCount)
{
    const int bondIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (bondIndex >= bondCount)
    {
        return;
    }

    const math::Real softCoefficient = 1.0 - bonds.damageFactor_[bondIndex];
    if (softCoefficient <= 0.0)
    {
        clearFailedBondState(bonds, bondIndex);
        return;
    }

    const int masterIndex = bonds.masterParticleIndex_[bondIndex];
    const int slaveIndex = bonds.slaveParticleIndex_[bondIndex];
    const math::Quaternion masterOrientation = masterParticles.orientation_[masterIndex];
    const math::Quaternion slaveOrientation = slaveParticles.orientation_[slaveIndex];
    const math::Vec3 masterEndpointNormal = math::rotateUnit(masterOrientation, bonds.masterEndpointLocalNormal_[bondIndex]);
    const math::Vec3 masterEndpointTangent1 = math::rotateUnit(masterOrientation, bonds.masterEndpointLocalTangent1_[bondIndex]);
    const math::Vec3 masterEndpointTangent2 = math::rotateUnit(masterOrientation, bonds.masterEndpointLocalTangent2_[bondIndex]);
    const math::Vec3 slaveEndpointNormal = math::rotateUnit(slaveOrientation, bonds.slaveEndpointLocalNormal_[bondIndex]);
    const math::Vec3 slaveEndpointTangent1 = math::rotateUnit(slaveOrientation, bonds.slaveEndpointLocalTangent1_[bondIndex]);
    const math::Vec3 slaveEndpointTangent2 = math::rotateUnit(slaveOrientation, bonds.slaveEndpointLocalTangent2_[bondIndex]);
    const math::Vec3 masterEndpointPosition = masterParticles.position_[masterIndex] + math::rotateUnit(masterOrientation, bonds.masterEndpointLocalPosition_[bondIndex]);
    const math::Vec3 slaveEndpointPosition = slaveParticles.position_[slaveIndex] + math::rotateUnit(slaveOrientation, bonds.slaveEndpointLocalPosition_[bondIndex]);
    bonds.point_[bondIndex] = 0.5 * (masterEndpointPosition + slaveEndpointPosition);

    math::Vec3 force;
    math::Vec3 masterEndpointTorque;
    math::Vec3 slaveEndpointTorque;
    math::Vec3 normal;
    math::Real axialDisplacement;
    math::Real normalElasticEnergy;
    math::Real shearElasticEnergy;
    math::Real bendingElasticEnergy;
    math::Real torsionalElasticEnergy;
    execution::calculateBondForceAndTorque(force,
                                           masterEndpointTorque,
                                           slaveEndpointTorque,
                                           normal,
                                           axialDisplacement,
                                           normalElasticEnergy,
                                           shearElasticEnergy,
                                           bendingElasticEnergy,
                                           torsionalElasticEnergy,
                                           masterEndpointPosition,
                                           slaveEndpointPosition,
                                           masterEndpointNormal,
                                           masterEndpointTangent1,
                                           masterEndpointTangent2,
                                           slaveEndpointNormal,
                                           slaveEndpointTangent1,
                                           slaveEndpointTangent2,
                                           softCoefficient * bonds.coefficientB1_[bondIndex],
                                           softCoefficient * bonds.coefficientB2_[bondIndex],
                                           softCoefficient * bonds.coefficientB3_[bondIndex],
                                           softCoefficient * bonds.coefficientB4_[bondIndex],
                                           bonds.equivalentLength_[bondIndex]);
    bonds.normal_[bondIndex] = normal;
    bonds.normalElasticEnergy_[bondIndex] = normalElasticEnergy;
    bonds.shearElasticEnergy_[bondIndex] = shearElasticEnergy;
    bonds.bendingElasticEnergy_[bondIndex] = bendingElasticEnergy;
    bonds.torsionalElasticEnergy_[bondIndex] = torsionalElasticEnergy;

    const math::Vec3 masterTorque = masterEndpointTorque + math::cross(masterEndpointPosition - masterParticles.position_[masterIndex], force);
    const math::Vec3 slaveTorque = slaveEndpointTorque + math::cross(slaveEndpointPosition - slaveParticles.position_[slaveIndex], -force);
    bonds.force_[bondIndex] = force;
    bonds.masterTorque_[bondIndex] = masterTorque;
    bonds.slaveTorque_[bondIndex] = slaveTorque;
    atomicAddVec3(masterParticles.force_, masterIndex, force);
    atomicAddVec3(masterParticles.torque_, masterIndex, masterTorque);
    atomicAddVec3(slaveParticles.force_, slaveIndex, -force);
    atomicAddVec3(slaveParticles.torque_, slaveIndex, slaveTorque);

    const bool failed = execution::updateBKDamage(bonds.damageFactor_[bondIndex],
                                                  bonds.maximumEnergyReleaseRatio_[bondIndex],
                                                  axialDisplacement,
                                                  normalElasticEnergy,
                                                  shearElasticEnergy,
                                                  bendingElasticEnergy,
                                                  torsionalElasticEnergy,
                                                  bonds.crossSectionArea_[bondIndex],
                                                  bonds.modeICriticalEnergy_[bondIndex],
                                                  bonds.modeIICriticalEnergy_[bondIndex],
                                                  bonds.modeMixityExponent_[bondIndex],
                                                  bonds.damageInitiationRatio_[bondIndex]);
    if (failed)
    {
        clearNewlyFailedBondEnergy(bonds, bondIndex);
    }
}

} // namespace

void launchBondForce(bondDeviceView bonds, bondParticleDeviceView masterParticles, bondParticleDeviceView slaveParticles, int bondCount, cudaStream_t stream)
{
    if (bondCount <= 0)
    {
        return;
    }
    const int gridSize = (bondCount + blockSize - 1) / blockSize;
    bondForceKernel<<<gridSize, blockSize, 0, stream>>>(bonds, masterParticles, slaveParticles, bondCount);
    host_device_detail::checkCuda(cudaGetLastError(), "bondForceKernel launch");
}

} // namespace fundem::cuda::detail
