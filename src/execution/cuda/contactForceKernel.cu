#include "contactForceKernel.cuh"

#include "atomicKernel.cuh"
#include "data/HostAoSDeviceSoA.h"
#include "execution/contactFunctions.h"

namespace fundem::cuda::detail
{
namespace
{

constexpr int blockSize = 256;

__device__ materialPairCoefficients pairCoefficients(const contactMaterialPairDeviceView& pairs, int index) noexcept
{
    materialPairCoefficients result;
    result.normalStiffnessPerUnitArea_ = pairs.normalStiffnessPerUnitArea_[index];
    result.slidingStiffnessPerUnitArea_ = pairs.slidingStiffnessPerUnitArea_[index];
    result.normalStiffness_ = pairs.normalStiffness_[index];
    result.slidingStiffness_ = pairs.slidingStiffness_[index];
    result.rollingStiffness_ = pairs.rollingStiffness_[index];
    result.torsionalStiffness_ = pairs.torsionalStiffness_[index];
    result.slidingFrictionCoefficient_ = pairs.slidingFrictionCoefficient_[index];
    result.rollingFrictionCoefficient_ = pairs.rollingFrictionCoefficient_[index];
    result.torsionalFrictionCoefficient_ = pairs.torsionalFrictionCoefficient_[index];
    result.normalDampingCoefficient_ = pairs.normalDampingCoefficient_[index];
    result.slidingDampingCoefficient_ = pairs.slidingDampingCoefficient_[index];
    result.rollingDampingCoefficient_ = pairs.rollingDampingCoefficient_[index];
    result.torsionalDampingCoefficient_ = pairs.torsionalDampingCoefficient_[index];
    result.linearStiffnessEnabled_ = pairs.linearStiffnessEnabled_[index];
    result.dampingEnabled_ = pairs.dampingEnabled_[index];
    return result;
}

__global__ void contactForceKernel(contactDeviceView contacts,
                                   contactParticleDeviceView masterParticles,
                                   contactParticleDeviceView slaveParticles,
                                   contactMaterialDeviceView materials,
                                   contactMaterialPairDeviceView materialPairs,
                                   execution::contactType type,
                                   math::Real timeStep,
                                   int contactCount)
{
    const int contactIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (contactIndex >= contactCount)
    {
        return;
    }

    const int masterIndex = contacts.masterParticleIndex_[contactIndex];
    const int slaveIndex = contacts.slaveParticleIndex_[contactIndex];
    const int masterMaterialIndex = masterParticles.materialIndex_[masterIndex];
    const int slaveMaterialIndex = slaveParticles.materialIndex_[slaveIndex];
    const math::Real effectiveRadius = contacts.effectiveRadius_[contactIndex];
    const bool nodalContact = type == execution::contactType::levelSetLevelSet;
    const math::Real effectiveMass = contacts.effectiveMass_[contactIndex];
    const math::Real dampingWeight = nodalContact ? execution::nodalDampingWeight(effectiveMass, masterParticles.inverseMass_[masterIndex], slaveParticles.inverseMass_[slaveIndex]) : 1.0;
    const materialPairCoefficients pair = pairCoefficients(materialPairs, materialPairContainer::index(masterMaterialIndex, slaveMaterialIndex));
    const execution::contactResponse response = execution::resolveContactResponse(pair,
                                                                                  type,
                                                                                  materials.youngModulus_[masterMaterialIndex],
                                                                                  materials.poissonRatio_[masterMaterialIndex],
                                                                                  materials.restitutionCoefficient_[masterMaterialIndex],
                                                                                  materials.youngModulus_[slaveMaterialIndex],
                                                                                  materials.poissonRatio_[slaveMaterialIndex],
                                                                                  materials.restitutionCoefficient_[slaveMaterialIndex],
                                                                                  effectiveMass,
                                                                                  effectiveRadius,
                                                                                  contacts.overlap_[contactIndex],
                                                                                  contacts.area_[contactIndex],
                                                                                  dampingWeight);

    const math::Vec3 point = contacts.point_[contactIndex];
    const math::Vec3 normal = contacts.normal_[contactIndex];
    const math::Vec3 relativeVelocity = execution::relativeVelocityAtContactPoint(masterParticles.position_[masterIndex],
                                                                                  masterParticles.velocity_[masterIndex],
                                                                                  masterParticles.angularVelocity_[masterIndex],
                                                                                  slaveParticles.position_[slaveIndex],
                                                                                  slaveParticles.velocity_[slaveIndex],
                                                                                  slaveParticles.angularVelocity_[slaveIndex],
                                                                                  point);
    math::Real normalForceMagnitude = 0.0;
    math::Vec3 slidingSpringDeformation = contacts.slidingSpringDeformation_[contactIndex];
    const math::Vec3 force = execution::calculateLinearContactForce(normalForceMagnitude,
                                                                    slidingSpringDeformation,
                                                                    relativeVelocity,
                                                                    normal,
                                                                    contacts.overlap_[contactIndex],
                                                                    timeStep,
                                                                    response.normalStiffness_,
                                                                    response.slidingStiffness_,
                                                                    pair.slidingFrictionCoefficient_,
                                                                    response.normalDampingCoefficient_,
                                                                    response.slidingDampingCoefficient_);
    contacts.normalForceMagnitude_[contactIndex] = normalForceMagnitude;
    contacts.slidingSpringDeformation_[contactIndex] = slidingSpringDeformation;
    contacts.normalElasticEnergy_[contactIndex] = response.normalElasticEnergy_;
    contacts.slidingElasticEnergy_[contactIndex] = execution::springElasticEnergy(response.slidingStiffness_, slidingSpringDeformation);

    math::Vec3 rotationalTorque = math::Vec3::zero();
    if (nodalContact)
    {
        contacts.rollingSpringDeformation_[contactIndex] = math::Vec3::zero();
        contacts.torsionalSpringDeformation_[contactIndex] = math::Vec3::zero();
        contacts.rollingElasticEnergy_[contactIndex] = 0.0;
        contacts.torsionalElasticEnergy_[contactIndex] = 0.0;
    }
    else
    {
        math::Vec3 rollingSpringDeformation = contacts.rollingSpringDeformation_[contactIndex];
        math::Vec3 torsionalSpringDeformation = contacts.torsionalSpringDeformation_[contactIndex];
        rotationalTorque = execution::calculateSphereRotationalContactTorque(rollingSpringDeformation,
                                                                             torsionalSpringDeformation,
                                                                             masterParticles.angularVelocity_[masterIndex] - slaveParticles.angularVelocity_[slaveIndex],
                                                                             normal,
                                                                             timeStep,
                                                                             response.rollingStiffness_,
                                                                             response.torsionalStiffness_,
                                                                             pair.rollingFrictionCoefficient_,
                                                                             pair.torsionalFrictionCoefficient_,
                                                                             response.rollingDampingCoefficient_,
                                                                             response.torsionalDampingCoefficient_,
                                                                             effectiveRadius,
                                                                             normalForceMagnitude);
        contacts.rollingSpringDeformation_[contactIndex] = rollingSpringDeformation;
        contacts.torsionalSpringDeformation_[contactIndex] = torsionalSpringDeformation;
        contacts.rollingElasticEnergy_[contactIndex] = execution::springElasticEnergy(response.rollingStiffness_, rollingSpringDeformation);
        contacts.torsionalElasticEnergy_[contactIndex] = execution::springElasticEnergy(response.torsionalStiffness_, torsionalSpringDeformation);
    }

    contacts.force_[contactIndex] = force;
    contacts.torque_[contactIndex] = rotationalTorque;
    atomicAddVec3(masterParticles.force_, masterIndex, force);
    atomicAddVec3(masterParticles.torque_, masterIndex, math::cross(point - masterParticles.position_[masterIndex], force) + rotationalTorque);
    atomicAddVec3(slaveParticles.force_, slaveIndex, -force);
    atomicAddVec3(slaveParticles.torque_, slaveIndex, math::cross(point - slaveParticles.position_[slaveIndex], -force) - rotationalTorque);
}

} // namespace

void launchContactForce(contactDeviceView contacts,
                        contactParticleDeviceView masterParticles,
                        contactParticleDeviceView slaveParticles,
                        contactMaterialDeviceView materials,
                        contactMaterialPairDeviceView materialPairs,
                        execution::contactType type,
                        math::Real timeStep,
                        int contactCount,
                        cudaStream_t stream)
{
    if (contactCount <= 0)
    {
        return;
    }
    const int gridSize = (contactCount + blockSize - 1) / blockSize;
    contactForceKernel<<<gridSize, blockSize, 0, stream>>>(contacts, masterParticles, slaveParticles, materials, materialPairs, type, timeStep, contactCount);
    host_device_detail::checkCuda(cudaGetLastError(), "contactForceKernel launch");
}

} // namespace fundem::cuda::detail
