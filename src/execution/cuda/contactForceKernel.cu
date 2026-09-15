#include "contactForceKernel.cuh"

#include "atomicKernel.cuh"
#include "data/HostAoSDeviceSoA.h"
#include "execution/contactFunctions.h"

namespace fundem::cuda::detail
{
namespace
{

constexpr int blockSize = 256;

__global__ void contactForceKernel(contactDeviceView contacts,
                                   contactParticleDeviceView masterParticles,
                                   contactParticleDeviceView slaveParticles,
                                   contactMaterialDeviceView materials,
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
    const bool masterInfiniteMass = masterParticles.inverseMass_[masterIndex] == 0.0;
    const bool slaveInfiniteMass = slaveParticles.inverseMass_[slaveIndex] == 0.0;
    const math::Real effectiveRadius = contacts.effectiveRadius_[contactIndex];
    const bool nodalContact = materials.type_[masterMaterialIndex] == materialType::levelSet;
    const bool sphereLevelSetContact = !nodalContact && materials.type_[slaveMaterialIndex] == materialType::levelSet;
    const math::Real stiffnessScale = nodalContact ? contacts.area_[contactIndex] : 1.0;

    const math::Real normalStiffness = sphereLevelSetContact ? materials.normalStiffness_[masterMaterialIndex]
                                                             : execution::effectiveStiffness(materials.normalStiffness_[masterMaterialIndex] * stiffnessScale,
                                                                                             materials.normalStiffness_[slaveMaterialIndex] * stiffnessScale,
                                                                                             masterInfiniteMass,
                                                                                             slaveInfiniteMass);
    const math::Real slidingStiffness = sphereLevelSetContact ? materials.slidingStiffness_[masterMaterialIndex]
                                                              : execution::effectiveStiffness(materials.slidingStiffness_[masterMaterialIndex] * stiffnessScale,
                                                                                              materials.slidingStiffness_[slaveMaterialIndex] * stiffnessScale,
                                                                                              masterInfiniteMass,
                                                                                              slaveInfiniteMass);
    const math::Real rollingStiffness =
        sphereLevelSetContact ? materials.rollingStiffness_[masterMaterialIndex]
                              : execution::effectiveStiffness(materials.rollingStiffness_[masterMaterialIndex], materials.rollingStiffness_[slaveMaterialIndex], masterInfiniteMass, slaveInfiniteMass);
    const math::Real torsionalStiffness =
        sphereLevelSetContact
            ? materials.torsionalStiffness_[masterMaterialIndex]
            : execution::effectiveStiffness(materials.torsionalStiffness_[masterMaterialIndex], materials.torsionalStiffness_[slaveMaterialIndex], masterInfiniteMass, slaveInfiniteMass);
    const math::Real slidingFrictionCoefficient = execution::harmonicMean(materials.slidingFrictionCoefficient_[masterMaterialIndex], materials.slidingFrictionCoefficient_[slaveMaterialIndex]);
    const math::Real rollingFrictionCoefficient = execution::harmonicMean(materials.rollingFrictionCoefficient_[masterMaterialIndex], materials.rollingFrictionCoefficient_[slaveMaterialIndex]);
    const math::Real torsionalFrictionCoefficient = execution::harmonicMean(materials.torsionalFrictionCoefficient_[masterMaterialIndex], materials.torsionalFrictionCoefficient_[slaveMaterialIndex]);
    const math::Real restitutionCoefficient = execution::harmonicMean(materials.restitutionCoefficient_[masterMaterialIndex], materials.restitutionCoefficient_[slaveMaterialIndex]);

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
                                                                    normalStiffness,
                                                                    slidingStiffness,
                                                                    slidingFrictionCoefficient,
                                                                    restitutionCoefficient,
                                                                    contacts.effectiveMass_[contactIndex]);
    contacts.normalForceMagnitude_[contactIndex] = normalForceMagnitude;
    contacts.slidingSpringDeformation_[contactIndex] = slidingSpringDeformation;
    contacts.normalElasticEnergy_[contactIndex] = contacts.overlap_[contactIndex] > 0.0 ? execution::springElasticEnergy(normalStiffness, contacts.overlap_[contactIndex]) : 0.0;
    contacts.slidingElasticEnergy_[contactIndex] = execution::springElasticEnergy(slidingStiffness, slidingSpringDeformation);

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
                                                                             rollingStiffness,
                                                                             torsionalStiffness,
                                                                             rollingFrictionCoefficient,
                                                                             torsionalFrictionCoefficient,
                                                                             restitutionCoefficient,
                                                                             contacts.effectiveMass_[contactIndex],
                                                                             effectiveRadius,
                                                                             normalForceMagnitude);
        contacts.rollingSpringDeformation_[contactIndex] = rollingSpringDeformation;
        contacts.torsionalSpringDeformation_[contactIndex] = torsionalSpringDeformation;
        contacts.rollingElasticEnergy_[contactIndex] = execution::springElasticEnergy(rollingStiffness, rollingSpringDeformation);
        contacts.torsionalElasticEnergy_[contactIndex] = execution::springElasticEnergy(torsionalStiffness, torsionalSpringDeformation);
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
                        math::Real timeStep,
                        int contactCount,
                        cudaStream_t stream)
{
    if (contactCount <= 0)
    {
        return;
    }
    const int gridSize = (contactCount + blockSize - 1) / blockSize;
    contactForceKernel<<<gridSize, blockSize, 0, stream>>>(contacts, masterParticles, slaveParticles, materials, timeStep, contactCount);
    host_device_detail::checkCuda(cudaGetLastError(), "contactForceKernel launch");
}

} // namespace fundem::cuda::detail
