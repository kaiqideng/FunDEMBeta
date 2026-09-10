#include "bond.h"

#include "execution/bondFunctions.h"

namespace fundem
{

bool bond::calculateForce() noexcept
{
    if (!isInContainer() || !particleCopiesUpdated_)
    {
        force_ = Vec3::zero();
        masterTorque_ = Vec3::zero();
        slaveTorque_ = Vec3::zero();
        return false;
    }
    particleCopiesUpdated_ = false;

    const rigidBody& master = masterParticle_;
    const rigidBody& slave = slaveParticle_;
    const Real softCoefficient = 1.0 - damageFactor_;
    if (softCoefficient <= 0.0)
    {
        normal_ = Vec3::zero();
        force_ = Vec3::zero();
        masterTorque_ = Vec3::zero();
        slaveTorque_ = Vec3::zero();
        normalElasticEnergy_ = 0.0;
        shearElasticEnergy_ = 0.0;
        bendingElasticEnergy_ = 0.0;
        torsionalElasticEnergy_ = 0.0;
        return true;
    }

    const math::Quaternion& masterOrientation = master.orientation();
    const math::Quaternion& slaveOrientation = slave.orientation();
    const Vec3 masterEndpointNormal = math::rotateUnit(masterOrientation, masterEndpointLocalNormal_);
    const Vec3 masterEndpointTangent1 = math::rotateUnit(masterOrientation, masterEndpointLocalTangent1_);
    const Vec3 masterEndpointTangent2 = math::rotateUnit(masterOrientation, masterEndpointLocalTangent2_);
    const Vec3 slaveEndpointNormal = math::rotateUnit(slaveOrientation, slaveEndpointLocalNormal_);
    const Vec3 slaveEndpointTangent1 = math::rotateUnit(slaveOrientation, slaveEndpointLocalTangent1_);
    const Vec3 slaveEndpointTangent2 = math::rotateUnit(slaveOrientation, slaveEndpointLocalTangent2_);
    const Vec3 masterEndpointPosition = master.position() + math::rotateUnit(masterOrientation, masterEndpointLocalPosition_);
    const Vec3 slaveEndpointPosition = slave.position() + math::rotateUnit(slaveOrientation, slaveEndpointLocalPosition_);
    point_ = 0.5 * (masterEndpointPosition + slaveEndpointPosition);

    Vec3 masterEndpointTorque;
    Vec3 slaveEndpointTorque;
    Real axialDisplacement;
    execution::calculateBondForceAndTorque(force_,
                                           masterEndpointTorque,
                                           slaveEndpointTorque,
                                           normal_,
                                           axialDisplacement,
                                           normalElasticEnergy_,
                                           shearElasticEnergy_,
                                           bendingElasticEnergy_,
                                           torsionalElasticEnergy_,
                                           masterEndpointPosition,
                                           slaveEndpointPosition,
                                           masterEndpointNormal,
                                           masterEndpointTangent1,
                                           masterEndpointTangent2,
                                           slaveEndpointNormal,
                                           slaveEndpointTangent1,
                                           slaveEndpointTangent2,
                                           softCoefficient * coefficientB1_,
                                           softCoefficient * coefficientB2_,
                                           softCoefficient * coefficientB3_,
                                           softCoefficient * coefficientB4_,
                                           equivalentLength_);

    masterTorque_ = masterEndpointTorque + math::cross(masterEndpointPosition - master.position(), force_);
    slaveTorque_ = slaveEndpointTorque + math::cross(slaveEndpointPosition - slave.position(), -force_);

    if (execution::updateBKDamage(damageFactor_,
                                  maximumEnergyReleaseRatio_,
                                  axialDisplacement,
                                  normalElasticEnergy_,
                                  shearElasticEnergy_,
                                  bendingElasticEnergy_,
                                  torsionalElasticEnergy_,
                                  fractureArea_,
                                  modeICriticalEnergy_,
                                  modeIICriticalEnergy_,
                                  modeMixityExponent_,
                                  damageInitiationRatio_))
    {
        normal_ = Vec3::zero();
        normalElasticEnergy_ = 0.0;
        shearElasticEnergy_ = 0.0;
        bendingElasticEnergy_ = 0.0;
        torsionalElasticEnergy_ = 0.0;
    }
    return true;
}

} // namespace fundem
