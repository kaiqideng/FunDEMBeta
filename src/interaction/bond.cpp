#include "bond.h"

#include "execution/bondFunctions.h"

namespace fundem
{

bond::Vec3 bond::particleTorque(bool master) const noexcept
{
    if (!isInContainer() || damageFactor_ >= 1.0 || normal_ == Vec3::zero())
        return Vec3::zero();

    const void* container = master ? masterParticleContainer_ : slaveParticleContainer_;
    const int index = master ? masterParticleIndex_ : slaveParticleIndex_;
    if (container == nullptr || index < 0 || type_ == bondType::undefined)
        return Vec3::zero();

    const rigidBody* particle = nullptr;
    if (type_ == bondType::sphereSphere || (master && type_ == bondType::sphereLSParticle))
    {
        const auto& particles = *static_cast<const particleContainer*>(container);
        if (static_cast<std::size_t>(index) < particles.hostSize())
            particle = &particles.host()[index];
    }
    else
    {
        const auto& particles = *static_cast<const LSParticleContainer*>(container);
        if (static_cast<std::size_t>(index) < particles.hostSize())
            particle = &particles.host()[index];
    }
    if (particle == nullptr)
        return Vec3::zero();

    const Vec3 endpointNormal = math::rotateUnit(particle->orientation(), master ? masterEndpointLocalNormal_ : slaveEndpointLocalNormal_);
    const Vec3 endpointOffset = math::rotateUnit(particle->orientation(), master ? masterEndpointLocalPosition_ : slaveEndpointLocalPosition_);
    const Vec3 moment = 0.5 * coefficientB2_ * (1.0 - damageFactor_) * math::cross(normal_, endpointNormal) + torque_ + math::cross(endpointOffset, force_);
    return master ? moment : -moment;
}

bool bond::calculateForce(bool updateDamage) noexcept
{
    if (!isInContainer() || !particleCopiesUpdated_)
    {
        force_ = Vec3::zero();
        torque_ = Vec3::zero();
        normal_ = Vec3::zero();
        return false;
    }
    particleCopiesUpdated_ = false;

    const rigidBody& master = masterParticle_;
    const rigidBody& slave = slaveParticle_;
    if (damageFactor_ >= 1.0)
    {
        normal_ = Vec3::zero();
        force_ = Vec3::zero();
        torque_ = Vec3::zero();
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
                                           torque_,
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
                                           coefficientB1_,
                                           coefficientB2_,
                                           coefficientB3_,
                                           coefficientB4_,
                                           bendingStiffness_,
                                           equivalentLength_);

    const bool failed = updateDamage && fractureMode_ != fundem::fractureMode::unbreakable && execution::updateBondDamage(damageFactor_,
                                                                                                          damageStarted_,
                                                                                                          onsetEffectiveStress_,
                                                                                                          onsetEffectiveDisplacement_,
                                                                                                          maximumEffectiveDisplacement_,
                                                                                                          failureEffectiveDisplacement_,
                                                                                                          effectiveFractureEnergy_,
                                                                                                          force_,
                                                                                                          torque_,
                                                                                                          normal_,
                                                                                                          normalElasticEnergy_,
                                                                                                          shearElasticEnergy_,
                                                                                                          bendingElasticEnergy_,
                                                                                                          torsionalElasticEnergy_,
                                                                                                          crossSectionArea_,
                                                                                                          modeICriticalEnergy_,
                                                                                                          modeIICriticalEnergy_,
                                                                                                          modeMixityExponent_,
                                                                                                          tensileStrength_,
                                                                                                          cohesion_,
                                                                                                          internalFrictionCoefficient_,
                                                                                                          fractureMode_ == fundem::fractureMode::damage);

    execution::applyBondDamageToResponse(force_,
                                         torque_,
                                         masterEndpointTorque,
                                         slaveEndpointTorque,
                                         normalElasticEnergy_,
                                         shearElasticEnergy_,
                                         bendingElasticEnergy_,
                                         torsionalElasticEnergy_,
                                         normal_,
                                         axialDisplacement,
                                         normalStiffness_,
                                         damageFactor_);
    if (failed)
        normal_ = Vec3::zero();
    return true;
}

} // namespace fundem
