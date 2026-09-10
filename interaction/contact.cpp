#include "contact.h"

#include "execution/contactFunctions.h"

namespace fundem
{

bool contact::calculateForce(Real timeStep) noexcept
{
    force_ = Vec3::zero();
    torque_ = Vec3::zero();
    normalForceMagnitude_ = 0.0;
    normalElasticEnergy_ = 0.0;
    slidingElasticEnergy_ = 0.0;
    rollingElasticEnergy_ = 0.0;
    torsionalElasticEnergy_ = 0.0;

    if (!isInContainer() || !particleCopiesUpdated_)
    {
        return false;
    }
    particleCopiesUpdated_ = false;

    const particle& master = masterParticle_;
    const particle& slave = slaveParticle_;
    const material& masterMaterial = master.getMaterial();
    const material& slaveMaterial = slave.getMaterial();
    const bool masterInfiniteMass = master.inverseMass() == 0.0;
    const bool slaveInfiniteMass = slave.inverseMass() == 0.0;
    const bool nodalContact = masterMaterial.isLevelSet();
    const bool sphereLevelSetContact = !nodalContact && slaveMaterial.isLevelSet();
    const Real stiffnessScale = nodalContact ? area_ : 1.0;
    const Real normalStiffness =
        sphereLevelSetContact
            ? masterMaterial.normalStiffness()
            : execution::effectiveStiffness(masterMaterial.normalStiffness() * stiffnessScale, slaveMaterial.normalStiffness() * stiffnessScale, masterInfiniteMass, slaveInfiniteMass);
    const Real slidingStiffness =
        sphereLevelSetContact
            ? masterMaterial.slidingStiffness()
            : execution::effectiveStiffness(masterMaterial.slidingStiffness() * stiffnessScale, slaveMaterial.slidingStiffness() * stiffnessScale, masterInfiniteMass, slaveInfiniteMass);
    const Real rollingStiffness = sphereLevelSetContact ? masterMaterial.rollingStiffness()
                                                        : execution::effectiveStiffness(masterMaterial.rollingStiffness(), slaveMaterial.rollingStiffness(), masterInfiniteMass, slaveInfiniteMass);
    const Real torsionalStiffness = sphereLevelSetContact
                                        ? masterMaterial.torsionalStiffness()
                                        : execution::effectiveStiffness(masterMaterial.torsionalStiffness(), slaveMaterial.torsionalStiffness(), masterInfiniteMass, slaveInfiniteMass);
    const Real slidingFrictionCoefficient = execution::harmonicMean(masterMaterial.slidingFrictionCoefficient(), slaveMaterial.slidingFrictionCoefficient());
    const Real rollingFrictionCoefficient = execution::harmonicMean(masterMaterial.rollingFrictionCoefficient(), slaveMaterial.rollingFrictionCoefficient());
    const Real torsionalFrictionCoefficient = execution::harmonicMean(masterMaterial.torsionalFrictionCoefficient(), slaveMaterial.torsionalFrictionCoefficient());
    const Real restitutionCoefficient = execution::harmonicMean(masterMaterial.restitutionCoefficient(), slaveMaterial.restitutionCoefficient());
    const Vec3 relativeVelocity =
        execution::relativeVelocityAtContactPoint(master.position(), master.velocity(), master.angularVelocity(), slave.position(), slave.velocity(), slave.angularVelocity(), point_);
    force_ = execution::calculateLinearContactForce(normalForceMagnitude_,
                                                    slidingSpringDeformation_,
                                                    relativeVelocity,
                                                    normal_,
                                                    overlap_,
                                                    timeStep,
                                                    normalStiffness,
                                                    slidingStiffness,
                                                    slidingFrictionCoefficient,
                                                    restitutionCoefficient,
                                                    effectiveMass_);
    normalElasticEnergy_ = overlap_ > 0.0 ? execution::springElasticEnergy(normalStiffness, overlap_) : 0.0;
    slidingElasticEnergy_ = execution::springElasticEnergy(slidingStiffness, slidingSpringDeformation_);

    if (nodalContact)
    {
        torque_ = Vec3::zero();
        rollingSpringDeformation_ = Vec3::zero();
        torsionalSpringDeformation_ = Vec3::zero();
        return true;
    }
    torque_ = execution::calculateSphereRotationalContactTorque(rollingSpringDeformation_,
                                                                torsionalSpringDeformation_,
                                                                master.angularVelocity() - slave.angularVelocity(),
                                                                normal_,
                                                                timeStep,
                                                                rollingStiffness,
                                                                torsionalStiffness,
                                                                rollingFrictionCoefficient,
                                                                torsionalFrictionCoefficient,
                                                                restitutionCoefficient,
                                                                effectiveMass_,
                                                                effectiveRadius_,
                                                                normalForceMagnitude_);
    rollingElasticEnergy_ = execution::springElasticEnergy(rollingStiffness, rollingSpringDeformation_);
    torsionalElasticEnergy_ = execution::springElasticEnergy(torsionalStiffness, torsionalSpringDeformation_);
    return true;
}

} // namespace fundem
