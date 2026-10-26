#include "contact.h"

#include "execution/contactFunctions.h"

namespace fundem
{

bool contact::calculateForce(Real timeStep, const materialPairContainer& materialPairs, execution::contactType type) noexcept
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
    const int masterMaterialIndex = master.materialIndex();
    const int slaveMaterialIndex = slave.materialIndex();
    if (masterMaterialIndex < 0 || slaveMaterialIndex < 0 || masterMaterialIndex >= materialPairs.materialCount() || slaveMaterialIndex >= materialPairs.materialCount())
        return false;
    const materialPair& pair = materialPairs.host()[materialPairContainer::index(masterMaterialIndex, slaveMaterialIndex)];
    const bool nodalContact = type == execution::contactType::levelSetLevelSet;
    if (nodalContact && !pair.hasLevelSetStiffness())
        return false;
    const materialPairCoefficients coefficients = pair.coefficients();
    const execution::contactResponse response = execution::resolveContactResponse(coefficients,
                                                                                  type,
                                                                                  masterMaterial.youngModulus(),
                                                                                  masterMaterial.poissonRatio(),
                                                                                  masterMaterial.restitutionCoefficient(),
                                                                                  slaveMaterial.youngModulus(),
                                                                                  slaveMaterial.poissonRatio(),
                                                                                  slaveMaterial.restitutionCoefficient(),
                                                                                  effectiveMass_,
                                                                                  effectiveRadius_,
                                                                                  overlap_,
                                                                                  area_,
                                                                                  nodalContact ? execution::nodalDampingWeight(effectiveMass_, master.inverseMass(), slave.inverseMass()) : 1.0);
    const Vec3 relativeVelocity =
        execution::relativeVelocityAtContactPoint(master.position(), master.velocity(), master.angularVelocity(), slave.position(), slave.velocity(), slave.angularVelocity(), point_);
    force_ = execution::calculateLinearContactForce(normalForceMagnitude_,
                                                    slidingSpringDeformation_,
                                                    relativeVelocity,
                                                    normal_,
                                                    overlap_,
                                                    timeStep,
                                                    response.normalStiffness_,
                                                    response.slidingStiffness_,
                                                    coefficients.slidingFrictionCoefficient_,
                                                    response.normalDampingCoefficient_,
                                                    response.slidingDampingCoefficient_);
    normalElasticEnergy_ = response.normalElasticEnergy_;
    slidingElasticEnergy_ = execution::springElasticEnergy(response.slidingStiffness_, slidingSpringDeformation_);

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
                                                                response.rollingStiffness_,
                                                                response.torsionalStiffness_,
                                                                coefficients.rollingFrictionCoefficient_,
                                                                coefficients.torsionalFrictionCoefficient_,
                                                                response.rollingDampingCoefficient_,
                                                                response.torsionalDampingCoefficient_,
                                                                effectiveRadius_,
                                                                normalForceMagnitude_);
    rollingElasticEnergy_ = execution::springElasticEnergy(response.rollingStiffness_, rollingSpringDeformation_);
    torsionalElasticEnergy_ = execution::springElasticEnergy(response.torsionalStiffness_, torsionalSpringDeformation_);
    return true;
}

} // namespace fundem
