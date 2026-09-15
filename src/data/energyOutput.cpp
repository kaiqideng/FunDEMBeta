#include "energyOutput.h"

#include "datWriter.h"
#include "execution/motionIntegration.h"

namespace fundem
{
namespace
{

constexpr int parallelParticleThreshold = 1024;
constexpr int parallelInteractionThreshold = 256;

} // namespace

energyRecord::Real energyRecord::totalEnergy() const noexcept
{
    return kineticEnergy_ + gravitationalPotentialEnergy_ + contactNormalElasticEnergy_ + contactSlidingElasticEnergy_ + contactRollingElasticEnergy_ + contactTorsionalElasticEnergy_ +
           bondNormalElasticEnergy_ + bondShearElasticEnergy_ + bondBendingElasticEnergy_ + bondTorsionalElasticEnergy_;
}

template <class ParticleStorage> void addRigidBodyEnergy(energyRecord& result, const ParticleStorage& particles, const math::Vec3& gravity) noexcept
{
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    math::Real kineticEnergy = 0.0;
    math::Real gravitationalPotentialEnergy = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(+ : kineticEnergy, gravitationalPotentialEnergy) schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const auto& value = particleHost[particleIndex];
        kineticEnergy += execution::kineticEnergy(value.velocity(), value.angularVelocity(), value.orientation(), value.inverseMass(), value.inertiaTensor());
        gravitationalPotentialEnergy += execution::gravitationalPotentialEnergy(value.position(), value.inverseMass(), gravity);
    }
    result.kineticEnergy_ += kineticEnergy;
    result.gravitationalPotentialEnergy_ += gravitationalPotentialEnergy;
}

void addParticleEnergy(energyRecord& result, const particleContainer& particles, const math::Vec3& gravity) noexcept { addRigidBodyEnergy(result, particles, gravity); }

void addParticleEnergy(energyRecord& result, const LSParticleContainer& particles, const math::Vec3& gravity) noexcept { addRigidBodyEnergy(result, particles, gravity); }

void addParticleEnergy(energyRecord& result, const SPHParticleContainer& particles, const math::Vec3& gravity) noexcept
{
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    math::Real kineticEnergy = 0.0;
    math::Real gravitationalPotentialEnergy = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(+ : kineticEnergy, gravitationalPotentialEnergy) schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const SPHParticle& value = particleHost[particleIndex];
        kineticEnergy += value.inverseMass() > 0.0 ? 0.5 * math::normSquared(value.velocity()) / value.inverseMass() : 0.0;
        gravitationalPotentialEnergy += execution::gravitationalPotentialEnergy(value.position(), value.inverseMass(), gravity);
    }
    result.kineticEnergy_ += kineticEnergy;
    result.gravitationalPotentialEnergy_ += gravitationalPotentialEnergy;
}

void addContactEnergy(energyRecord& result, const contactContainer& contacts) noexcept
{
    const auto& contactHost = contacts.host();
    const int contactCount = static_cast<int>(contactHost.size());
    math::Real normalEnergy = 0.0;
    math::Real slidingEnergy = 0.0;
    math::Real rollingEnergy = 0.0;
    math::Real torsionalEnergy = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(+ : normalEnergy, slidingEnergy, rollingEnergy, torsionalEnergy) schedule(static) if (contactCount >= parallelInteractionThreshold)
#endif
    for (int contactIndex = 0; contactIndex < contactCount; ++contactIndex)
    {
        const contact& value = contactHost[contactIndex];
        normalEnergy += value.normalElasticEnergy();
        slidingEnergy += value.slidingElasticEnergy();
        rollingEnergy += value.rollingElasticEnergy();
        torsionalEnergy += value.torsionalElasticEnergy();
    }
    result.contactNormalElasticEnergy_ += normalEnergy;
    result.contactSlidingElasticEnergy_ += slidingEnergy;
    result.contactRollingElasticEnergy_ += rollingEnergy;
    result.contactTorsionalElasticEnergy_ += torsionalEnergy;
}

void addBondEnergy(energyRecord& result, const bondContainer& bonds) noexcept
{
    const auto& bondHost = bonds.host();
    const int bondCount = static_cast<int>(bondHost.size());
    math::Real normalEnergy = 0.0;
    math::Real shearEnergy = 0.0;
    math::Real bendingEnergy = 0.0;
    math::Real torsionalEnergy = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(+ : normalEnergy, shearEnergy, bendingEnergy, torsionalEnergy) schedule(static) if (bondCount >= parallelInteractionThreshold)
#endif
    for (int bondIndex = 0; bondIndex < bondCount; ++bondIndex)
    {
        const bond& value = bondHost[bondIndex];
        normalEnergy += value.normalElasticEnergy();
        shearEnergy += value.shearElasticEnergy();
        bendingEnergy += value.bendingElasticEnergy();
        torsionalEnergy += value.torsionalElasticEnergy();
    }
    result.bondNormalElasticEnergy_ += normalEnergy;
    result.bondShearElasticEnergy_ += shearEnergy;
    result.bondBendingElasticEnergy_ += bendingEnergy;
    result.bondTorsionalElasticEnergy_ += torsionalEnergy;
}

void appendEnergyDAT(const std::string& fileName, math::Real time, const energyRecord& energy, bool appendExisting)
{
    datWriter writer(fileName,
                     {"time",
                      "kineticEnergy",
                      "gravitationalPotentialEnergy",
                      "contactNormalElasticEnergy",
                      "contactSlidingElasticEnergy",
                      "contactRollingElasticEnergy",
                      "contactTorsionalElasticEnergy",
                      "bondNormalElasticEnergy",
                      "bondShearElasticEnergy",
                      "bondBendingElasticEnergy",
                      "bondTorsionalElasticEnergy",
                      "totalEnergy"},
                     appendExisting);
    writer.appendRow({time,
                      energy.kineticEnergy_,
                      energy.gravitationalPotentialEnergy_,
                      energy.contactNormalElasticEnergy_,
                      energy.contactSlidingElasticEnergy_,
                      energy.contactRollingElasticEnergy_,
                      energy.contactTorsionalElasticEnergy_,
                      energy.bondNormalElasticEnergy_,
                      energy.bondShearElasticEnergy_,
                      energy.bondBendingElasticEnergy_,
                      energy.bondTorsionalElasticEnergy_,
                      energy.totalEnergy()});
}

} // namespace fundem
