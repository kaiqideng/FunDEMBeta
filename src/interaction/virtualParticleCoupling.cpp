#include "virtualParticleCoupling.h"

#include "execution/sphFunctions.h"
#include "execution/sphCouplingFunctions.h"
#include "execution/cpu/parallelPolicy.h"

#include <algorithm>
#include <stdexcept>

namespace fundem
{

void virtualParticleCoupling::reset() noexcept
{
    ownerOffsets_.clear();
    virtualParticleIndices_.clear();
    velocitySums_.clear();
    accelerationSums_.clear();
    previousVelocities_.clear();
    forcesByOwner_.clear();
    torquesByOwner_.clear();
}

void virtualParticleCoupling::initialize(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity)
{
    auto& virtualParticleHost = virtualParticles.host();
    const auto& LSParticleHost = LSParticles.host();
    const int virtualParticleCount = static_cast<int>(virtualParticleHost.size());
    const int LSParticleCount = static_cast<int>(LSParticleHost.size());

    ownerOffsets_.assign(LSParticleCount + 1, 0);
    for (const virtualParticle& value : virtualParticleHost)
    {
        ++ownerOffsets_[value.ownerLSParticleIndex() + 1];
    }
    for (int ownerIndex = 0; ownerIndex < LSParticleCount; ++ownerIndex)
    {
        ownerOffsets_[ownerIndex + 1] += ownerOffsets_[ownerIndex];
    }

    virtualParticleIndices_.resize(virtualParticleCount);
    std::vector<int> writeOffsets = ownerOffsets_;
    for (int virtualParticleIndex = 0; virtualParticleIndex < virtualParticleCount; ++virtualParticleIndex)
    {
        const int ownerIndex = virtualParticleHost[virtualParticleIndex].ownerLSParticleIndex();
        virtualParticleIndices_[writeOffsets[ownerIndex]++] = virtualParticleIndex;
    }

    velocitySums_.assign(virtualParticleCount, Vec3::zero());
    accelerationSums_.assign(virtualParticleCount, Vec3::zero());
    previousVelocities_.assign(virtualParticleCount, Vec3::zero());
    forcesByOwner_.assign(LSParticleCount, Vec3::zero());
    torquesByOwner_.assign(LSParticleCount, Vec3::zero());
    prepareKinematics(virtualParticles, LSParticles, gravity, 0);
    resetKinematicsReference(virtualParticles, LSParticles);
}

void virtualParticleCoupling::accumulateKinematics(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, math::Real timeStep)
{
    const auto& virtualParticleHost = virtualParticles.host();
    const auto& LSParticleHost = LSParticles.host();
    const int virtualParticleCount = static_cast<int>(virtualParticleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (virtualParticleCount >= cpu::parallelParticleThreshold)
#endif
    for (int virtualParticleIndex = 0; virtualParticleIndex < virtualParticleCount; ++virtualParticleIndex)
    {
        const virtualParticle& value = virtualParticleHost[virtualParticleIndex];
        const LSParticle& owner = LSParticleHost[value.ownerLSParticleIndex()];
        if (owner.inverseMass() <= 0.0 && timeStep > 0.0)
        {
            const Vec3 velocity = owner.velocity() + math::cross(owner.angularVelocity(), math::rotateUnit(owner.orientation(), value.localPosition()));
            execution::accumulatePrescribedSPHKinematics(velocitySums_[virtualParticleIndex], accelerationSums_[virtualParticleIndex], previousVelocities_[virtualParticleIndex], velocity, timeStep);
            continue;
        }
        execution::accumulateVirtualParticleVelocityAndAcceleration(velocitySums_[virtualParticleIndex],
                                                                    accelerationSums_[virtualParticleIndex],
                                                                    value.localPosition(),
                                                                    owner.orientation(),
                                                                    owner.velocity(),
                                                                    owner.angularVelocity(),
                                                                    owner.force(),
                                                                    owner.torque(),
                                                                    owner.inverseMass(),
                                                                    owner.inertiaTensor(),
                                                                    owner.inverseInertiaTensor(),
                                                                    gravity);
    }
}

void virtualParticleCoupling::resetKinematicsReference(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles)
{
    const auto& samples = virtualParticles.host();
    previousVelocities_.resize(samples.size());
    for (std::size_t index = 0; index < samples.size(); ++index)
    {
        const virtualParticle& sample = samples[index];
        const LSParticle& owner = LSParticles.host()[sample.ownerLSParticleIndex()];
        previousVelocities_[index] = owner.velocity() + math::cross(owner.angularVelocity(), math::rotateUnit(owner.orientation(), sample.localPosition()));
    }
}

void virtualParticleCoupling::prepareKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, int sampleCount)
{
    if (sampleCount < 0)
    {
        throw std::invalid_argument("The virtual-particle kinematics sample count must not be negative.");
    }
    auto& virtualParticleHost = virtualParticles.host();
    const auto& LSParticleHost = LSParticles.host();
    const int virtualParticleCount = static_cast<int>(virtualParticleHost.size());
    const math::Real inverseSampleCount = sampleCount > 0 ? 1.0 / math::Real(sampleCount) : 0.0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (virtualParticleCount >= cpu::parallelParticleThreshold)
#endif
    for (int virtualParticleIndex = 0; virtualParticleIndex < virtualParticleCount; ++virtualParticleIndex)
    {
        virtualParticle& value = virtualParticleHost[virtualParticleIndex];
        const LSParticle& owner = LSParticleHost[value.ownerLSParticleIndex()];
        Vec3 position;
        Vec3 normal;
        execution::calculateVirtualParticlePositionAndNormal(position, normal, value.localPosition(), value.localNormal(), owner.position(), owner.orientation());
        value.setPosition(position);
        value.setNormal(normal);
        Vec3 velocity = inverseSampleCount * velocitySums_[virtualParticleIndex];
        Vec3 acceleration = inverseSampleCount * accelerationSums_[virtualParticleIndex];
        if (sampleCount == 0)
        {
            execution::accumulateVirtualParticleVelocityAndAcceleration(velocity, acceleration, value.localPosition(), owner.orientation(), owner.velocity(), owner.angularVelocity(),
                                                                        owner.force(), owner.torque(), owner.inverseMass(), owner.inertiaTensor(), owner.inverseInertiaTensor(), gravity);
        }
        value.setVelocity(velocity);
        value.setAcceleration(acceleration);
    }
}

void virtualParticleCoupling::collectForceAndTorque(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles)
{
    const auto& virtualParticleHost = virtualParticles.host();
    const auto& LSParticleHost = LSParticles.host();
    const int ownerCount = ownerOffsets_.empty() ? 0 : static_cast<int>(ownerOffsets_.size()) - 1;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (ownerCount > 1 && static_cast<int>(virtualParticleHost.size()) >= cpu::parallelParticleThreshold)
#endif
    for (int ownerIndex = 0; ownerIndex < ownerCount; ++ownerIndex)
    {
        Vec3 force = Vec3::zero();
        Vec3 torque = Vec3::zero();
        for (int offset = ownerOffsets_[ownerIndex]; offset < ownerOffsets_[ownerIndex + 1]; ++offset)
        {
            const virtualParticle& value = virtualParticleHost[virtualParticleIndices_[offset]];
            force += value.force();
            torque += math::cross(value.position() - LSParticleHost[ownerIndex].position(), value.force());
        }
        forcesByOwner_[ownerIndex] = force;
        torquesByOwner_[ownerIndex] = torque;
    }
}

void virtualParticleCoupling::applyForceAndTorque(LSParticleContainer& LSParticles, math::Real scale) const
{
    auto& LSParticleHost = LSParticles.host();
    const int LSParticleCount = static_cast<int>(LSParticleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (LSParticleCount >= cpu::parallelParticleThreshold)
#endif
    for (int ownerIndex = 0; ownerIndex < LSParticleCount; ++ownerIndex)
    {
        LSParticleHost[ownerIndex].addForce(scale * forcesByOwner_[ownerIndex]);
        LSParticleHost[ownerIndex].addTorque(scale * torquesByOwner_[ownerIndex]);
    }
}

void virtualParticleCoupling::clearKinematics() noexcept
{
    std::fill(velocitySums_.begin(), velocitySums_.end(), Vec3::zero());
    std::fill(accelerationSums_.begin(), accelerationSums_.end(), Vec3::zero());
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

void virtualParticleCoupling::prepareDeviceKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, int sampleCount, cudaStream_t stream)
{
    prepareKinematics(virtualParticles, LSParticles, gravity, sampleCount);
    virtualParticles.copyHostFieldToDeviceAsync<pointMass::positionField>(stream);
    virtualParticles.copyHostFieldToDeviceAsync<virtualParticle::normalField>(stream);
    virtualParticles.copyHostFieldToDeviceAsync<pointMass::velocityField>(stream);
    virtualParticles.copyHostFieldToDeviceAsync<virtualParticle::accelerationField>(stream);
}

void virtualParticleCoupling::copyForceAndTorqueFromDevice(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream)
{
    if (virtualParticles.deviceSize() == 0)
    {
        return;
    }
    virtualParticles.copyDeviceFieldToHost<pointMass::forceField>(stream);
    collectForceAndTorque(virtualParticles, LSParticles);
}

#else

void virtualParticleCoupling::prepareDeviceKinematics(virtualParticleContainer&, const LSParticleContainer&, const Vec3&, int, cudaStream_t)
{
    throw std::runtime_error("The virtual-particle device coupling is unavailable in this FunDEM build.");
}

void virtualParticleCoupling::copyForceAndTorqueFromDevice(virtualParticleContainer&, const LSParticleContainer&, cudaStream_t)
{
    throw std::runtime_error("The virtual-particle device coupling is unavailable in this FunDEM build.");
}

#endif

} // namespace fundem
