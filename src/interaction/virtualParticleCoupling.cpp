#include "virtualParticleCoupling.h"

#include "execution/sphFunctions.h"
#include "execution/sphCouplingFunctions.h"
#include "execution/cpu/parallelPolicy.h"

#include <algorithm>
#include <stdexcept>

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "data/HostAoSDeviceSoA.h"
#endif

namespace fundem
{

void virtualParticleCoupling::reset() noexcept
{
    ownerOffsets_.clear();
    virtualParticleIndices_.clear();
    currentState_ = {};
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

    currentState_.velocitySums_.assign(virtualParticleCount, Vec3::zero());
    currentState_.accelerationSums_.assign(virtualParticleCount, Vec3::zero());
    currentState_.previousVelocities_.assign(virtualParticleCount, Vec3::zero());
    currentState_.forcesByOwner_.assign(LSParticleCount, Vec3::zero());
    currentState_.torquesByOwner_.assign(LSParticleCount, Vec3::zero());
    currentState_.forceIncrementsByOwner_.assign(LSParticleCount, Vec3::zero());
    currentState_.torqueIncrementsByOwner_.assign(LSParticleCount, Vec3::zero());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (virtualParticleCount >= cpu::parallelParticleThreshold)
#endif
    for (int virtualParticleIndex = 0; virtualParticleIndex < virtualParticleCount; ++virtualParticleIndex)
    {
        virtualParticle& value = virtualParticleHost[virtualParticleIndex];
        const LSParticle& owner = LSParticleHost[value.ownerLSParticleIndex()];
        Vec3 position;
        Vec3 normal;
        Vec3 velocity = Vec3::zero();
        Vec3 acceleration = Vec3::zero();
        execution::calculateVirtualParticlePositionAndNormal(position, normal, value.localPosition(), value.localNormal(), owner.position(), owner.orientation());
        execution::accumulateVirtualParticleVelocityAndAcceleration(velocity,
                                                                    acceleration,
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
        value.setPosition(position);
        value.setNormal(normal);
        value.setVelocity(velocity);
        value.setAcceleration(acceleration);
        currentState_.previousVelocities_[virtualParticleIndex] = velocity;
    }
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
            execution::accumulatePrescribedSPHKinematics(currentState_.velocitySums_[virtualParticleIndex], currentState_.accelerationSums_[virtualParticleIndex],
                                                        currentState_.previousVelocities_[virtualParticleIndex], velocity, timeStep);
            continue;
        }
        execution::accumulateVirtualParticleVelocityAndAcceleration(currentState_.velocitySums_[virtualParticleIndex],
                                                                    currentState_.accelerationSums_[virtualParticleIndex],
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
    currentState_.previousVelocities_.resize(samples.size());
    for (std::size_t index = 0; index < samples.size(); ++index)
    {
        const virtualParticle& sample = samples[index];
        const LSParticle& owner = LSParticles.host()[sample.ownerLSParticleIndex()];
        currentState_.previousVelocities_[index] = owner.velocity() + math::cross(owner.angularVelocity(), math::rotateUnit(owner.orientation(), sample.localPosition()));
    }
}

void virtualParticleCoupling::averageKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, int sampleCount)
{
    if (sampleCount <= 0)
    {
        throw std::invalid_argument("The virtual-particle kinematics sample count must be positive.");
    }
    auto& virtualParticleHost = virtualParticles.host();
    const auto& LSParticleHost = LSParticles.host();
    const int virtualParticleCount = static_cast<int>(virtualParticleHost.size());
    const math::Real inverseSampleCount = 1.0 / math::Real(sampleCount);
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
        value.setVelocity(inverseSampleCount * currentState_.velocitySums_[virtualParticleIndex]);
        value.setAcceleration(inverseSampleCount * currentState_.accelerationSums_[virtualParticleIndex]);
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
        currentState_.forceIncrementsByOwner_[ownerIndex] = force - currentState_.forcesByOwner_[ownerIndex];
        currentState_.torqueIncrementsByOwner_[ownerIndex] = torque - currentState_.torquesByOwner_[ownerIndex];
        currentState_.forcesByOwner_[ownerIndex] = force;
        currentState_.torquesByOwner_[ownerIndex] = torque;
    }
}

void virtualParticleCoupling::applyForceAndTorque(LSParticleContainer& LSParticles) const
{
    auto& LSParticleHost = LSParticles.host();
    const int LSParticleCount = static_cast<int>(LSParticleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (LSParticleCount >= cpu::parallelParticleThreshold)
#endif
    for (int ownerIndex = 0; ownerIndex < LSParticleCount; ++ownerIndex)
    {
        LSParticleHost[ownerIndex].addForce(currentState_.forcesByOwner_[ownerIndex]);
        LSParticleHost[ownerIndex].addTorque(currentState_.torquesByOwner_[ownerIndex]);
    }
}

void virtualParticleCoupling::clearKinematics() noexcept
{
    std::fill(currentState_.velocitySums_.begin(), currentState_.velocitySums_.end(), Vec3::zero());
    std::fill(currentState_.accelerationSums_.begin(), currentState_.accelerationSums_.end(), Vec3::zero());
}

void virtualParticleCoupling::applyImpulseCorrection(LSParticleContainer& LSParticles, math::Real correctionTime) const
{
    auto& owners = LSParticles.host();
    const int ownerCount = static_cast<int>(owners.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (ownerCount >= cpu::parallelParticleThreshold)
#endif
    for (int ownerIndex = 0; ownerIndex < ownerCount; ++ownerIndex)
    {
        LSParticle& owner = owners[ownerIndex];
        Vec3 velocityCorrection;
        Vec3 angularVelocityCorrection;
        execution::calculateSPHCouplingVelocityCorrection(velocityCorrection, angularVelocityCorrection,
                                                         currentState_.forceIncrementsByOwner_[ownerIndex], currentState_.torqueIncrementsByOwner_[ownerIndex],
                                                         owner.orientation(), owner.inverseMass(), owner.inverseInertiaTensor(), correctionTime);
        owner.setVelocity(owner.velocity() + velocityCorrection);
        owner.setAngularVelocity(owner.angularVelocity() + angularVelocityCorrection);
    }
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

void virtualParticleCoupling::prepareDeviceKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, int sampleCount, cudaStream_t stream)
{
    averageKinematics(virtualParticles, LSParticles, sampleCount);
    virtualParticles.copyHostFieldToDeviceAsync<pointMass::positionField>(stream);
    virtualParticles.copyHostFieldToDeviceAsync<virtualParticle::normalField>(stream);
    virtualParticles.copyHostFieldToDeviceAsync<pointMass::velocityField>(stream);
    virtualParticles.copyHostFieldToDeviceAsync<virtualParticle::accelerationField>(stream);
}

void virtualParticleCoupling::copyForceAndTorqueFromDevice(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream)
{
    const int virtualParticleCount = static_cast<int>(virtualParticles.deviceSize());
    const int ownerCount = ownerOffsets_.empty() ? 0 : static_cast<int>(ownerOffsets_.size()) - 1;
    if (virtualParticleCount == 0)
    {
        return;
    }

    std::vector<Vec3> virtualParticleForces(virtualParticleCount);
    host_device_detail::checkCuda(cudaMemcpyAsync(virtualParticleForces.data(), virtualParticles.device<pointMass::forceField>(), virtualParticleCount * sizeof(Vec3), cudaMemcpyDeviceToHost, stream),
                                  "cudaMemcpyAsync virtual-particle force to host");
    host_device_detail::synchronize(stream);

    const auto& virtualParticleHost = virtualParticles.host();
    const auto& LSParticleHost = LSParticles.host();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (ownerCount > 1 && virtualParticleCount >= cpu::parallelParticleThreshold)
#endif
    for (int ownerIndex = 0; ownerIndex < ownerCount; ++ownerIndex)
    {
        Vec3 force = Vec3::zero();
        Vec3 torque = Vec3::zero();
        for (int offset = ownerOffsets_[ownerIndex]; offset < ownerOffsets_[ownerIndex + 1]; ++offset)
        {
            const int virtualParticleIndex = virtualParticleIndices_[offset];
            const Vec3& virtualParticleForce = virtualParticleForces[virtualParticleIndex];
            force += virtualParticleForce;
            torque += math::cross(virtualParticleHost[virtualParticleIndex].position() - LSParticleHost[ownerIndex].position(), virtualParticleForce);
        }
        currentState_.forceIncrementsByOwner_[ownerIndex] = force - currentState_.forcesByOwner_[ownerIndex];
        currentState_.torqueIncrementsByOwner_[ownerIndex] = torque - currentState_.torquesByOwner_[ownerIndex];
        currentState_.forcesByOwner_[ownerIndex] = force;
        currentState_.torquesByOwner_[ownerIndex] = torque;
    }
}

#else

void virtualParticleCoupling::prepareDeviceKinematics(virtualParticleContainer&, const LSParticleContainer&, int, cudaStream_t)
{
    throw std::runtime_error("The virtual-particle device coupling is unavailable in this FunDEM build.");
}

void virtualParticleCoupling::copyForceAndTorqueFromDevice(const virtualParticleContainer&, const LSParticleContainer&, cudaStream_t)
{
    throw std::runtime_error("The virtual-particle device coupling is unavailable in this FunDEM build.");
}

#endif

} // namespace fundem
