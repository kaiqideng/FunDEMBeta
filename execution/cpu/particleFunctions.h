/**
 * @file particleFunctions.h
 * @brief Provides OpenMP-enabled CPU loops for integration and force assembly.
 */
#pragma once

#include "execution/cpu/parallelPolicy.h"
#include "interaction/interactionContainer.h"

#include <numeric>
#include <vector>

namespace fundem::cpu
{

/** Compressed interaction indices incident to each particle. */
struct interactionAdjacency {
    std::vector<int> offsets_;            ///< Particle range offsets with one sentinel.
    std::vector<int> interactionIndices_; ///< Flattened incident contact or bond indices.
};

/** Builds a compressed particle-to-interaction adjacency from active interaction flags. */
template <class IndexGetter> interactionAdjacency buildInteractionAdjacency(int particleCount, const std::vector<unsigned char>& active, IndexGetter indexGetter)
{
    interactionAdjacency result;
    result.offsets_.assign(particleCount + 1, 0);
    const int interactionCount = static_cast<int>(active.size());
    for (int interactionIndex = 0; interactionIndex < interactionCount; ++interactionIndex)
    {
        if (active[interactionIndex] != 0)
        {
            ++result.offsets_[indexGetter(interactionIndex) + 1];
        }
    }
    std::partial_sum(result.offsets_.begin(), result.offsets_.end(), result.offsets_.begin());
    result.interactionIndices_.resize(result.offsets_.back());
    std::vector<int> writeOffsets(result.offsets_.begin(), result.offsets_.end() - 1);
    for (int interactionIndex = 0; interactionIndex < interactionCount; ++interactionIndex)
    {
        if (active[interactionIndex] != 0)
        {
            result.interactionIndices_[writeOffsets[indexGetter(interactionIndex)]++] = interactionIndex;
        }
    }
    return result;
}

/** Clears force and torque on every host particle, using OpenMP above the configured threshold. */
template <class ParticleStorage> void clearForceAndTorque(ParticleStorage& particles)
{
    auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        auto& value = particleHost[particleIndex];
        value.setForce(math::Vec3::zero());
        value.setTorque(math::Vec3::zero());
    }
}

/** Advances host linear and angular velocities in parallel when profitable. */
template <class ParticleStorage> void integrateVelocity(ParticleStorage& particles, const math::Vec3& gravity, math::Real timeStep)
{
    auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        particleHost[particleIndex].updateVelocityAndAngularVelocity(gravity, timeStep);
    }
}

/** Advances host positions and orientations in parallel when profitable. */
template <class ParticleStorage> void integratePosition(ParticleStorage& particles, math::Real timeStep)
{
    auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        particleHost[particleIndex].updatePositionAndOrientation(timeStep);
    }
}

/** Evaluates contacts and race-freely assembles their forces and torques onto host particles. */
template <class MasterStorage, class SlaveStorage> void addContactForceAndTorque(contactContainer& contacts, MasterStorage& masterParticles, SlaveStorage& slaveParticles, math::Real historyTimeStep)
{
    auto& contactHost = contacts.host();
    auto& masterParticleHost = masterParticles.host();
    auto& slaveParticleHost = slaveParticles.host();
    const int contactCount = static_cast<int>(contactHost.size());
    if (contactCount < parallelInteractionThreshold)
    {
        for (contact& value : contactHost)
        {
            if (!value.calculateForce(historyTimeStep))
            {
                continue;
            }
            auto& master = masterParticleHost[value.masterParticleIndex()];
            auto& slave = slaveParticleHost[value.slaveParticleIndex()];
            master.addForce(value.force());
            master.addTorque(math::cross(value.point() - master.position(), value.force()) + value.torque());
            slave.addForce(-value.force());
            slave.addTorque(math::cross(value.point() - slave.position(), -value.force()) - value.torque());
        }
        return;
    }

    std::vector<unsigned char> calculated(contactCount, 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int contactIndex = 0; contactIndex < contactCount; ++contactIndex)
    {
        calculated[contactIndex] = contactHost[contactIndex].calculateForce(historyTimeStep) ? 1 : 0;
    }

    const interactionAdjacency masterAdjacency =
        buildInteractionAdjacency(static_cast<int>(masterParticleHost.size()), calculated, [&](int contactIndex) { return contactHost[contactIndex].masterParticleIndex(); });
    const interactionAdjacency slaveAdjacency =
        buildInteractionAdjacency(static_cast<int>(slaveParticleHost.size()), calculated, [&](int contactIndex) { return contactHost[contactIndex].slaveParticleIndex(); });
#if defined(_OPENMP)
#pragma omp parallel
    {
#pragma omp for schedule(static)
#endif
        for (int particleIndex = 0; particleIndex < static_cast<int>(masterParticleHost.size()); ++particleIndex)
        {
            auto& master = masterParticleHost[particleIndex];
            math::Vec3 force = math::Vec3::zero();
            math::Vec3 torque = math::Vec3::zero();
            for (int adjacencyIndex = masterAdjacency.offsets_[particleIndex]; adjacencyIndex < masterAdjacency.offsets_[particleIndex + 1]; ++adjacencyIndex)
            {
                const contact& value = contactHost[masterAdjacency.interactionIndices_[adjacencyIndex]];
                force += value.force();
                torque += math::cross(value.point() - master.position(), value.force()) + value.torque();
            }
            master.addForce(force);
            master.addTorque(torque);
        }
#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
        for (int particleIndex = 0; particleIndex < static_cast<int>(slaveParticleHost.size()); ++particleIndex)
        {
            auto& slave = slaveParticleHost[particleIndex];
            math::Vec3 force = math::Vec3::zero();
            math::Vec3 torque = math::Vec3::zero();
            for (int adjacencyIndex = slaveAdjacency.offsets_[particleIndex]; adjacencyIndex < slaveAdjacency.offsets_[particleIndex + 1]; ++adjacencyIndex)
            {
                const contact& value = contactHost[slaveAdjacency.interactionIndices_[adjacencyIndex]];
                force -= value.force();
                torque += math::cross(value.point() - slave.position(), -value.force()) - value.torque();
            }
            slave.addForce(force);
            slave.addTorque(torque);
        }
#if defined(_OPENMP)
    }
#endif
}

/** Refreshes and evaluates bonds, then race-freely assembles their loads onto host particles. */
template <class MasterStorage, class SlaveStorage> void addBondForceAndTorque(bondContainer& bonds, MasterStorage& masterParticles, SlaveStorage& slaveParticles)
{
    auto& bondHost = bonds.host();
    auto& masterParticleHost = masterParticles.host();
    auto& slaveParticleHost = slaveParticles.host();
    const int bondCount = static_cast<int>(bondHost.size());
    if (bondCount < parallelInteractionThreshold)
    {
        for (bond& value : bondHost)
        {
            if (!value.updateMasterSlaveParticle(masterParticles, slaveParticles) || !value.calculateForce())
            {
                continue;
            }
            auto& master = masterParticleHost[value.masterParticleIndex()];
            auto& slave = slaveParticleHost[value.slaveParticleIndex()];
            master.addForce(value.force());
            master.addTorque(value.masterTorque());
            slave.addForce(-value.force());
            slave.addTorque(value.slaveTorque());
        }
        return;
    }

    std::vector<unsigned char> calculated(bondCount, 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (int bondIndex = 0; bondIndex < bondCount; ++bondIndex)
    {
        bond& value = bondHost[bondIndex];
        calculated[bondIndex] = value.updateMasterSlaveParticle(masterParticles, slaveParticles) && value.calculateForce() ? 1 : 0;
    }

    const interactionAdjacency masterAdjacency =
        buildInteractionAdjacency(static_cast<int>(masterParticleHost.size()), calculated, [&](int bondIndex) { return bondHost[bondIndex].masterParticleIndex(); });
    const interactionAdjacency slaveAdjacency =
        buildInteractionAdjacency(static_cast<int>(slaveParticleHost.size()), calculated, [&](int bondIndex) { return bondHost[bondIndex].slaveParticleIndex(); });
#if defined(_OPENMP)
#pragma omp parallel
    {
#pragma omp for schedule(static)
#endif
        for (int particleIndex = 0; particleIndex < static_cast<int>(masterParticleHost.size()); ++particleIndex)
        {
            auto& master = masterParticleHost[particleIndex];
            math::Vec3 force = math::Vec3::zero();
            math::Vec3 torque = math::Vec3::zero();
            for (int adjacencyIndex = masterAdjacency.offsets_[particleIndex]; adjacencyIndex < masterAdjacency.offsets_[particleIndex + 1]; ++adjacencyIndex)
            {
                const bond& value = bondHost[masterAdjacency.interactionIndices_[adjacencyIndex]];
                force += value.force();
                torque += value.masterTorque();
            }
            master.addForce(force);
            master.addTorque(torque);
        }
#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
        for (int particleIndex = 0; particleIndex < static_cast<int>(slaveParticleHost.size()); ++particleIndex)
        {
            auto& slave = slaveParticleHost[particleIndex];
            math::Vec3 force = math::Vec3::zero();
            math::Vec3 torque = math::Vec3::zero();
            for (int adjacencyIndex = slaveAdjacency.offsets_[particleIndex]; adjacencyIndex < slaveAdjacency.offsets_[particleIndex + 1]; ++adjacencyIndex)
            {
                const bond& value = bondHost[slaveAdjacency.interactionIndices_[adjacencyIndex]];
                force -= value.force();
                torque += value.slaveTorque();
            }
            slave.addForce(force);
            slave.addTorque(torque);
        }
#if defined(_OPENMP)
    }
#endif
}

/** Same-container convenience overload for bond load assembly. */
template <class ParticleStorage> void addBondForceAndTorque(bondContainer& bonds, ParticleStorage& particles) { addBondForceAndTorque(bonds, particles, particles); }

} // namespace fundem::cpu
