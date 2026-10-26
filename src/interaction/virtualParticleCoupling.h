/**
 * @file virtualParticleCoupling.h
 * @brief Collects virtual-boundary kinematics and SPH reaction loads by LS owner.
 */
#pragma once

#include "data/CudaTypes.h"
#include "particle/LSParticle.h"
#include "particle/virtualParticle.h"

#include <vector>

namespace fundem
{

/** Owns host-side virtual-particle kinematics and resultant SPH wall loads. */
class virtualParticleCoupling
{
public:
    using Vec3 = math::Vec3;

    /** Releases owner mappings and all accumulator state. */
    void reset() noexcept;
    /** Generates virtual particles and stable owner-to-sample ranges. */
    void initialize(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity);
    /** Samples complete left-endpoint LS kinematics, before force assembly clears the loads. */
    void accumulateKinematics(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, math::Real timeStep = 0.0);
    /** Starts the prescribed-motion difference at the complete initialized endpoint. */
    void resetKinematicsReference(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles);
    /** Refreshes the current boundary pose, using completed-interval averages or instantaneous kinematics when @p sampleCount is zero. */
    void prepareKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, int sampleCount);
    /** Reduces virtual-particle reaction forces and torques by LS owner. */
    void collectForceAndTorque(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles);
    /** Averages host kinematics and uploads only the fields required by device SPH. */
    void prepareDeviceKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, int sampleCount, cudaStream_t stream);
    /** Downloads the boundary force field and reduces it through the common host owner mapping. */
    void copyForceAndTorqueFromDevice(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream);
    /** Adds scaled held reaction loads without advancing contact history or rigid-body state. */
    void applyForceAndTorque(LSParticleContainer& LSParticles, math::Real scale = 1.0) const;
    /** Clears sampled kinematics while retaining the owner mapping. */
    void clearKinematics() noexcept;

private:
    std::vector<int> ownerOffsets_;           ///< Owner-to-sample offsets with one sentinel.
    std::vector<int> virtualParticleIndices_; ///< Flattened virtual-particle indices.
    std::vector<Vec3> velocitySums_;           ///< Accumulated owner-sample velocities.
    std::vector<Vec3> accelerationSums_;       ///< Accumulated owner-sample accelerations.
    std::vector<Vec3> previousVelocities_;     ///< Previous world sample velocities for prescribed acceleration.
    std::vector<Vec3> forcesByOwner_;          ///< Held world fluid reaction force per LS owner.
    std::vector<Vec3> torquesByOwner_;         ///< Held world fluid reaction torque per LS owner.
};

} // namespace fundem
