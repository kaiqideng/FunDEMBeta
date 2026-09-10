/**
 * @file virtualParticleCoupling.h
 * @brief Aggregates SPH boundary kinematics and resultant loads by LS owner.
 */
#pragma once

#include "data/CudaTypes.h"
#include "particle/LSParticle.h"
#include "particle/virtualParticle.h"

#include <utility>
#include <vector>

namespace fundem
{

/** Owns host-side virtual-particle kinematics and resultant SPH wall loads. */
class virtualParticleCoupling
{
public:
    using Vec3 = math::Vec3;

    /** Capturable host accumulator state used by deferred SPH snapshots. */
    struct state {
        std::vector<Vec3> velocitySums_;     ///< Accumulated owner-sample velocities.
        std::vector<Vec3> accelerationSums_; ///< Accumulated owner-sample accelerations.
        std::vector<Vec3> previousVelocities_; ///< Previous world sample velocities for prescribed acceleration.
        std::vector<Vec3> forcesByOwner_;    ///< Fluid reaction force per LS owner.
        std::vector<Vec3> torquesByOwner_;   ///< Fluid reaction torque per LS owner.
        std::vector<Vec3> forceIncrementsByOwner_; ///< New minus held world force.
        std::vector<Vec3> torqueIncrementsByOwner_; ///< New minus held world torque.
    };

    /** Releases owner mappings and all accumulator state. */
    void reset() noexcept;
    /** Generates virtual particles and stable owner-to-sample ranges. */
    void initialize(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity);
    /** Samples complete left-endpoint LS kinematics, before force assembly clears the loads. */
    void accumulateKinematics(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const Vec3& gravity, math::Real timeStep = 0.0);
    /** Starts the prescribed-motion difference at the complete initialized endpoint. */
    void resetKinematicsReference(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles);
    /** Writes accumulated average kinematics to virtual particles. */
    void averageKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, int sampleCount);
    /** Reduces virtual-particle reaction forces and torques by LS owner. */
    void collectForceAndTorque(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles);
    /** Averages host kinematics and uploads only the fields required by device SPH. */
    void prepareDeviceKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, int sampleCount, cudaStream_t stream);
    /** Copies device wall forces and reduces them into host owner loads. */
    void copyForceAndTorqueFromDevice(const virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream);
    /** Adds reduced fluid reaction loads to LS rigid bodies. */
    void applyForceAndTorque(LSParticleContainer& LSParticles) const;
    /** Matches the refreshed wall-load impulse after the final DEM half kick; never used by observation-only flushes. */
    void applyImpulseCorrection(LSParticleContainer& LSParticles, math::Real correctionTime) const;
    /** Clears sampled kinematics while retaining the owner mapping. */
    void clearKinematics() noexcept;

    /** Copies all accumulators for a non-destructive multirate snapshot. */
    state captureState() const { return currentState_; }
    /** Restores previously captured accumulators. */
    void restoreState(state value) { currentState_ = std::move(value); }

private:
    std::vector<int> ownerOffsets_;           ///< Owner-to-sample offsets with one sentinel.
    std::vector<int> virtualParticleIndices_; ///< Flattened virtual-particle indices.
    state currentState_;                      ///< Current kinematic and reaction accumulators.
};

} // namespace fundem
