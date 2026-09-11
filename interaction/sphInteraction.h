/**
 * @file sphInteraction.h
 * @brief Declares CPU WCSPH interaction stages and compact neighbor storage.
 */
#pragma once

#include "particle/LSParticle.h"
#include "particle/SPHParticle.h"
#include "particle/virtualParticle.h"

#include <vector>

namespace fundem::cpu
{

/** Fluid-state extrema and validity counters produced by a host reduction. */
struct SPHStateStatistics {
    math::Real maximumVelocity_{0.0};     ///< Largest finite fluid speed.
    math::Real maximumAcceleration_{0.0}; ///< Largest finite acceleration including gravity.
    math::Real minimumDensity_{0.0};      ///< Smallest finite positive density.
    math::Real maximumDensity_{0.0};      ///< Largest finite density.
    int invalidValueCount_{0};            ///< Number of particles containing invalid state.
};

/** Virtual-boundary kinematic extrema and validity counters produced by a host reduction. */
struct SPHKinematicsStatistics {
    math::Real maximumVelocity_{0.0};     ///< Largest finite boundary speed.
    math::Real maximumAcceleration_{0.0}; ///< Largest finite boundary acceleration.
    int invalidValueCount_{0};            ///< Number of invalid boundary samples.
};

/** Immutable physical values shared by all CPU interaction loops in one solver stage. */
struct SPHInteractionParameters {
    math::Real smoothingLength_{0.0};        ///< Kernel smoothing length for this stage.
    math::Real particleMass_{0.0};           ///< Uniform fluid-particle mass.
    math::Real referenceDensity_{0.0};       ///< Equation-of-state reference density.
    math::Real latticeKernelSum3D_{0.0};     ///< Reference density reinitialization normalization.
    math::Real soundSpeed_{0.0};             ///< Artificial speed of sound.
    math::Real dynamicViscosity_{0.0};       ///< Uniform dynamic viscosity.
    math::Vec3 gravity_{math::Vec3::zero()}; ///< Uniform gravitational acceleration.
};

/** Builds compact host neighborhoods and evaluates the CPU WCSPH stages. */
class SPHInteraction
{
public:
    /** Initializes density, pressure, force history, and finite particle mass. */
    void initializeParticles(SPHParticleContainer& particles, math::Real referenceDensity, math::Real soundSpeed) const;
    /**
     * Rebuilds compact fluid-fluid and fluid-boundary neighborhoods using a
     * uniform host background grid, at advection boundaries or when the search skin expires.
     */
    void buildNeighborhood(const SPHParticleContainer& particles,
                           const virtualParticleContainer& virtualParticles,
                           const math::Vec3& minimumBoundary,
                           const math::Vec3& maximumBoundary,
                           math::Real cellSize,
                           math::Real searchRadius);
    /** Marks free-surface particles from kernel-support deficiency. */
    void updateFreeSurface(SPHParticleContainer& particles, const virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const;
    /** Reinitializes density using the current fixed neighborhood. */
    void reinitializeDensity(SPHParticleContainer& particles, const virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const;
    /** Advances density and pressure by one acoustic substep. */
    void updateDensity(SPHParticleContainer& particles, const virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters, math::Real timeStep) const;
    /** Evaluates viscosity, gravity-independent prior force, and wall reaction. */
    void updatePriorForceAndBoundaryForce(SPHParticleContainer& particles, virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const;
    /** Evaluates pressure force and corresponding wall reaction. */
    void updatePressureForceAndBoundaryForce(SPHParticleContainer& particles, virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const;
    /** Integrates fluid velocity for one acoustic substep. */
    void integrateVelocity(SPHParticleContainer& particles, const math::Vec3& gravity, math::Real timeStep) const;
    /** Integrates fluid position for one acoustic substep. */
    void integratePosition(SPHParticleContainer& particles, math::Real timeStep) const;

    /** Reduces fluid extrema and invalid-state counts without persistent storage. */
    SPHStateStatistics stateStatistics(const SPHParticleContainer& particles, const math::Vec3& gravity) const;
    /** Reduces boundary kinematic extrema and invalid-state counts. */
    SPHKinematicsStatistics kinematicsStatistics(const virtualParticleContainer& particles) const;

private:
    /** Compressed-row neighbor indices for one owner population. */
    struct neighborList {
        const int* begin(int particleIndex) const noexcept { return indices_.empty() ? nullptr : indices_.data() + offsets_[particleIndex]; }
        const int* end(int particleIndex) const noexcept { return indices_.empty() ? nullptr : indices_.data() + offsets_[particleIndex + 1]; }

        std::vector<int> offsets_; ///< Owner-to-neighbor range offsets with one sentinel.
        std::vector<int> indices_; ///< Flattened neighbor indices.
    };

    neighborList particleNeighbors_;        ///< Fluid-fluid neighborhoods.
    neighborList virtualParticleNeighbors_; ///< Fluid-boundary neighborhoods.
};

} // namespace fundem::cpu
