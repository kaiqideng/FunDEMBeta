/**
 * @file sphInteractionKernel.cuh
 * @brief Declares streamed CUDA launch interfaces for WCSPH interactions.
 */
#pragma once

#include "particle/LSParticle.h"
#include "particle/SPHParticle.h"
#include "particle/spatialGrid.h"
#include "particle/virtualParticle.h"
#include "interaction/SPHNeighborhood.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{

/** Reduces actual fluid and boundary displacement from the cached grid reference. */
math::Real maximumSPHNeighborDisplacementSquared(const SPHParticleContainer& particles,
                                                 const virtualParticleContainer& boundaries,
                                                 const SPHNeighborPositionContainer& reference,
                                                 cudaStream_t stream);

/** Device-reduced fluid extrema and invalid-state count. */
struct SPHStateStatistics {
    math::Real maximumVelocity_{0.0};     ///< Largest finite fluid speed.
    math::Real maximumAcceleration_{0.0}; ///< Largest finite acceleration.
    math::Real minimumDensity_{0.0};      ///< Smallest finite positive density.
    math::Real maximumDensity_{0.0};      ///< Largest finite density.
    int invalidValueCount_{0};            ///< Number of invalid particles.
};

/** Device-reduced virtual-boundary kinematic extrema. */
struct SPHKinematicsStatistics {
    math::Real maximumVelocity_{0.0};     ///< Largest finite boundary speed.
    math::Real maximumAcceleration_{0.0}; ///< Largest finite boundary acceleration.
    int invalidValueCount_{0};            ///< Number of invalid boundary samples.
};

/** Initializes all SPH particle state asynchronously on @p stream. */
void launchInitializeSPHParticles(SPHParticleContainer& particles, math::Real referenceDensity, math::Real soundSpeed, cudaStream_t stream = nullptr);
/** Updates boundary sample positions and normals from their owning LS bodies. */
void launchUpdateVirtualParticlePositionAndNormal(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream = nullptr);
/** Clears sampled velocity, acceleration, force, and normal accumulators. */
void launchClearVirtualParticleKinematics(virtualParticleContainer& virtualParticles, cudaStream_t stream = nullptr);
/** Accumulates owner-body kinematics at every virtual boundary sample. */
void launchAccumulateVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const math::Vec3& gravity, math::Real timeStep, cudaStream_t stream = nullptr);
/** Compatibility overload for instantaneous initialization without a difference interval. */
inline void launchAccumulateVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const math::Vec3& gravity, cudaStream_t stream = nullptr)
{
    launchAccumulateVirtualParticleKinematics(virtualParticles, LSParticles, gravity, math::Real{0.0}, stream);
}
/** Resets the prescribed-motion difference after initial external motion hooks have run. */
void launchResetVirtualParticleKinematicsReference(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream = nullptr);
/** Converts accumulated boundary kinematics to the interval average. */
void launchAverageVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, math::Real inverseSampleCount, cudaStream_t stream = nullptr);
/** Atomically scatters sampled fluid loads to their owning LS bodies. */
void launchAddVirtualParticleForceAndTorque(LSParticleContainer& LSParticles, const virtualParticleContainer& virtualParticles, cudaStream_t stream = nullptr);
/** Retains the current world force before a fluid update for impulse matching. */
void launchCacheVirtualParticleCouplingForce(virtualParticleContainer& virtualParticles, cudaStream_t stream = nullptr);
/** Refreshes held world torques at the fluid evaluation pose and records load increments. */
void launchRefreshVirtualParticleCouplingLoads(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream = nullptr);
/** Applies endpoint load-impulse corrections after the final DEM velocity half kick. */
void launchApplyVirtualParticleImpulseCorrection(LSParticleContainer& LSParticles, const virtualParticleContainer& virtualParticles, math::Real correctionTime, cudaStream_t stream = nullptr);

/** Detects free-surface SPH particles using fluid and boundary neighborhoods. */
void launchUpdateSPHFreeSurface(SPHParticleContainer& particles,
                                const virtualParticleContainer& virtualParticles,
                                const spatialGridContainer::const_device_type& particleSpatialGrid,
                                const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                math::Real smoothingLength,
                                math::Real neighborSearchRadius,
                                math::Real particleMass,
                                cudaStream_t stream = nullptr);
/** Applies kernel-summation density reinitialization and updates pressure. */
void launchReinitializeSPHDensity(SPHParticleContainer& particles,
                                  const virtualParticleContainer& virtualParticles,
                                  const spatialGridContainer::const_device_type& particleSpatialGrid,
                                  const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                  math::Real smoothingLength,
                                  math::Real neighborSearchRadius,
                                  math::Real particleMass,
                                  math::Real referenceDensity,
                                  math::Real latticeKernelSum3D,
                                  math::Real soundSpeed,
                                  cudaStream_t stream = nullptr);
/** Advances density with the Riemann continuity rate and refreshes pressure. */
void launchUpdateSPHDensity(SPHParticleContainer& particles,
                            const virtualParticleContainer& virtualParticles,
                            const spatialGridContainer::const_device_type& particleSpatialGrid,
                            const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                            math::Real smoothingLength,
                            math::Real neighborSearchRadius,
                            math::Real particleMass,
                            math::Real referenceDensity,
                            math::Real soundSpeed,
                            math::Real timeStep,
                            const math::Vec3& gravity,
                            cudaStream_t stream = nullptr);
/** Evaluates viscous/prior acceleration and the matching boundary reaction. */
void launchUpdateSPHPriorForceAndBoundaryForce(SPHParticleContainer& particles,
                                               virtualParticleContainer& virtualParticles,
                                               const spatialGridContainer::const_device_type& particleSpatialGrid,
                                               const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                               math::Real smoothingLength,
                                               math::Real neighborSearchRadius,
                                               math::Real particleMass,
                                               math::Real dynamicViscosity,
                                               cudaStream_t stream = nullptr);
/** Evaluates pressure acceleration and the matching boundary reaction. */
void launchUpdateSPHPressureForceAndBoundaryForce(SPHParticleContainer& particles,
                                                  virtualParticleContainer& virtualParticles,
                                                  const spatialGridContainer::const_device_type& particleSpatialGrid,
                                                  const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                                  math::Real smoothingLength,
                                                  math::Real neighborSearchRadius,
                                                  math::Real particleMass,
                                                  math::Real referenceDensity,
                                                  math::Real soundSpeed,
                                                  const math::Vec3& gravity,
                                                  cudaStream_t stream = nullptr);
/** Integrates SPH velocity using the current split-force state. */
void launchIntegrateSPHVelocity(SPHParticleContainer& particles, const math::Vec3& gravity, math::Real timeStep, cudaStream_t stream = nullptr);
/** Integrates SPH particle positions on @p stream. */
void launchIntegrateSPHPosition(SPHParticleContainer& particles, math::Real timeStep, cudaStream_t stream = nullptr);
/** Reduces current fluid statistics and synchronizes only the returned host result. */
SPHStateStatistics calculateSPHStateStatistics(const SPHParticleContainer& particles, const math::Vec3& gravity, cudaStream_t stream = nullptr);
/** Reduces current boundary statistics and synchronizes only the returned host result. */
SPHKinematicsStatistics calculateSPHKinematicsStatistics(const virtualParticleContainer& particles, cudaStream_t stream = nullptr);

} // namespace fundem::cuda
