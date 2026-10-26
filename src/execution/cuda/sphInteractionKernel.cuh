/**
 * @file sphInteractionKernel.cuh
 * @brief Declares streamed CUDA launch interfaces for WCSPH interactions.
 */
#pragma once

#include "interaction/SPHInteractionTypes.h"
#include "particle/LSParticle.h"
#include "particle/SPHParticle.h"
#include "particle/spatialGrid.h"
#include "particle/virtualParticle.h"
#include "interaction/SPHNeighborhood.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{

using SPHStateStatistics = fundem::SPHStateStatistics;
using SPHKinematicsStatistics = fundem::SPHKinematicsStatistics;

/** Reduces actual fluid and boundary displacement from the cached grid reference. */
math::Real maximumSPHNeighborDisplacementSquared(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries, const SPHNeighborPositionContainer& reference, cudaStream_t stream);

/** Initializes all SPH particle state asynchronously on @p stream. */
void launchInitializeSPHParticles(SPHParticleContainer& particles, math::Real referenceDensity, math::Real soundSpeed, cudaStream_t stream = nullptr);
/** Updates boundary sample positions and normals from their owning LS bodies. */
void launchUpdateVirtualParticlePositionAndNormal(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream = nullptr);
/** Clears sampled velocity and acceleration accumulators. */
void launchClearVirtualParticleKinematics(virtualParticleContainer& virtualParticles, cudaStream_t stream = nullptr);
/** Accumulates owner-body kinematics at every virtual boundary sample. */
void launchAccumulateVirtualParticleKinematics(virtualParticleContainer& virtualParticles,
                                               const LSParticleContainer& LSParticles,
                                               const math::Vec3& gravity,
                                               math::Real timeStep,
                                               cudaStream_t stream = nullptr);
/** Compatibility overload for instantaneous initialization without a difference interval. */
inline void launchAccumulateVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const math::Vec3& gravity, cudaStream_t stream = nullptr)
{
    launchAccumulateVirtualParticleKinematics(virtualParticles, LSParticles, gravity, math::Real{0.0}, stream);
}
/** Resets the prescribed-motion difference after initial external motion hooks have run. */
void launchResetVirtualParticleKinematicsReference(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream = nullptr);
/** Refreshes the current boundary pose using completed-interval averages, or instantaneous kinematics when @p sampleCount is zero. */
void launchPrepareVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const math::Vec3& gravity, int sampleCount, cudaStream_t stream = nullptr);
/** Atomically scatters scaled held fluid loads without advancing contact history or rigid-body state. */
void launchAddVirtualParticleForceAndTorque(LSParticleContainer& LSParticles, const virtualParticleContainer& virtualParticles, math::Real scale = 1.0, cudaStream_t stream = nullptr);
/** Refreshes held world torques at the fluid evaluation pose. */
void launchRefreshVirtualParticleCouplingLoads(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream = nullptr);

/** Detects free-surface SPH particles using fluid and boundary neighborhoods. */
void launchUpdateSPHFreeSurface(SPHParticleContainer& particles,
                                const virtualParticleContainer& virtualParticles,
                                const spatialGridContainer::const_device_type& particleSpatialGrid,
                                const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                const SPHInteractionParameters& parameters,
                                math::Real neighborSearchRadius,
                                cudaStream_t stream = nullptr);
/** Applies kernel-summation density reinitialization and updates pressure. */
void launchReinitializeSPHDensity(SPHParticleContainer& particles,
                                  const virtualParticleContainer& virtualParticles,
                                  const spatialGridContainer::const_device_type& particleSpatialGrid,
                                  const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                  const SPHInteractionParameters& parameters,
                                  math::Real neighborSearchRadius,
                                  cudaStream_t stream = nullptr);
/** Advances density with the Riemann continuity rate and refreshes pressure. */
void launchUpdateSPHDensity(SPHParticleContainer& particles,
                            const virtualParticleContainer& virtualParticles,
                            const spatialGridContainer::const_device_type& particleSpatialGrid,
                            const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                            const SPHInteractionParameters& parameters,
                            math::Real neighborSearchRadius,
                            math::Real timeStep,
                            cudaStream_t stream = nullptr);
/** Evaluates viscous/prior acceleration and the matching boundary reaction. */
void launchUpdateSPHPriorForceAndBoundaryForce(SPHParticleContainer& particles,
                                               virtualParticleContainer& virtualParticles,
                                               const spatialGridContainer::const_device_type& particleSpatialGrid,
                                               const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                               const SPHInteractionParameters& parameters,
                                               math::Real neighborSearchRadius,
                                               cudaStream_t stream = nullptr);
/** Evaluates pressure acceleration and the matching boundary reaction. */
void launchUpdateSPHPressureForceAndBoundaryForce(SPHParticleContainer& particles,
                                                  virtualParticleContainer& virtualParticles,
                                                  const spatialGridContainer::const_device_type& particleSpatialGrid,
                                                  const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                                  const SPHInteractionParameters& parameters,
                                                  math::Real neighborSearchRadius,
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
