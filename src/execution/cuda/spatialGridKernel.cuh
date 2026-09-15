/**
 * @file spatialGridKernel.cuh
 * @brief Declares streamed CUDA launch interfaces for spatial-grid construction.
 */
#pragma once

#include "particle/rigidBody.h"
#include "particle/SPHParticle.h"
#include "particle/spatialGrid.h"
#include "particle/virtualParticle.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{
namespace detail
{

/** Builds and sorts one spatial grid without synchronizing @p stream. */
void launchBuildSpatialGrid(spatialGridContainer& spatialGrid, const math::Vec3* particlePositions, int particleCount, cudaStream_t stream);

} // namespace detail

/** Builds a rigid-particle grid from a device-backed particle storage. */
template <class ParticleStorage> void launchBuildSpatialGrid(spatialGridContainer& spatialGrid, const ParticleStorage& particles, cudaStream_t stream = nullptr)
{
    detail::launchBuildSpatialGrid(spatialGrid, particles.template device<rigidBody::positionField>(), static_cast<int>(particles.deviceSize()), stream);
}

/** Builds the SPH background grid on @p stream. */
inline void launchBuildSPHSpatialGrid(spatialGridContainer& spatialGrid, const SPHParticleContainer& particles, cudaStream_t stream = nullptr)
{
    detail::launchBuildSpatialGrid(spatialGrid, particles.device<pointMass::positionField>(), static_cast<int>(particles.deviceSize()), stream);
}

/** Builds the virtual-boundary-particle background grid on @p stream. */
inline void launchBuildVirtualParticleSpatialGrid(spatialGridContainer& spatialGrid, const virtualParticleContainer& particles, cudaStream_t stream = nullptr)
{
    detail::launchBuildSpatialGrid(spatialGrid, particles.device<pointMass::positionField>(), static_cast<int>(particles.deviceSize()), stream);
}

} // namespace fundem::cuda
