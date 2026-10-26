/**
 * @file sphJetKernel.cuh
 * @brief Declares the streamed CUDA launch interface for SPH jet velocity constraints.
 */
#pragma once

#include "execution/sphJetFunctions.h"
#include "particle/SPHParticle.h"
#include "particle/virtualParticle.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{

/** Updates effective wall volumes and returns whether a mask changed; fences @p stream. */
bool launchUpdateSPHJetOpenings(virtualParticleContainer& virtualParticles,
                                 const execution::SPHJetOpeningContainer& openings,
                                 math::Real nominalVolume,
                                 math::Real clearance,
                                 cudaStream_t stream = nullptr);

/** Applies a prescribed velocity to jet-pipe particles in the selected range. */
void launchApplySPHJetVelocity(SPHParticleContainer& particles,
                               int particleBegin,
                               int particleCount,
                               bool active,
                               const math::Vec3& outletCenter,
                               const math::Vec3& unitDirection,
                               math::Real radius,
                               math::Real pipeLength,
                               const math::Vec3& velocity,
                               math::Real referenceDensity,
                               cudaStream_t stream = nullptr);

} // namespace fundem::cuda
