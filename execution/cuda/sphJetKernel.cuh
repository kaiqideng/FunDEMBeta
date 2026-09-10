/**
 * @file sphJetKernel.cuh
 * @brief Declares the streamed CUDA launch interface for SPH jet velocity constraints.
 */
#pragma once

#include "particle/SPHParticle.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{

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
                               cudaStream_t stream = nullptr);

} // namespace fundem::cuda
