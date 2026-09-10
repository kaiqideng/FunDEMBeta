/**
 * @file structureIceLoadKernel.cuh
 * @brief Declares streamed CUDA launchers for ice hydrodynamics and structure-load accumulation.
 */
#pragma once

#include "particle/LSParticle.h"
#include "particle/particle.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{

/** Adds buoyancy and water drag to every spherical ice element. */
void launchSphereBuoyancyAndDrag(particle::device_type particles,
                                 const math::Vec3& gravity,
                                 const math::Vec3& waterVelocity,
                                 math::Real waterDensity,
                                 math::Real waterLevel,
                                 math::Real dragCoefficient,
                                 cudaStream_t stream);

/** Adds each LS structure force to its matching interval accumulator. */
void launchStructureForceAccumulation(math::Vec3* accumulatedForces, LSParticle::device_type structures, cudaStream_t stream);

/** Clears one structure-force accumulator on the supplied stream. */
void launchStructureForceClear(math::Vec3* accumulatedForces, int structureCount, cudaStream_t stream);

} // namespace fundem::cuda
