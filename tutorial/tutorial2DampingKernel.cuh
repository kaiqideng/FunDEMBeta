/**
 * @file tutorial2DampingKernel.cuh
 * @brief Declares the streamed CUDA damping helper for the bonded cloth tutorial.
 */
#pragma once

#include "particle/particle.h"

namespace fundem::tutorial
{

/** Adds mass-proportional linear damping and isotropic rotational drag to finite spheres. */
void launchClothDamping(particle::device_type particles, math::Real dampingRate, math::Real rotationalDragCoefficient, cudaStream_t stream);

} // namespace fundem::tutorial
