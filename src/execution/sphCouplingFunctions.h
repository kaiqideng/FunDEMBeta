/**
 * @file sphCouplingFunctions.h
 * @brief Defines prescribed-boundary kinematics shared by SPH CPU and CUDA coupling.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem::execution
{

/** Samples prescribed motion without adding a second centripetal term. */
FUNDEM_MATH_HD inline void accumulatePrescribedSPHKinematics(math::Vec3& velocitySum,
                                                            math::Vec3& accelerationSum,
                                                            math::Vec3& previousVelocity,
                                                            const math::Vec3& currentVelocity,
                                                            math::Real timeStep) noexcept
{
    velocitySum += currentVelocity;
    if (timeStep > 0.0)
        accelerationSum += (currentVelocity - previousVelocity) / timeStep;
    previousVelocity = currentVelocity;
}

} // namespace fundem::execution
