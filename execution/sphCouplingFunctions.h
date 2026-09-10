/**
 * @file sphCouplingFunctions.h
 * @brief Defines endpoint impulse matching for explicitly held SPH wall loads.
 */
#pragma once

#include "math/Quaternion.h"

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

/**
 * Matches the reaction impulse after a completed acoustic interval. The rigid
 * body has already received its final velocity-Verlet half kick. With held
 * loads F0 and refreshed loads F1, correctionTime is acousticDt - DEMDt / 2.
 * This is an impulse at the endpoint: it does not repeat gravity, external
 * loads, or gyroscopic integration, and it does not change the endpoint pose.
 */
FUNDEM_MATH_HD inline void calculateSPHCouplingVelocityCorrection(math::Vec3& velocityCorrection,
                                                                 math::Vec3& angularVelocityCorrection,
                                                                 const math::Vec3& forceIncrement,
                                                                 const math::Vec3& torqueIncrement,
                                                                 const math::Quaternion& orientation,
                                                                 math::Real inverseMass,
                                                                 const math::Mat3& inverseInertiaTensor,
                                                                 math::Real correctionTime) noexcept
{
    velocityCorrection = math::Vec3::zero();
    angularVelocityCorrection = math::Vec3::zero();
    if (inverseMass <= 0.0 || correctionTime <= 0.0)
        return;
    velocityCorrection = (correctionTime * inverseMass) * forceIncrement;
    angularVelocityCorrection = correctionTime * math::rotateUnit(orientation, inverseInertiaTensor * math::inverseRotateUnit(orientation, torqueIncrement));
}

} // namespace fundem::execution
