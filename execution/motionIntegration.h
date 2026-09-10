/**
 * @file motionIntegration.h
 * @brief Implements host/device point-mass and rigid-body integration formulas.
 */
#pragma once

#include "math/Quaternion.h"

namespace fundem::execution
{

using math::Mat3;
using math::Quaternion;
using math::Real;
using math::Vec3;

/** Returns translational plus rotational kinetic energy; fixed bodies contribute zero. */
FUNDEM_MATH_HD inline Real kineticEnergy(const Vec3& velocity, const Vec3& angularVelocity, const Quaternion& orientation, Real inverseMass, const Mat3& inertiaTensor) noexcept
{
    if (inverseMass <= 0.0)
    {
        return 0.0;
    }
    const Vec3 bodyAngularVelocity = math::inverseRotateUnit(orientation, angularVelocity);
    return 0.5 * (math::normSquared(velocity) / inverseMass + math::dot(bodyAngularVelocity, inertiaTensor * bodyAngularVelocity));
}

/** Returns gravitational potential energy relative to the coordinate origin. */
FUNDEM_MATH_HD inline Real gravitationalPotentialEnergy(const Vec3& position, Real inverseMass, const Vec3& gravity) noexcept
{
    return inverseMass > 0.0 ? math::dot(-gravity, position) / inverseMass : 0.0;
}

/** Advances linear velocity with the accumulated force and gravity. */
FUNDEM_MATH_HD inline void integrateVelocity(Vec3& velocity, const Vec3& force, Real inverseMass, const Vec3& gravity, Real timeStep) noexcept
{
    if (inverseMass > 0.0)
    {
        velocity += timeStep * (inverseMass * force + gravity);
    }
}

/** Advances world-frame angular velocity using the Euler rigid-body equation. */
FUNDEM_MATH_HD inline void integrateAngularVelocity(Vec3& angularVelocity,
                                                    const Vec3& torque,
                                                    const Quaternion& orientation,
                                                    const Mat3& inertiaTensor,
                                                    const Mat3& inverseInertiaTensor,
                                                    Real timeStep) noexcept
{
    const Vec3 angularMomentum = math::rotateUnit(orientation, inertiaTensor * math::inverseRotateUnit(orientation, angularVelocity));
    const Vec3 angularAcceleration = math::rotateUnit(orientation, inverseInertiaTensor * math::inverseRotateUnit(orientation, torque - math::cross(angularVelocity, angularMomentum)));
    angularVelocity += timeStep * angularAcceleration;
}

/** Advances position by one explicit velocity step. */
FUNDEM_MATH_HD inline void integratePosition(Vec3& position, const Vec3& velocity, Real timeStep) noexcept { position += timeStep * velocity; }

/** Advances and renormalizes orientation from world-frame angular velocity. */
FUNDEM_MATH_HD inline void integrateQuaternion(Quaternion& orientation, const Vec3& angularVelocity, Real timeStep) noexcept
{
    Quaternion updatedOrientation;
    if (math::tryIntegrateWorldAngularVelocity(orientation, angularVelocity, timeStep, updatedOrientation))
    {
        orientation = updatedOrientation;
    }
}

} // namespace fundem::execution
