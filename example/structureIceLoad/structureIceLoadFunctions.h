/**
 * @file structureIceLoadFunctions.h
 * @brief Defines CPU/GPU-shared buoyancy and drag formulas for spherical ice elements.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem::execution
{

/** Returns the volume of a sphere below a horizontal water surface. */
FUNDEM_MATH_HD inline math::Real sphereSubmergedVolume(math::Real radius, math::Real centerZ, math::Real waterLevel) noexcept
{
    if (radius <= math::defaultTolerance)
    {
        return 0.0;
    }

    const math::Real submergedHeight = waterLevel - (centerZ - radius);
    if (submergedHeight <= 0.0)
    {
        return 0.0;
    }
    if (submergedHeight >= 2.0 * radius)
    {
        return (4.0 / 3.0) * math::pi * radius * radius * radius;
    }
    return math::pi * submergedHeight * submergedHeight * (radius - submergedHeight / 3.0);
}

/** Returns Archimedes buoyancy opposite to the gravity vector. */
FUNDEM_MATH_HD inline math::Vec3 sphereBuoyancyForce(math::Real submergedVolume, math::Real waterDensity, const math::Vec3& gravity) noexcept { return -waterDensity * submergedVolume * gravity; }

/** Returns quadratic translational drag based on the submerged-area ratio. */
FUNDEM_MATH_HD inline math::Vec3 sphereWaterDragForce(const math::Vec3& particleVelocity,
                                                      const math::Vec3& waterVelocity,
                                                      math::Real projectedArea,
                                                      math::Real waterDensity,
                                                      math::Real dragCoefficient) noexcept
{
    const math::Vec3 relativeVelocity = particleVelocity - waterVelocity;
    const math::Real speed = math::norm(relativeVelocity);
    return -0.5 * waterDensity * dragCoefficient * projectedArea * speed * relativeVelocity;
}

/** Returns quadratic rotational drag using the sphere radius as the moment arm. */
FUNDEM_MATH_HD inline math::Vec3 sphereWaterDragTorque(const math::Vec3& angularVelocity, math::Real radius, math::Real projectedArea, math::Real waterDensity, math::Real dragCoefficient) noexcept
{
    const math::Real angularSpeed = math::norm(angularVelocity);
    return -0.5 * waterDensity * dragCoefficient * projectedArea * radius * radius * angularSpeed * angularVelocity;
}

/** Adds buoyancy plus translational and rotational water drag to one spherical ice element. */
FUNDEM_MATH_HD inline void addSphereBuoyancyAndDrag(math::Vec3& force,
                                                    math::Vec3& torque,
                                                    const math::Vec3& position,
                                                    const math::Vec3& velocity,
                                                    const math::Vec3& angularVelocity,
                                                    math::Real radius,
                                                    const math::Vec3& gravity,
                                                    const math::Vec3& waterVelocity,
                                                    math::Real waterDensity,
                                                    math::Real waterLevel,
                                                    math::Real dragCoefficient) noexcept
{
    const math::Real submergedVolume = sphereSubmergedVolume(radius, position.z, waterLevel);
    if (submergedVolume <= 0.0)
    {
        return;
    }

    force += sphereBuoyancyForce(submergedVolume, waterDensity, gravity);
    if (dragCoefficient <= 0.0)
    {
        return;
    }

    const math::Real sphereVolume = (4.0 / 3.0) * math::pi * radius * radius * radius;
    const math::Real projectedArea = submergedVolume / sphereVolume * math::pi * radius * radius;
    force += sphereWaterDragForce(velocity, waterVelocity, projectedArea, waterDensity, dragCoefficient);
    torque += sphereWaterDragTorque(angularVelocity, radius, projectedArea, waterDensity, dragCoefficient);
}

} // namespace fundem::execution
