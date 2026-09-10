/**
 * @file sphJetFunctions.h
 * @brief Provides shared host/device geometry for SPH jet velocity constraints.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem::execution
{

using math::Real;
using math::Vec3;

/** Returns whether @p position lies in the upstream cylindrical jet pipe. */
FUNDEM_MATH_HD inline bool isInsideSPHJetPipe(const Vec3& position, const Vec3& outletCenter, const Vec3& unitDirection, Real radius, Real pipeLength) noexcept
{
    const Vec3 outletOffset = outletCenter - position;
    const Real axialDistance = math::dot(outletOffset, unitDirection);
    if (axialDistance <= 0.0 || axialDistance > pipeLength)
    {
        return false;
    }

    const Vec3 radialOffset = outletOffset - axialDistance * unitDirection;
    return math::normSquared(radialOffset) <= radius * radius;
}

} // namespace fundem::execution
