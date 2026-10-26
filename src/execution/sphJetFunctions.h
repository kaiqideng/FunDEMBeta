/**
 * @file sphJetFunctions.h
 * @brief Provides shared host/device geometry for SPH jet velocity constraints.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "math/Vector3.h"

namespace fundem::execution
{

using math::Real;
using math::Vec3;

/** Fixed world-space inlet opening, independent of emission completion. */
struct SPHJetOpening {
    /** Finite upstream cylinder; never removes a downstream target wall. */
    struct geometry {
        Vec3 outletCenter_{Vec3::zero()};  ///< World outlet center.
        Vec3 unitDirection_{Vec3::zero()}; ///< Unit downstream axis.
        Real radius_{0.0};                ///< Inlet radius.
        Real pipeLength_{0.0};            ///< Upstream buffer length.
    };
    geometry value_; ///< Geometry copied unchanged to the device.
    struct valueField { inline static constexpr auto member = &SPHJetOpening::value_; };
    using DeviceLayout = deviceLayout<valueField>;
};
using SPHJetOpeningContainer = hostAoSDeviceSoA<SPHJetOpening>;

/** Tests the finite inlet cylinder expanded by a sampling clearance. */
FUNDEM_MATH_HD inline bool isInsideSPHJetOpening(const Vec3& position, const SPHJetOpening::geometry& opening, Real clearance) noexcept
{
    const Vec3 outletOffset = opening.outletCenter_ - position;
    const Real axialDistance = math::dot(outletOffset, opening.unitDirection_);
    if (axialDistance < -clearance || axialDistance > opening.pipeLength_ + clearance)
        return false;
    const Vec3 radialOffset = outletOffset - axialDistance * opening.unitDirection_;
    const Real radius = opening.radius_ + clearance;
    return math::normSquared(radialOffset) <= radius * radius;
}

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
