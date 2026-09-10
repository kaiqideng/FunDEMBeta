/**
 * @file contactDetection.h
 * @brief Provides geometry and effective-property helpers shared by contact searches.
 */
#pragma once

#include "math/Indexing.h"
#include "math/Quaternion.h"
#include "data/CudaTypes.h"

namespace fundem::execution
{

using math::Quaternion;
using math::Real;
using math::Vec3;

/** Returns the pair effective mass from two inverse masses. */
FUNDEM_MATH_HD inline Real effectiveMass(Real firstInverseMass, Real secondInverseMass) noexcept
{
    const Real inverseMassSum = firstInverseMass + secondInverseMass;
    return inverseMassSum > 0.0 ? 1.0 / inverseMassSum : 0.0;
}

/** Returns the reduced radius of two positive radii. */
FUNDEM_MATH_HD inline Real effectiveRadius(Real firstRadius, Real secondRadius) noexcept
{
    const Real radiusSum = firstRadius + secondRadius;
    return radiusSum > 0.0 ? firstRadius * secondRadius / radiusSum : 0.0;
}

/** Estimates circular sphere-contact area as pi times effective radius and overlap. */
FUNDEM_MATH_HD inline Real sphereContactArea(Real contactEffectiveRadius, Real overlap) noexcept
{
    return contactEffectiveRadius > 0.0 && overlap > 0.0 ? math::pi * contactEffectiveRadius * overlap : 0.0;
}

/** Tests two bounding spheres without computing contact geometry. */
FUNDEM_MATH_HD inline bool boundingSpheresOverlap(const Vec3& firstPosition, Real firstRadius, const Vec3& secondPosition, Real secondRadius) noexcept
{
    const Real radiusSum = firstRadius + secondRadius;
    return math::distanceSquared(firstPosition, secondPosition) <= radiusSum * radiusSum;
}

/** Computes the clamped 3-by-3-by-3 cell range surrounding a grid position. */
FUNDEM_MATH_HD inline void clampNeighborSpatialGridRange(int3& spatialGridStart, int3& spatialGridEnd, const int3& centralSpatialGridPosition, const int3& spatialGridSize) noexcept
{
    spatialGridStart.x = centralSpatialGridPosition.x > 0 ? centralSpatialGridPosition.x - 1 : 0;
    spatialGridStart.y = centralSpatialGridPosition.y > 0 ? centralSpatialGridPosition.y - 1 : 0;
    spatialGridStart.z = centralSpatialGridPosition.z > 0 ? centralSpatialGridPosition.z - 1 : 0;
    spatialGridEnd.x = centralSpatialGridPosition.x < spatialGridSize.x - 1 ? centralSpatialGridPosition.x + 1 : spatialGridSize.x - 1;
    spatialGridEnd.y = centralSpatialGridPosition.y < spatialGridSize.y - 1 ? centralSpatialGridPosition.y + 1 : spatialGridSize.y - 1;
    spatialGridEnd.z = centralSpatialGridPosition.z < spatialGridSize.z - 1 ? centralSpatialGridPosition.z + 1 : spatialGridSize.z - 1;
}

/** Trilinearly interpolates a signed-distance value within one grid cell. */
FUNDEM_MATH_HD inline Real interpolateLevelSetValue(const Vec3& interpolation, Real phi000, Real phi100, Real phi010, Real phi110, Real phi001, Real phi101, Real phi011, Real phi111) noexcept
{
    const Real wx0 = 1.0 - interpolation.x;
    const Real wy0 = 1.0 - interpolation.y;
    const Real wz0 = 1.0 - interpolation.z;
    return phi000 * wx0 * wy0 * wz0 + phi100 * interpolation.x * wy0 * wz0 + phi010 * wx0 * interpolation.y * wz0 + phi110 * interpolation.x * interpolation.y * wz0 +
           phi001 * wx0 * wy0 * interpolation.z + phi101 * interpolation.x * wy0 * interpolation.z + phi011 * wx0 * interpolation.y * interpolation.z +
           phi111 * interpolation.x * interpolation.y * interpolation.z;
}

/** Differentiates the trilinear interpolant and returns its local-space gradient. */
FUNDEM_MATH_HD inline Vec3 interpolateLevelSetGradient(const Vec3& interpolation,
                                                       Real gridNodeInverseSpacing,
                                                       Real phi000,
                                                       Real phi100,
                                                       Real phi010,
                                                       Real phi110,
                                                       Real phi001,
                                                       Real phi101,
                                                       Real phi011,
                                                       Real phi111) noexcept
{
    const Real wx0 = 1.0 - interpolation.x;
    const Real wx1 = interpolation.x;
    const Real wy0 = 1.0 - interpolation.y;
    const Real wy1 = interpolation.y;
    const Real wz0 = 1.0 - interpolation.z;
    const Real wz1 = interpolation.z;
    const Real gradientX = (phi100 - phi000) * wy0 * wz0 + (phi110 - phi010) * wy1 * wz0 + (phi101 - phi001) * wy0 * wz1 + (phi111 - phi011) * wy1 * wz1;
    const Real gradientY = (phi010 - phi000) * wx0 * wz0 + (phi110 - phi100) * wx1 * wz0 + (phi011 - phi001) * wx0 * wz1 + (phi111 - phi101) * wx1 * wz1;
    const Real gradientZ = (phi001 - phi000) * wx0 * wy0 + (phi101 - phi100) * wx1 * wy0 + (phi011 - phi010) * wx0 * wy1 + (phi111 - phi110) * wx1 * wy1;
    return gridNodeInverseSpacing * Vec3{gradientX, gradientY, gradientZ};
}

/** Reads a signed distance from a device-style scalar grid. */
FUNDEM_MATH_HD constexpr Real levelSetSignedDistance(Real value) noexcept { return value; }

/** Reads a signed distance from a host grid-node object. */
template <typename GridNodeValue> FUNDEM_MATH_HD constexpr Real levelSetSignedDistance(const GridNodeValue& value) noexcept { return value.signedDistance_; }

/**
 * Queries a point or sphere against a level-set grid.
 * @return True when the query overlaps the represented solid and a finite world-space normal is available.
 */
template <typename GridNodeValue>
FUNDEM_MATH_HD inline bool detectLevelSetContact(Real& overlap,
                                                 Vec3& normal,
                                                 const GridNodeValue* gridNodes,
                                                 const Vec3& queryPosition,
                                                 const Vec3& levelSetPosition,
                                                 const Quaternion& levelSetOrientation,
                                                 const Vec3& gridNodeOrigin,
                                                 Real gridNodeInverseSpacing,
                                                 const int3& gridNodeSize,
                                                 int signedDistanceOffset,
                                                 Real queryRadius = 0.0) noexcept
{
    overlap = 0.0;
    normal = Vec3::zero();

    const Vec3 queryLocalPosition = math::inverseRotateUnit(levelSetOrientation, queryPosition - levelSetPosition);
    const Vec3 levelSetGridPosition = gridNodeInverseSpacing * (queryLocalPosition - gridNodeOrigin);
    const int x0 = static_cast<int>(math::detail::floor(levelSetGridPosition.x));
    const int y0 = static_cast<int>(math::detail::floor(levelSetGridPosition.y));
    const int z0 = static_cast<int>(math::detail::floor(levelSetGridPosition.z));
    if (x0 < 0 || y0 < 0 || z0 < 0 || x0 >= gridNodeSize.x - 1 || y0 >= gridNodeSize.y - 1 || z0 >= gridNodeSize.z - 1)
    {
        return false;
    }

    const int x1 = x0 + 1;
    const int y1 = y0 + 1;
    const int z1 = z0 + 1;
    const Vec3 interpolation{levelSetGridPosition.x - Real(x0), levelSetGridPosition.y - Real(y0), levelSetGridPosition.z - Real(z0)};
    const int sizeX = gridNodeSize.x;
    const int sizeY = gridNodeSize.y;
    const Real phi000 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x0, y0, z0, sizeX, sizeY)]);
    const Real phi100 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x1, y0, z0, sizeX, sizeY)]);
    const Real phi010 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x0, y1, z0, sizeX, sizeY)]);
    const Real phi110 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x1, y1, z0, sizeX, sizeY)]);
    const Real phi001 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x0, y0, z1, sizeX, sizeY)]);
    const Real phi101 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x1, y0, z1, sizeX, sizeY)]);
    const Real phi011 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x0, y1, z1, sizeX, sizeY)]);
    const Real phi111 = levelSetSignedDistance(gridNodes[signedDistanceOffset + math::linearIndex(x1, y1, z1, sizeX, sizeY)]);

    overlap = queryRadius - interpolateLevelSetValue(interpolation, phi000, phi100, phi010, phi110, phi001, phi101, phi011, phi111);
    if (overlap < 0.0)
    {
        overlap = 0.0;
        return false;
    }

    normal = math::rotateUnit(levelSetOrientation, interpolateLevelSetGradient(interpolation, gridNodeInverseSpacing, phi000, phi100, phi010, phi110, phi001, phi101, phi011, phi111));
    if (!math::tryNormalize(normal))
    {
        overlap = 0.0;
        normal = Vec3::zero();
        return false;
    }
    return true;
}

/** Places a surface-node/level-set contact halfway through the penetration segment. */
FUNDEM_MATH_HD inline Vec3 surfaceNodeLevelSetContactPoint(const Vec3& surfaceNodePosition, const Vec3& normal, Real overlap) noexcept { return surfaceNodePosition + 0.5 * overlap * normal; }

/** Returns the sphere-side contact point for a sphere/level-set overlap. */
FUNDEM_MATH_HD inline Vec3 sphereLevelSetContactPoint(const Vec3& spherePosition, Real sphereRadius, const Vec3& normal, Real overlap) noexcept
{
    return spherePosition - (sphereRadius - overlap) * normal;
}

/** Computes sphere/sphere contact point, slave-to-master normal, overlap, area, and effective radius. */
FUNDEM_MATH_HD inline bool detectSphereContact(Vec3& point,
                                               Vec3& normal,
                                               Real& overlap,
                                               Real& area,
                                               Real& contactEffectiveRadius,
                                               const Vec3& masterPosition,
                                               Real masterRadius,
                                               const Vec3& slavePosition,
                                               Real slaveRadius) noexcept
{
    point = Vec3::zero();
    normal = Vec3::zero();
    overlap = 0.0;
    area = 0.0;
    contactEffectiveRadius = 0.0;

    const Vec3 centerDifference = masterPosition - slavePosition;
    const Real centerDistanceSquared = math::normSquared(centerDifference);
    const Real radiusSum = masterRadius + slaveRadius;
    if (radiusSum <= 0.0 || centerDistanceSquared > radiusSum * radiusSum || centerDistanceSquared <= 0.0)
    {
        return false;
    }

    const Real centerDistance = math::detail::sqrt(centerDistanceSquared);
    normal = centerDifference / centerDistance;
    overlap = radiusSum - centerDistance;
    point = slavePosition + (slaveRadius - 0.5 * overlap) * normal;
    contactEffectiveRadius = effectiveRadius(masterRadius, slaveRadius);
    area = sphereContactArea(contactEffectiveRadius, overlap);
    return true;
}

} // namespace fundem::execution
