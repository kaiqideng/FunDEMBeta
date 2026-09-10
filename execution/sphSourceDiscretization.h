/**
 * @file sphSourceDiscretization.h
 * @brief Provides regular cylindrical SPH-source discretization.
 */
#pragma once

#include "math/Vector3.h"

#include <cmath>
#include <limits>

namespace fundem::execution
{

using math::Real;

/** Integer lattice dimensions and particle count for one cylindrical SPH jet. */
struct SPHJetDiscretization {
    int radialIndex_{0};
    int axialCount_{0};
    int particleCount_{0};
    Real radialLimitSquared_{0.0};

    /** Tests one dimensionless transverse lattice coordinate against the circular inlet. */
    bool contains(int firstIndex, int secondIndex) const noexcept
    {
        const Real first = static_cast<Real>(firstIndex);
        const Real second = static_cast<Real>(secondIndex);
        return first * first + second * second <= radialLimitSquared_;
    }
};

/** Computes the regular lattice used to populate one cylindrical SPH source. */
inline bool trySPHJetDiscretization(Real spacing, Real radius, Real length, SPHJetDiscretization& result) noexcept
{
    result = {};
    if (!math::isFinite(spacing) || spacing <= 0.0 || !math::isFinite(radius) || radius < 0.5 * spacing || !math::isFinite(length) || length < spacing)
        return false;

    const Real radialRatio = (radius - 0.5 * spacing) / spacing;
    const Real axialCount = std::floor(length / spacing);
    const int maximumCount = std::numeric_limits<int>::max();
    // A valid lattice must fit its int-indexed particle container before any cast or loop.
    if (!math::isFinite(radialRatio) || !math::isFinite(axialCount) || axialCount > maximumCount || radialRatio > (maximumCount - 1) / 2)
        return false;
    const Real inscribedHalfWidth = std::floor(radialRatio / std::sqrt(2.0));
    if ((2.0 * inscribedHalfWidth + 1.0) * (2.0 * inscribedHalfWidth + 1.0) > maximumCount / axialCount)
        return false;
    result.radialIndex_ = static_cast<int>(std::floor(radialRatio));
    result.axialCount_ = static_cast<int>(axialCount);
    result.radialLimitSquared_ = radialRatio * radialRatio;
    int crossSectionParticleCount = 0;
    for (int y = -result.radialIndex_; y <= result.radialIndex_; ++y)
    {
        // Count one circular lattice row in O(1), correcting sqrt rounding with
        // the same predicate used by particle generation.
        const Real extentSquared = result.radialLimitSquared_ - static_cast<Real>(y) * y;
        int maximumX = static_cast<int>(std::floor(std::sqrt(extentSquared > 0.0 ? extentSquared : 0.0)));
        if (maximumX > result.radialIndex_)
            maximumX = result.radialIndex_;
        while (maximumX >= 0 && !result.contains(maximumX, y))
            --maximumX;
        while (maximumX < result.radialIndex_ && result.contains(maximumX + 1, y))
            ++maximumX;
        const int rowCount = 2 * maximumX + 1;
        if (crossSectionParticleCount > maximumCount / result.axialCount_ - rowCount)
        {
            result = {};
            return false;
        }
        crossSectionParticleCount += rowCount;
    }
    result.particleCount_ = crossSectionParticleCount * result.axialCount_;
    return true;
}

} // namespace fundem::execution
