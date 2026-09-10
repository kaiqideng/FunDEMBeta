/**
 * @file sphSourceDiscretization.h
 * @brief Provides regular cylindrical SPH-source discretization.
 */
#pragma once

#include "math/Vector3.h"

#include <cmath>

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
    result.radialIndex_ = static_cast<int>(std::floor(radialRatio));
    result.axialCount_ = static_cast<int>(std::floor(length / spacing));
    result.radialLimitSquared_ = radialRatio * radialRatio;
    int crossSectionParticleCount = 0;
    for (int y = -result.radialIndex_; y <= result.radialIndex_; ++y)
    {
        for (int x = -result.radialIndex_; x <= result.radialIndex_; ++x)
        {
            if (result.contains(x, y))
                ++crossSectionParticleCount;
        }
    }
    result.particleCount_ = crossSectionParticleCount * result.axialCount_;
    return true;
}

} // namespace fundem::execution
