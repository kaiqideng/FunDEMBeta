/**
 * @file PeriodicDomain.h
 * @brief Shared host/device coordinate transforms for a minimum-image periodic cell.
 */
#pragma once

#include "Vector3.h"

namespace fundem::math
{

/** Optional periodic coordinates; callers own the boundary policy and validate interaction ranges. */
struct periodicDomain
{
    Vec3 minimum_{Vec3::zero()}; ///< Lower corner of the fundamental cell.
    Vec3 length_{1.0, 1.0, 1.0}; ///< Positive cell lengths; disabled axes are unchanged.
    bool periodicAxes_[3]{false, false, false}; ///< Independently enabled Cartesian directions.

    FUNDEM_MATH_HD bool enabled() const noexcept { return periodicAxes_[0] || periodicAxes_[1] || periodicAxes_[2]; }

    /** Checks enabled-axis metadata without changing particle coordinates. */
    FUNDEM_MATH_HD bool valid() const noexcept
    {
        for (int axis = 0; axis < 3; ++axis)
            if (periodicAxes_[axis] && (!math::isFinite(minimum_[axis]) || !math::isFinite(length_[axis]) || length_[axis] <= 0.0 || !math::isFinite(minimum_[axis] + length_[axis])))
                return false;
        return true;
    }

    /** Maps a coordinate into [minimum, minimum + length) only in enabled directions. */
    FUNDEM_MATH_HD Vec3 wrapPosition(Vec3 position) const noexcept
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            if (!periodicAxes_[axis])
                continue;
            Real coordinate = ::fmod(position[axis] - minimum_[axis], length_[axis]);
            if (coordinate < 0.0)
                coordinate += length_[axis];
            position[axis] = minimum_[axis] + coordinate;
            if (position[axis] >= minimum_[axis] + length_[axis])
                position[axis] = minimum_[axis];
        }
        return position;
    }

    /** Returns the nearest periodic displacement; physical support must remain smaller than half a period. */
    FUNDEM_MATH_HD Vec3 minimumImage(Vec3 displacement) const noexcept
    {
        for (int axis = 0; axis < 3; ++axis)
            if (periodicAxes_[axis])
                displacement[axis] -= length_[axis] * ::round(displacement[axis] / length_[axis]);
        return displacement;
    }

    /** Translates a query neighbor to the image nearest its reference point. */
    FUNDEM_MATH_HD Vec3 imagePosition(const Vec3& point, const Vec3& reference) const noexcept
    {
        Vec3 image = point;
        for (int axis = 0; axis < 3; ++axis)
            if (periodicAxes_[axis])
                image[axis] -= length_[axis] * ::round((point[axis] - reference[axis]) / length_[axis]);
        return image;
    }
};

} // namespace fundem::math
