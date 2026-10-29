/**
 * @file sphSpatialGridFunctions.h
 * @brief Unique periodic cell ranges shared by host and device SPH searches.
 */
#pragma once

#include "math/PeriodicDomain.h"

namespace fundem::execution
{

/** Maps a visited cell back into its periodic axis, without duplicating cells. */
FUNDEM_MATH_HD inline int SPHGridCell(int coordinate, int size, bool periodic) noexcept
{
    if (!periodic)
        return coordinate;
    const int cell = coordinate % size;
    return cell < 0 ? cell + size : cell;
}

/** Returns a finite cell span; a wide periodic query visits the entire axis once. */
FUNDEM_MATH_HD inline void SPHGridAxisRange(int& begin, int& end, math::Real coordinate, math::Real extent, int size, bool periodic) noexcept
{
    if (periodic)
    {
        if (2.0 * extent + 1.0 >= size)
        {
            begin = 0;
            end = size - 1;
            return;
        }
        begin = static_cast<int>(math::detail::floor(coordinate - extent));
        end = static_cast<int>(math::detail::floor(coordinate + extent));
        if (end - begin >= size)
            end = begin + size - 1;
        return;
    }
    const math::Real lower = coordinate - extent, upper = coordinate + extent;
    begin = lower <= 0.0 ? 0 : (lower >= size ? size - 1 : static_cast<int>(math::detail::floor(lower)));
    end = upper <= 0.0 ? 0 : (upper >= size ? size - 1 : static_cast<int>(math::detail::floor(upper)));
}

} // namespace fundem::execution
