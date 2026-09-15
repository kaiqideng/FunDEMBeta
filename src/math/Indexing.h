/**
 * @file Indexing.h
 * @brief Provides compact multidimensional-to-linear index conversion helpers.
 */
#pragma once

#include "Constants.h"

namespace fundem::math
{

/** Converts a 3D grid coordinate to an x-major linear index. */
FUNDEM_MATH_HD constexpr int linearIndex(int x, int y, int z, int sizeX, int sizeY) noexcept { return (z * sizeY + y) * sizeX + x; }

} // namespace fundem::math
