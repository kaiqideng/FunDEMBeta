/**
 * @file parallelPolicy.h
 * @brief Shares work-size thresholds for OpenMP particle and interaction loops.
 */
#pragma once

namespace fundem::cpu
{

constexpr int parallelParticleThreshold = 1024;   ///< Minimum particle count for parallel bulk loops.
constexpr int parallelInteractionThreshold = 256; ///< Minimum interaction count for staged parallel assembly.

} // namespace fundem::cpu
