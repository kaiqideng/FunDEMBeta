/**
 * @file parallelPolicy.h
 * @brief Shares work-size thresholds for OpenMP particle and interaction loops.
 */
#pragma once

namespace fundem::cpu
{

constexpr int parallelParticleThreshold = 1024;   ///< Minimum particle count for parallel bulk loops.
constexpr int parallelInteractionThreshold = 256; ///< Minimum interaction count for staged parallel assembly.
constexpr int parallelSurfaceNodeThreshold = 4096; ///< LS narrow-phase work threshold, counted across candidate pairs.
constexpr int surfaceNodeWorkBlockSize = 1024;     ///< Maximum nodes per independently scheduled LS query block.

} // namespace fundem::cpu
