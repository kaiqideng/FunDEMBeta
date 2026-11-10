/**
 * @file benchmarkIO.h
 * @brief Reads the original ON-DEM coordinates and writes MercuryCG-compatible frames.
 */
#pragma once

#include "periodicShear.h"

#include <filesystem>
#include <vector>

namespace fundem::validation
{

/** Original coordinate bounds and particles in unchanged source-file order. */
struct initialConfiguration {
    Vec3 minimum_{Vec3::zero()};           ///< Lower domain corner in the source coordinates.
    Vec3 maximum_{Vec3::zero()};           ///< Upper domain corner in the source coordinates.
    std::vector<initialSphere> particles_; ///< Wall and mobile spheres with their original IDs.
};

/**
 * Reads the labeled five-column original configuration without sorting its IDs.
 * Benchmark mode requires 3480 unit-diameter spheres and the published domain,
 * with IDs 1--240 lower-wall, 241--480 upper-wall, and 481--3480 mobile.
 * The file contains no density; density is supplied by periodicShear, not inferred.
 * Small test configurations may disable benchmark checks; their wall roles are
 * determined by centers at zMin + radius or zMax - radius.
 */
initialConfiguration readInitialConfiguration(const std::filesystem::path& path, bool requireBenchmark = true);

/**
 * Writes ShearFlow.data.NNNN and ShearFlow.fstat.NNNN without overwriting frames.
 * All spheres, including wall spheres, occupy their unchanged zero-based data
 * row indices. Finite/finite contacts have two oppositely directed fstat rows;
 * finite/wall contacts have one row with the finite endpoint first. Normal and
 * tangent express the force on the first endpoint: force = fn*n + ft*t.
 * Tangential displacement is -||spring|| to match the published example's
 * signed convention; this alone does not establish identical history models.
 */
void writeMercuryFrame(const periodicShear& simulation, const std::filesystem::path& directory, int frame);

} // namespace fundem::validation
