/**
 * @file SPHContinuationState.h
 * @brief Physical fluid and inlet state at a represented completed DEM time.
 */
#pragma once

#include "particle/SPHJet.h"
#include "particle/SPHParticle.h"

#include <vector>

namespace fundem
{

/** One finite source's supply and release progress, independent of a search grid. */
struct SPHJetContinuation {
    SPHJet value;
    math::Real startTime{0.0};
    int axialCount{0};
    int bufferLayerCount{0};
    int firstInletLayer{0};
    std::vector<math::Vec3> offsets;
    std::vector<int> particleIndices;
};

/**
 * Accepted represented-time state for a new acoustic interval. A selected
 * subsystem cannot retain a future endpoint computed with excluded neighbours.
 * Restore preserves its physical fields and inlet progress, rebuilding current
 * loads/searches without claiming an uninterrupted leading-interval restart.
 */
struct SPHContinuationState {
    math::Real representedTime{0.0};
    std::vector<SPHParticle> particles;
    std::vector<SPHJetContinuation> jets;
};

} // namespace fundem
