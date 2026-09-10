/**
 * @file SPHNeighborhood.h
 * @brief Tracks actual displacement while SPH search structures are reused.
 */
#pragma once

#include "particle/SPHParticle.h"
#include "particle/virtualParticle.h"

namespace fundem
{

/** Position at the most recent neighborhood build, independent of physical particle state. */
struct SPHNeighborPosition {
    math::Vec3 position_{math::Vec3::zero()}; ///< World position at the last search build.
    struct positionField { inline static constexpr auto member = &SPHNeighborPosition::position_; };
    using DeviceLayout = deviceLayout<positionField>;
};
using SPHNeighborPositionContainer = hostAoSDeviceSoA<SPHNeighborPosition, SPHNeighborPosition::DeviceLayout>;

/** Conservatively bounds relative motion by twice the largest fluid/boundary displacement. */
class SPHNeighborhood
{
public:
    /** Forgets search validity without discarding allocated storage. */
    void invalidate() noexcept { valid_ = false; }
    bool valid() const noexcept { return valid_; }
    const SPHNeighborPositionContainer& positions() const noexcept { return positions_; }
    double deviceMemoryGB() const noexcept { return positions_.deviceMemoryGB(); }
    /** Captures host positions after a neighborhood rebuild. */
    void captureHost(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries);
    /** Captures device positions on the same stream as the grid build. */
    void captureDevice(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries, cudaStream_t stream);
    /** Checks the host skin bound before a stage uses cached neighbors. */
    bool needsRebuild(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries, math::Real skin) const;
    /** Tests a reduced displacement against the relative-motion allowance. */
    bool exceedsSkin(math::Real maximumSquaredDisplacement, math::Real skin) const noexcept
    {
        return !valid_ || !math::isFinite(maximumSquaredDisplacement) || skin < 0.0 || 4.0 * maximumSquaredDisplacement > skin * skin;
    }

private:
    SPHNeighborPositionContainer positions_; ///< Fluid positions followed by boundary positions.
    bool valid_{false};                     ///< Whether reference positions match the cached search.
};

} // namespace fundem
