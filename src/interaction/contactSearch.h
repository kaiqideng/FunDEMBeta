/**
 * @file contactSearch.h
 * @brief Declares CPU contact searches for sphere and level-set particle pairs.
 */
#pragma once

#include "contact.h"
#include "particle/LSParticle.h"

#include <vector>

namespace fundem
{

/** Ordered candidate pair with master preceding slave. */
class particlePair
{
public:
    int masterParticleIndex_{-1}; ///< Master index in its particle container.
    int slaveParticleIndex_{-1};  ///< Slave index in its particle container.
};

/**
 * CPU broad- and narrow-phase contact search for sphere-sphere, sphere-LS, and
 * LS-LS interactions. Infinite-mass bodies are ordered as slaves.
 */
class contactSearch
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    contactSearch() = default;
    /** Creates a search with a requested grid size; non-positive selects automatic sizing. */
    explicit contactSearch(Real cellSize) noexcept { setCellSize(cellSize); }

    /** Requests a positive grid cell size; zero selects automatic sizing. */
    bool setCellSize(Real value) noexcept;
    void setDirectSearch(bool value) noexcept { directSearch_ = value; }

    Real requestedCellSize() const noexcept { return requestedCellSize_; }
    Real effectiveCellSize() const noexcept { return effectiveCellSize_; }
    bool directSearch() const noexcept { return directSearch_; }
    int candidatePairCount() const noexcept { return static_cast<int>(candidatePairs_.size()); }
    int occupiedCellCount() const noexcept { return static_cast<int>(cells_.size()); }
    const std::vector<particlePair>& candidatePairs() const noexcept { return candidatePairs_; }

    /** Rebuilds sphere-sphere contacts and returns their count. */
    int findContacts(contactContainer& contacts, const particleContainer& spheres);
    /** Rebuilds sphere-LS contacts with LS particles as slaves. */
    int findContacts(contactContainer& contacts, const particleContainer& spheres, const LSParticleContainer& LSParticles);
    /** Rebuilds LS-LS contacts using node-work-based parallel blocks and stable pair/node ordering. */
    int findContacts(contactContainer& contacts, const LSParticleContainer& LSParticles);

    /** Clears cached grid entries, occupied cells, and candidate pairs. */
    void clear() noexcept;

private:
    /** A contiguous part of one master surface; no particle or node indices are reordered. */
    struct surfaceNodeWorkBlock {
        int pairIndex_{0}; ///< Index into the current ordered candidate pairs.
        int nodeBegin_{0}; ///< Inclusive master surface-node index.
        int nodeEnd_{0};   ///< Exclusive master surface-node index.
    };
    /** Integer coordinate of one CPU background-grid cell. */
    struct gridCell {
        int x_{0}; ///< Cell coordinate along x.
        int y_{0}; ///< Cell coordinate along y.
        int z_{0}; ///< Cell coordinate along z.

        bool operator==(const gridCell& other) const noexcept { return x_ == other.x_ && y_ == other.y_ && z_ == other.z_; }
        bool operator<(const gridCell& other) const noexcept { return x_ != other.x_ ? x_ < other.x_ : (y_ != other.y_ ? y_ < other.y_ : z_ < other.z_); }
    };

    /** Particle index paired with its cell before sorting. */
    struct gridParticle {
        gridCell cell_;         ///< Cell containing the particle center.
        int particleIndex_{-1}; ///< Original particle index.
    };

    /** Compact range of sorted particles belonging to one occupied cell. */
    struct occupiedCell {
        gridCell cell_;                   ///< Occupied cell coordinate.
        int sortedParticleIndexBegin_{0}; ///< Inclusive grid-particle range begin.
        int sortedParticleIndexEnd_{0};   ///< Exclusive grid-particle range end.
    };

    /** Clears transient results before one search while retaining configuration. */
    void beginSearch() noexcept;
    /** Chooses the effective grid scale for the current maximum radius. */
    void configureGrid(Real maximumRadius) noexcept;
    /** Builds sorted grid entries and compact occupied-cell ranges. */
    void buildGrid(const std::vector<Vec3>& positions);
    /** Maps one world position to an integer grid coordinate. */
    gridCell calculateCell(const Vec3& position) const noexcept;
    /** Finds a compact occupied-cell record by binary search. */
    const occupiedCell* findCell(const gridCell& cell) const noexcept;
    /** Returns the number of cells needed to cover @p distance. */
    int searchCellRange(Real distance) const noexcept;
    /** Generates same-container broad-phase pairs with stable master/slave order. */
    void findCandidatePairs(const std::vector<Vec3>& positions, const std::vector<Real>& radii);
    /** Generates sphere-master/LS-slave broad-phase pairs. */
    void findSphereLevelSetCandidatePairs(const std::vector<Vec3>& spherePositions,
                                          const std::vector<Real>& sphereRadii,
                                          const std::vector<Vec3>& LSParticlePositions,
                                          const std::vector<Real>& LSParticleRadii);

    Real requestedCellSize_{0.0};              ///< User-requested cell size; zero selects automatic.
    Real effectiveCellSize_{1.0};              ///< Cell size selected for the current search.
    Real inverseCellSize_{1.0};                ///< Reciprocal effective cell size.
    Vec3 gridOrigin_{Vec3::zero()};            ///< Minimum occupied cell origin.
    bool directSearch_{false};                 ///< Bypass the grid and enumerate all pairs.
    std::vector<gridParticle> gridParticles_;  ///< Entries sorted lexicographically by cell.
    std::vector<occupiedCell> cells_;          ///< Compact ranges for occupied cells.
    std::vector<particlePair> candidatePairs_; ///< Broad-phase candidate pairs.
    std::vector<surfaceNodeWorkBlock> surfaceNodeWorkBlocks_; ///< Flattened pair/node work for the current LS search.
    std::vector<std::vector<contact>> surfaceNodeBlockContacts_; ///< Reused block-local buffers, merged in stable order.
};

} // namespace fundem
