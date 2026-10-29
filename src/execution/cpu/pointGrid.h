/**
 * @file pointGrid.h
 * @brief Declares reusable host occupied-cell queries and fixed-radius point pairs.
 */
#pragma once

#include "math/Vector3.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace fundem::cpu
{

/** Local indices in one point array or in the first/second arrays of a cross query.
 * A point pair has no particle radius, contact, material, or mass semantics.
 */
struct pointPair {
    int firstIndex_{-1};  ///< Index in the first input point array.
    int secondIndex_{-1}; ///< Index in the second input array (the same array for a self query).
};

/** Sorted occupied-cell index shared by CPU contact and center-distance searches.
 * The index owns its cell/index arrays, not the supplied positions. Concurrent
 * const queries are safe when build()/clear() are not running. Coordinates are
 * discretized with checked signed 64-bit cells. Outward-rounded world bounds
 * protect queries against origin-subtraction rounding. Unrepresentable cells or
 * query extents fall back to all indexed points rather than overflow or lose neighbors.
 */
class pointGrid
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    /** Rebuilds from finite points with finite positive cell size.
     * @throws std::invalid_argument On nonfinite points or nonpositive/nonfinite size.
     * @throws std::length_error When local point indices cannot fit in int.
     * The original point array and its ordering remain unchanged.
     */
    void build(const std::vector<Vec3>& positions, Real cellSize);

    /** Clears indexed entries while retaining allocated capacity. */
    void clear() noexcept;

    std::size_t occupiedCellCount() const noexcept { return cells_.size(); }

    /** Visits local indices in a conservative radius neighborhood, not exact pairs.
     * Normal neighborhoods retain the CPU contact search's z/y/x cell traversal
     * and increasing original index within each cell. Large sparse neighborhoods
     * scan only occupied cells in that same order. The visitor performs its own
     * exact distance/contact test. A zero radius is valid; invalid/nonfinite
     * query data throws std::invalid_argument. No target positions are reread.
     */
    template <class Visitor> void forEachCandidate(const Vec3& position, Real radius, Visitor&& visitor) const
    {
        const auto range = candidateRange(position, radius);
        if (range.empty_)
            return;
        const auto visitCell = [&](const occupiedCell& cell) {
            for (int index = cell.begin_; index < cell.end_; ++index)
                visitor(entries_[index].pointIndex_);
        };
        if (range.scanOccupied_)
        {
            for (const auto& cell : cells_)
                if (range.all_ || (cell.cell_.x_ >= range.lower_.x_ && cell.cell_.x_ <= range.upper_.x_ &&
                                   cell.cell_.y_ >= range.lower_.y_ && cell.cell_.y_ <= range.upper_.y_ &&
                                   cell.cell_.z_ >= range.lower_.z_ && cell.cell_.z_ <= range.upper_.z_))
                    visitCell(cell);
            return;
        }
        // Explicit end checks avoid incrementing an integer already at its limit.
        for (auto z = range.lower_.z_;; ++z)
        {
            for (auto y = range.lower_.y_;; ++y)
            {
                for (auto x = range.lower_.x_;; ++x)
                {
                    if (const auto* cell = findCell({x, y, z}))
                        visitCell(*cell);
                    if (x == range.upper_.x_)
                        break;
                }
                if (y == range.upper_.y_)
                    break;
            }
            if (z == range.upper_.z_)
                break;
        }
    }

private:
    struct gridCell {
        std::int64_t x_{0}, y_{0}, z_{0}; ///< Checked cell coordinates.
        bool operator==(const gridCell& other) const noexcept { return x_ == other.x_ && y_ == other.y_ && z_ == other.z_; }
        bool operator<(const gridCell& other) const noexcept
        {
            return z_ != other.z_ ? z_ < other.z_ : (y_ != other.y_ ? y_ < other.y_ : x_ < other.x_);
        }
    };
    struct gridPoint {
        gridCell cell_;   ///< Occupied cell containing this point.
        int pointIndex_; ///< Original local input index.
    };
    struct occupiedCell {
        gridCell cell_; ///< Cell coordinate, ordered z/y/x.
        int begin_;     ///< Inclusive sorted-entry index.
        int end_;       ///< Exclusive sorted-entry index.
    };
    struct queryRange {
        gridCell lower_, upper_; ///< Inclusive occupied-coordinate bounds.
        bool empty_{true};       ///< No occupied cells intersect this query.
        bool scanOccupied_{false}; ///< Scan occupied rather than empty cells.
        bool all_{false};        ///< Coordinate limits require conservative fallback.
    };

    bool calculateCell(gridCell& cell, const Vec3& position) const noexcept;
    const occupiedCell* findCell(const gridCell& cell) const noexcept;
    queryRange candidateRange(const Vec3& position, Real radius) const;

    Vec3 origin_{Vec3::zero()}; ///< Minimum point coordinate used for quantization.
    gridCell maximumCell_{};   ///< Largest occupied coordinate on each axis.
    Real inverseCellSize_{1.0}; ///< Reciprocal cell size, where representable.
    bool unboundedCells_{false}; ///< Build requires conservative all-point queries.
    std::vector<gridPoint> entries_; ///< Owned original-index entries sorted by cell.
    std::vector<occupiedCell> cells_; ///< Compact occupied-cell ranges.
};

/** Finds exact self pairs with i < j and std::hypot(delta) <= radius.
 * Radius must be finite and nonnegative. Zero selects coincident distinct points.
 * Results are ordered by (firstIndex_, secondIndex_) regardless of thread count.
 * The CPU-only query does not reorder inputs or require contact/particle data.
 */
std::vector<pointPair> findPointPairs(const std::vector<math::Vec3>& points, math::Real radius);

/** Cross-array fixed-radius counterpart; each returned index is local to its array.
 * All qualifying cross pairs are retained, including equal local index numbers.
 * Invalid inputs and index limits have the same contract as pointGrid::build().
 */
std::vector<pointPair> findPointPairs(const std::vector<math::Vec3>& first, const std::vector<math::Vec3>& second, math::Real radius);

} // namespace fundem::cpu
