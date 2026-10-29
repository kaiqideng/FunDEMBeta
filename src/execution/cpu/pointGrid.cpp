#include "pointGrid.h"

#include "parallelPolicy.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fundem::cpu
{
namespace
{
using Real = math::Real;
using Vec3 = math::Vec3;

void requireIndexable(std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("Point count exceeds the local index representation.");
}

void requireRadius(Real radius)
{
    if (!math::isFinite(radius) || radius < 0.0)
        throw std::invalid_argument("Point query radius must be finite and nonnegative.");
}

/** Conversion bounds are tested before casting; the upper Real bound rounds up. */
bool cellCoordinate(std::int64_t& result, Real value) noexcept
{
    if (!math::isFinite(value) || value < static_cast<Real>(std::numeric_limits<std::int64_t>::min()) ||
        value >= static_cast<Real>(std::numeric_limits<std::int64_t>::max()))
        return false;
    result = static_cast<std::int64_t>(std::floor(value));
    return true;
}

std::vector<pointPair> findPairs(const std::vector<Vec3>& first, const std::vector<Vec3>& second, Real radius, bool self)
{
    requireRadius(radius);
    requireIndexable(first.size());
    pointGrid grid;
    grid.build(second, radius > 0.0 ? radius : 1.0);
    if (!self)
        for (const auto& point : first)
            if (!math::isFinite(point))
                throw std::invalid_argument("Point positions must be finite.");
    if (first.empty() || second.empty())
        return {};

    const int count = static_cast<int>(first.size());
    std::vector<std::vector<pointPair>> localPairs(first.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (count >= parallelInteractionThreshold)
#endif
    for (int firstIndex = 0; firstIndex < count; ++firstIndex)
    {
        auto& pairs = localPairs[firstIndex];
        grid.forEachCandidate(first[firstIndex], radius, [&](int secondIndex) {
            if (self && secondIndex <= firstIndex)
                return;
            const auto difference = second[secondIndex] - first[firstIndex];
            if (std::hypot(difference.x, difference.y, difference.z) <= radius)
                pairs.push_back({firstIndex, secondIndex});
        });
        std::sort(pairs.begin(), pairs.end(), [](const auto& left, const auto& right) { return left.secondIndex_ < right.secondIndex_; });
    }
    std::size_t pairCount = 0;
    for (const auto& pairs : localPairs)
    {
        if (pairs.size() > std::numeric_limits<std::size_t>::max() - pairCount)
            throw std::length_error("Point-pair count exceeds the storage representation.");
        pairCount += pairs.size();
    }
    std::vector<pointPair> pairs;
    pairs.reserve(pairCount);
    for (auto& local : localPairs)
        pairs.insert(pairs.end(), local.begin(), local.end());
    return pairs;
}
} // namespace

void pointGrid::build(const std::vector<Vec3>& positions, Real cellSize)
{
    if (!math::isFinite(cellSize) || cellSize <= 0.0)
        throw std::invalid_argument("Point-grid cell size must be finite and positive.");
    requireIndexable(positions.size());
    clear();
    if (positions.empty())
        return;

    const int count = static_cast<int>(positions.size());
    Real minimumX = positions.front().x, minimumY = positions.front().y, minimumZ = positions.front().z;
    int invalid = 0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(min : minimumX, minimumY, minimumZ) reduction(| : invalid) schedule(static) if (count >= parallelParticleThreshold)
#endif
    for (int index = 0; index < count; ++index)
    {
        const auto& point = positions[index];
        invalid |= !math::isFinite(point);
        minimumX = std::min(minimumX, point.x);
        minimumY = std::min(minimumY, point.y);
        minimumZ = std::min(minimumZ, point.z);
    }
    if (invalid)
        throw std::invalid_argument("Point positions must be finite.");
    origin_ = {minimumX, minimumY, minimumZ};
    inverseCellSize_ = 1.0 / cellSize;
    entries_.resize(positions.size());
    int overflow = 0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(| : overflow) schedule(static) if (count >= parallelParticleThreshold)
#endif
    for (int index = 0; index < count; ++index)
    {
        gridCell cell;
        overflow |= !calculateCell(cell, positions[index]);
        entries_[index] = {cell, index};
    }
    unboundedCells_ = overflow != 0;
    if (unboundedCells_)
        for (auto& entry : entries_)
            entry.cell_ = {};
    std::sort(entries_.begin(), entries_.end(), [](const auto& left, const auto& right) {
        return left.cell_ == right.cell_ ? left.pointIndex_ < right.pointIndex_ : left.cell_ < right.cell_;
    });
    cells_.reserve(positions.size());
    for (int begin = 0; begin < count;)
    {
        int end = begin + 1;
        while (end < count && entries_[end].cell_ == entries_[begin].cell_)
            ++end;
        const auto cell = entries_[begin].cell_;
        cells_.push_back({cell, begin, end});
        maximumCell_.x_ = std::max(maximumCell_.x_, cell.x_);
        maximumCell_.y_ = std::max(maximumCell_.y_, cell.y_);
        maximumCell_.z_ = std::max(maximumCell_.z_, cell.z_);
        begin = end;
    }
}

bool pointGrid::calculateCell(gridCell& cell, const Vec3& position) const noexcept
{
    return cellCoordinate(cell.x_, (position.x - origin_.x) * inverseCellSize_) &&
           cellCoordinate(cell.y_, (position.y - origin_.y) * inverseCellSize_) &&
           cellCoordinate(cell.z_, (position.z - origin_.z) * inverseCellSize_);
}

const pointGrid::occupiedCell* pointGrid::findCell(const gridCell& cell) const noexcept
{
    const auto found = std::lower_bound(cells_.begin(), cells_.end(), cell, [](const occupiedCell& entry, const gridCell& value) { return entry.cell_ < value; });
    return found != cells_.end() && found->cell_ == cell ? &*found : nullptr;
}

pointGrid::queryRange pointGrid::candidateRange(const Vec3& position, Real radius) const
{
    requireRadius(radius);
    if (!math::isFinite(position))
        throw std::invalid_argument("Point query position must be finite.");
    queryRange range;
    if (cells_.empty())
        return range;
    range.empty_ = false;
    gridCell center;
    std::int64_t distance;
    if (unboundedCells_ || !calculateCell(center, position) || !cellCoordinate(distance, std::ceil(radius * inverseCellSize_)))
    {
        range.scanOccupied_ = range.all_ = true;
        return range;
    }
    const auto minimum = std::numeric_limits<std::int64_t>::min();
    const auto maximum = std::numeric_limits<std::int64_t>::max();
    const auto outsideLimits = [&](std::int64_t coordinate) { return coordinate < minimum + distance || coordinate > maximum - distance; };
    if (outsideLimits(center.x_) || outsideLimits(center.y_) || outsideLimits(center.z_))
    {
        range.scanOccupied_ = range.all_ = true;
        return range;
    }
    // A large origin offset can round two nearby world points into cells more
    // than ceil(radius/cellSize) apart. Map outward-rounded world bounds through
    // the same monotone quantization as build(), then union the original range.
    // Expanding radius too covers subtraction cancellation in the exact test.
    const Real infinity = std::numeric_limits<Real>::infinity();
    const Real expandedRadius = std::nextafter(radius, infinity);
    const Vec3 worldLower{std::nextafter(position.x - expandedRadius, -infinity), std::nextafter(position.y - expandedRadius, -infinity),
                          std::nextafter(position.z - expandedRadius, -infinity)};
    const Vec3 worldUpper{std::nextafter(position.x + expandedRadius, infinity), std::nextafter(position.y + expandedRadius, infinity),
                          std::nextafter(position.z + expandedRadius, infinity)};
    gridCell outwardLower, outwardUpper;
    if (!calculateCell(outwardLower, worldLower) || !calculateCell(outwardUpper, worldUpper))
    {
        range.scanOccupied_ = range.all_ = true;
        return range;
    }
    range.lower_ = {std::max(std::int64_t{0}, std::min(center.x_ - distance, outwardLower.x_)),
                     std::max(std::int64_t{0}, std::min(center.y_ - distance, outwardLower.y_)),
                     std::max(std::int64_t{0}, std::min(center.z_ - distance, outwardLower.z_))};
    range.upper_ = {std::min(maximumCell_.x_, std::max(center.x_ + distance, outwardUpper.x_)),
                     std::min(maximumCell_.y_, std::max(center.y_ + distance, outwardUpper.y_)),
                     std::min(maximumCell_.z_, std::max(center.z_ + distance, outwardUpper.z_))};
    range.empty_ = range.lower_.x_ > range.upper_.x_ || range.lower_.y_ > range.upper_.y_ || range.lower_.z_ > range.upper_.z_;
    if (!range.empty_)
    {
        const long double x = static_cast<long double>(range.upper_.x_) - range.lower_.x_ + 1.0L;
        const long double y = static_cast<long double>(range.upper_.y_) - range.lower_.y_ + 1.0L;
        const long double z = static_cast<long double>(range.upper_.z_) - range.lower_.z_ + 1.0L;
        range.scanOccupied_ = x * y * z > static_cast<long double>(cells_.size());
    }
    return range;
}

void pointGrid::clear() noexcept
{
    origin_ = Vec3::zero();
    maximumCell_ = {};
    inverseCellSize_ = 1.0;
    unboundedCells_ = false;
    entries_.clear();
    cells_.clear();
}

std::vector<pointPair> findPointPairs(const std::vector<Vec3>& points, Real radius) { return findPairs(points, points, radius, true); }

std::vector<pointPair> findPointPairs(const std::vector<Vec3>& first, const std::vector<Vec3>& second, Real radius)
{
    return findPairs(first, second, radius, false);
}

} // namespace fundem::cpu
