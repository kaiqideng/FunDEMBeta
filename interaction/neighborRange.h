/**
 * @file neighborRange.h
 * @brief Defines neighbor ranges, candidate entries, and rebuildable histories.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "math/Constants.h"

namespace fundem
{

/** Count and inclusive prefix sum associated with one search owner. */
class neighborRange
{
private:
    int count_{0};     ///< Number of neighbors or contacts owned by this entry.
    int prefixSum_{0}; ///< Inclusive total through this entry.

public:
    /** Mutable device arrays used for counting and prefix-sum construction. */
    struct device_type {
        int* count_{nullptr};     ///< Per-owner counts.
        int* prefixSum_{nullptr}; ///< Inclusive prefix sums.
    };

    struct countField {
        inline static constexpr auto member = &neighborRange::count_;
    };
    struct prefixSumField {
        inline static constexpr auto member = &neighborRange::prefixSum_;
    };

    using DeviceLayout = deviceLayout<countField, prefixSumField>;
};

using neighborRangeContainer = hostAoSDeviceSoA<neighborRange, neighborRange::DeviceLayout>;

/** Previous-step inclusive prefix sum used to locate contact history ranges. */
class neighborRangeHistory
{
public:
    neighborRangeHistory() = default;
    /** Stores one previous-step inclusive prefix sum. */
    explicit neighborRangeHistory(int prefixSum) noexcept : prefixSum_(prefixSum) {}

private:
    int prefixSum_{0}; ///< Previous-step inclusive prefix sum.

public:
    /** Mutable device view of previous-step prefix sums. */
    struct device_type {
        int* prefixSum_{nullptr}; ///< Previous-step inclusive prefix sums.
    };

    struct prefixSumField {
        inline static constexpr auto member = &neighborRangeHistory::prefixSum_;
    };

    using DeviceLayout = deviceLayout<prefixSumField>;
};

using neighborRangeHistoryContainer = hostAoSDeviceSoA<neighborRangeHistory, neighborRangeHistory::DeviceLayout>;

/** One candidate neighbor entry and its LS contact-area contribution. */
class particleNeighbor
{
private:
    int particleIndex_{-1};     ///< Neighbor particle-container index.
    math::Real patchArea_{0.0}; ///< Surface patch area assigned to this pair.

public:
    /** Mutable device view of candidate-neighbor entries. */
    struct device_type {
        int* particleIndex_{nullptr};    ///< Neighbor indices.
        math::Real* patchArea_{nullptr}; ///< Associated patch areas.
    };

    struct particleIndexField {
        inline static constexpr auto member = &particleNeighbor::particleIndex_;
    };
    struct patchAreaField {
        inline static constexpr auto member = &particleNeighbor::patchArea_;
    };

    using DeviceLayout = deviceLayout<particleIndexField, patchAreaField>;
};

using particleNeighborContainer = hostAoSDeviceSoA<particleNeighbor, particleNeighbor::DeviceLayout>;

} // namespace fundem
