/**
 * @file periodicContactForce.h
 * @brief Additive CPU contact force with standalone X/Y or caller-selected periodic axes.
 */
#pragma once

#include "particle/particle.h"

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace fundem::validation
{

using math::Real;
using math::Vec3;

/** Normal force release convention, independent of the core's ordinary DEM laws. */
enum class normalReleasePolicy
{
    clampNonAttractive,     ///< Default: suppress an attractive dashpot resultant.
    retainUntilSeparation ///< Retain the signed resultant while overlap is positive; selected by the snapshot diagnostic.
};

/** Direct contact-pair coefficients, independent of FunDEM's material mixing. */
struct contactForceParameters {
    Real normalStiffness_{2.0e5};     ///< Normal spring coefficient.
    Real slidingStiffness_{5.7143e4}; ///< Tangential spring coefficient.
    Real normalDamping_{50.0};        ///< Normal force/relative-velocity dashpot coefficient.
    Real slidingDamping_{14.2857};    ///< Tangential force/relative-velocity dashpot coefficient.
    Real friction_{0.3};              ///< Coulomb sliding coefficient.
    normalReleasePolicy normalRelease_{normalReleasePolicy::clampNonAttractive}; ///< Existing validation behavior remains the default.
};

/** Current contact with the force and spring expressed on its first endpoint. */
struct shearContact {
    int first_{-1};                       ///< Finite endpoint; ascending index for two finite grains.
    int second_{-1};                      ///< Other zero-based container index, never a source-file ID.
    Vec3 point_{Vec3::zero()};            ///< Common contact point, wrapped into the enabled periodic cell.
    Vec3 normal_{Vec3::zero()};           ///< Unit normal from the second image center toward the first.
    Vec3 tangentialForce_{Vec3::zero()};  ///< Tangential spring-plus-dashpot force on the first endpoint.
    Vec3 tangentialSpring_{Vec3::zero()}; ///< Accumulated relative tangential displacement.
    Real overlap_{0.0};                   ///< Positive overlap of the two bounding spheres.
    Real normalForce_{0.0};               ///< Signed normal force; nonnegative under the default release policy.
};

/** Optional external boundary service; empty callbacks preserve standalone validation behavior. */
struct periodicBoundaryFunctions {
    const void* context_{nullptr};
    Vec3 (*wrapPosition_)(const void*, Vec3) noexcept{nullptr};
    Vec3 (*minimumImage_)(const void*, Vec3) noexcept{nullptr};
    bool (*isPeriodic_)(const void*, int) noexcept{nullptr}; ///< Null retains the standalone X/Y-periodic convention.
};

/**
 * Standalone additive spring/dashpot force with periodic search and pair history.
 *
 * The caller owns clearing forces and integrating motion. apply() only adds
 * forces and torques, so it can compose with other external loads. Particle
 * count and index ordering must stay fixed; particles have the common positive
 * radius supplied at construction. Enabled periodic centers are wrapped by the caller before
 * evaluation. Infinite-mass pairs are skipped and a finite endpoint is first.
 * No native FunDEM contact model is invoked by this component.
 */
class periodicContactForce
{
public:
    /** Validates the fixed population, domain and direct coefficients once. */
    periodicContactForce(const Vec3& minimum, const Vec3& maximum, Real radius, int particleCount, contactForceParameters parameters = {}, periodicBoundaryFunctions boundary = {});

    /** Adds one complete contact pass; zero time step refreshes without spring advance. */
    void apply(particleContainer::host_container_type& particles, Real historyTimeStep);

    const std::vector<shearContact>& contacts() const noexcept { return contacts_; }
    const contactForceParameters& parameters() const noexcept { return parameters_; }
    const Vec3& minimumBoundary() const noexcept { return domainMinimum_; }
    Vec3 maximumBoundary() const noexcept { return domainMinimum_ + domainLength_; }

    /** Maps enabled axes to the half-open cell; the standalone default enables X/Y. */
    Vec3 wrapPosition(Vec3 position) const noexcept;
    /** Returns nearest images only along enabled axes. */
    Vec3 minimumImage(Vec3 displacement) const noexcept;

private:
    /** One stable endpoint-pair history, independent of periodic image changes. */
    struct contactHistory {
        Vec3 spring_{Vec3::zero()};   ///< Tangential displacement on the deterministic first endpoint.
        std::uint64_t generation_{0}; ///< Last contact-pass generation in which the pair overlapped.
    };

    /** Rebuilds the periodic linked cells without reordering particle storage. */
    void buildCells(const particleContainer::host_container_type& particles);
    /** Accumulates one overlapping pair with minimum-image contact arms. */
    void evaluatePair(particleContainer::host_container_type& particles, int first, int second, Real historyTimeStep);
    /** Safely clamps a real cell coordinate before converting it to a container index. */
    int cellIndex(const Vec3& position) const noexcept;
    /** Produces an order-independent stable key from two container indices. */
    static std::uint64_t pairKey(int first, int second) noexcept;

    contactForceParameters parameters_;                           ///< Explicit spring, dashpot and friction coefficients.
    Vec3 domainMinimum_;                                          ///< Immutable periodic-domain lower corner.
    Vec3 domainLength_;                                           ///< Immutable X/Y periods and Z cell-list span.
    periodicBoundaryFunctions boundary_;                          ///< Non-owning optional service; its provider must outlive this force component.
    std::array<bool, 3> periodicAxes_{};                           ///< Cached boundary selection; no virtual dispatch in the neighbor-cell loops.
    Vec3 cellWidth_;                                              ///< At least one diameter, except for a single short Z cell.
    int cellCountX_{1};                                           ///< Number of periodic X cells.
    int cellCountY_{1};                                           ///< Number of periodic Y cells.
    int cellCountZ_{1};                                           ///< Number of clamped, nonperiodic Z cells.
    std::vector<int> cellHeads_;                                  ///< First stable particle index in each linked cell.
    std::vector<int> next_;                                       ///< Next stable particle index in a cell.
    std::unordered_map<std::uint64_t, contactHistory> histories_; ///< Persistent contact springs.
    std::vector<shearContact> contacts_;                          ///< Complete current contact set, not display-sampled.
    std::uint64_t generation_{0};                                 ///< Contact-pass generation used to remove separated histories.
};

} // namespace fundem::validation
