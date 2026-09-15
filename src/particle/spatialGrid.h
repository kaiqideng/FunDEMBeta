/**
 * @file spatialGrid.h
 * @brief Defines single- and multi-level particle background-grid storage.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "math/Vector3.h"
#include "particle.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace fundem
{

inline constexpr int maximumSpatialGridLevelCount = 4;

/** Hash and original particle index for one grid-sorted entry. */
class spatialGridParticle
{
public:
    spatialGridParticle() = default;
    explicit spatialGridParticle(int index) noexcept : Index_(index) {}

    /** Mutable device arrays used while hashing and sorting particles. */
    struct device_type {
        int* particleHash_{nullptr};  ///< Hash key for each sorted entry.
        int* particleIndex_{nullptr}; ///< Original particle-container index.
    };

    /** Read-only view consumed by neighbor traversal. */
    struct const_device_type {
        const int* particleHash_{nullptr};  ///< Sorted hash keys.
        const int* particleIndex_{nullptr}; ///< Original particle-container indices.
    };

    struct hashField;
    struct indexField;

private:
    int hash_{-1};  ///< Cell hash assigned before sorting.
    int Index_{-1}; ///< Original particle-container index.

public:
    struct hashField {
        inline static constexpr auto member = &spatialGridParticle::hash_;
    };
    struct indexField {
        inline static constexpr auto member = &spatialGridParticle::Index_;
    };

    using DeviceLayout = deviceLayout<hashField, indexField>;
};

using spatialGridParticleContainer = hostAoSDeviceSoA<spatialGridParticle, spatialGridParticle::DeviceLayout>;

/** One occupied-cell range in the sorted spatial-grid particle array. */
class spatialGridCell
{
public:
    /** Mutable occupied-cell device arrays. */
    struct device_type {
        int* cellHash_{nullptr};                     ///< Hash of each occupied cell.
        int* cellSortedParticleIndexBegin_{nullptr}; ///< Inclusive begin in the sorted particle array.
        int* cellSortedParticleIndexEnd_{nullptr};   ///< Exclusive end in the sorted particle array.
    };

    /** Read-only occupied-cell device arrays. */
    struct const_device_type {
        const int* cellHash_{nullptr};                     ///< Hash of each occupied cell.
        const int* cellSortedParticleIndexBegin_{nullptr}; ///< Inclusive sorted-entry begin.
        const int* cellSortedParticleIndexEnd_{nullptr};   ///< Exclusive sorted-entry end.
    };

    struct hashField;
    struct sortedParticleIndexBeginField;
    struct sortedParticleIndexEndField;

private:
    int hash_{-1};                     ///< Occupied-cell hash.
    int sortedParticleIndexBegin_{-1}; ///< Inclusive sorted-entry begin.
    int sortedParticleIndexEnd_{-1};   ///< Exclusive sorted-entry end.

public:
    struct hashField {
        inline static constexpr auto member = &spatialGridCell::hash_;
    };
    struct sortedParticleIndexBeginField {
        inline static constexpr auto member = &spatialGridCell::sortedParticleIndexBegin_;
    };
    struct sortedParticleIndexEndField {
        inline static constexpr auto member = &spatialGridCell::sortedParticleIndexEnd_;
    };

    using DeviceLayout = deviceLayout<hashField, sortedParticleIndexBeginField, sortedParticleIndexEndField>;
};

using spatialGridCellContainer = hostAoSDeviceSoA<spatialGridCell, spatialGridCell::DeviceLayout>;

/** Geometry and packed storage ranges of one spatial-grid level. */
struct spatialGridLevel {
    math::Real maximumDiameter_{0.0};                ///< Largest particle diameter assigned to this level.
    math::Vec3 inverseCellSize_{math::Vec3::zero()}; ///< Reciprocal cell dimensions.
    int3 size_{1, 1, 1};                             ///< Cell counts along each axis.
    int particleBegin_{0};                           ///< Inclusive level range in sorted particle storage.
    int particleEnd_{0};                             ///< Exclusive level particle range.
    int occupiedCellBegin_{0};                       ///< Inclusive level range in occupied-cell storage.
    int occupiedCellEnd_{0};                         ///< Exclusive occupied-cell range.
    bool hasInfiniteMass_{false};                    ///< Whether this is the final fixed-body level.
};

/** Owns up to four size-banded background grids and their device search arrays. */
class spatialGridContainer
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    /** Mutable complete device view used by grid construction kernels. */
    struct device_type {
        int* particleHash_{nullptr};                              ///< Hash per level-ordered particle entry.
        int* particleIndex_{nullptr};                             ///< Original particle index per sorted entry.
        int* cellHash_{nullptr};                                  ///< Hash per occupied or dense cell.
        int* cellSortedParticleIndexBegin_{nullptr};              ///< Inclusive sorted range begins.
        int* cellSortedParticleIndexEnd_{nullptr};                ///< Exclusive sorted range ends.
        Vec3 minimumBoundary_{Vec3::zero()};                      ///< World origin used by every grid level.
        spatialGridLevel levels_[maximumSpatialGridLevelCount]{}; ///< Per-level dimensions and offsets.
        int levelCount_{1};                                       ///< Number of active grid levels.
    };

    /** Read-only complete device view used by search kernels. */
    struct const_device_type {
        const int* particleHash_{nullptr};                        ///< Hash per level-ordered particle entry.
        const int* particleIndex_{nullptr};                       ///< Original particle index per sorted entry.
        const int* cellHash_{nullptr};                            ///< Hash per occupied or dense cell.
        const int* cellSortedParticleIndexBegin_{nullptr};        ///< Inclusive sorted range begins.
        const int* cellSortedParticleIndexEnd_{nullptr};          ///< Exclusive sorted range ends.
        Vec3 minimumBoundary_{Vec3::zero()};                      ///< World origin used by every grid level.
        spatialGridLevel levels_[maximumSpatialGridLevelCount]{}; ///< Per-level dimensions and offsets.
        int levelCount_{1};                                       ///< Number of active grid levels.
    };

    spatialGridContainer() = default;
    spatialGridContainer(const spatialGridContainer&) = delete;
    spatialGridContainer& operator=(const spatialGridContainer&) = delete;
    spatialGridContainer(spatialGridContainer&&) noexcept = default;
    spatialGridContainer& operator=(spatialGridContainer&&) noexcept = default;

    const Vec3& minimumBoundary() const noexcept { return minimumBoundary_; }
    const Vec3& maximumBoundary() const noexcept { return maximumBoundary_; }
    int levelCount() const noexcept { return levelCount_; }
    const spatialGridLevel& level(int index) const noexcept { return levels_[index]; }
    int particleCount() const noexcept { return static_cast<int>(particles_.deviceSize()); }
    int particleCount(int levelIndex) const noexcept { return levels_[levelIndex].particleEnd_ - levels_[levelIndex].particleBegin_; }
    int cellCount(int levelIndex) const noexcept { return levels_[levelIndex].occupiedCellEnd_ - levels_[levelIndex].occupiedCellBegin_; }
    int cellCount() const noexcept
    {
        int count = 0;
        for (int levelIndex = 0; levelIndex < levelCount_; ++levelIndex)
        {
            count += cellCount(levelIndex);
        }
        return count;
    }
    bool usesDenseCells() const noexcept { return denseCells_; }
    double deviceMemoryGB() const noexcept { return particles_.deviceMemoryGB() + cells_.deviceMemoryGB(); }

    /**
     * Configures adaptive diameter bands, appends an infinite-mass final level
     * when needed, assigns particles, and uploads the unsorted entries.
     * @tparam ParticleStorage Host/device particle container type.
     * @param particles Source particle container.
     * @param minimumBoundary Lower search-domain bound.
     * @param maximumBoundary Upper search-domain bound; degenerate axes expand.
     * @param stream CUDA stream receiving the initial upload.
     */
    template <class ParticleStorage> void set(const ParticleStorage& particles, const Vec3& minimumBoundary, Vec3 maximumBoundary, cudaStream_t stream = nullptr)
    {
        if (!math::isFinite(minimumBoundary) || !math::isFinite(maximumBoundary))
        {
            throw std::invalid_argument("Invalid spatial-grid boundary.");
        }

        denseCells_ = false;
        Real minimumFiniteMassDiameter = 0.0;
        Real maximumFiniteMassDiameter = 0.0;
        Real maximumInfiniteMassDiameter = 0.0;
        bool hasInfiniteMassParticles = false;
        for (const auto& value : particles.host())
        {
            const particle& body = static_cast<const particle&>(value);
            const Real diameter = 2.0 * body.radius();
            if (body.inverseMass() == 0.0)
            {
                hasInfiniteMassParticles = true;
                if (math::isFinite(diameter) && diameter > maximumInfiniteMassDiameter)
                {
                    maximumInfiniteMassDiameter = diameter;
                }
                continue;
            }
            if (!math::isFinite(diameter) || diameter <= math::defaultTolerance)
            {
                continue;
            }
            minimumFiniteMassDiameter = minimumFiniteMassDiameter > 0.0 && minimumFiniteMassDiameter < diameter ? minimumFiniteMassDiameter : diameter;
            maximumFiniteMassDiameter = maximumFiniteMassDiameter > diameter ? maximumFiniteMassDiameter : diameter;
        }

        setLevels(minimumFiniteMassDiameter, maximumFiniteMassDiameter, maximumInfiniteMassDiameter, hasInfiniteMassParticles);
        setBoundary(minimumBoundary, maximumBoundary, maximumFiniteMassDiameter > maximumInfiniteMassDiameter ? maximumFiniteMassDiameter : maximumInfiniteMassDiameter);
        assignParticles(particles);
        setLevelGridValues();

        cells_.host().assign(particles.hostSize(), spatialGridCell{});
        cells_.resizeDevice(particles.hostSize());
        particles_.copyHostToDeviceAsync(stream);
    }

    /**
     * Configures one dense uniform grid, primarily for SPH.
     * @tparam ParticleStorage Host/device particle container type.
     * @param particles Source particle container.
     * @param cellSize Positive uniform cell size.
     * @param minimumBoundary Lower domain bound.
     * @param maximumBoundary Upper domain bound.
     * @param stream CUDA stream receiving the particle-index upload.
     */
    template <class ParticleStorage> void setUniform(const ParticleStorage& particles, Real cellSize, const Vec3& minimumBoundary, Vec3 maximumBoundary, cudaStream_t stream = nullptr)
    {
        if (!math::isFinite(cellSize) || cellSize <= math::defaultTolerance || !math::isFinite(minimumBoundary) || !math::isFinite(maximumBoundary))
        {
            throw std::invalid_argument("Invalid uniform spatial-grid configuration.");
        }

        denseCells_ = true;
        levels_ = {};
        levelCount_ = 1;
        levels_[0].maximumDiameter_ = cellSize;
        setBoundary(minimumBoundary, maximumBoundary, cellSize);

        particles_.host().clear();
        particles_.host().reserve(particles.hostSize());
        levels_[0].particleBegin_ = 0;
        levels_[0].occupiedCellBegin_ = 0;
        for (int particleIndex = 0; particleIndex < static_cast<int>(particles.hostSize()); ++particleIndex)
        {
            particles_.host().emplace_back(particleIndex);
        }
        levels_[0].particleEnd_ = static_cast<int>(particles_.hostSize());
        setLevelGridValues();

        const int3 size = levels_[0].size_;
        levels_[0].occupiedCellEnd_ = size.x * size.y * size.z;
        cells_.host().clear();
        cells_.resizeDevice(static_cast<std::size_t>(levels_[0].occupiedCellEnd_));
        particles_.copyHostToDeviceAsync(stream);
    }

    spatialGridParticle::device_type particles() noexcept { return {particles_.device<spatialGridParticle::hashField>(), particles_.device<spatialGridParticle::indexField>()}; }

    spatialGridParticle::const_device_type particles() const noexcept { return {particles_.device<spatialGridParticle::hashField>(), particles_.device<spatialGridParticle::indexField>()}; }

    spatialGridCell::device_type cells() noexcept
    {
        return {cells_.device<spatialGridCell::hashField>(), cells_.device<spatialGridCell::sortedParticleIndexBeginField>(), cells_.device<spatialGridCell::sortedParticleIndexEndField>()};
    }

    spatialGridCell::const_device_type cells() const noexcept
    {
        return {cells_.device<spatialGridCell::hashField>(), cells_.device<spatialGridCell::sortedParticleIndexBeginField>(), cells_.device<spatialGridCell::sortedParticleIndexEndField>()};
    }

    /** Builds a mutable device view without allocating or copying. */
    device_type mutableDevice() noexcept
    {
        device_type result;
        result.particleHash_ = particles_.device<spatialGridParticle::hashField>();
        result.particleIndex_ = particles_.device<spatialGridParticle::indexField>();
        result.cellHash_ = cells_.device<spatialGridCell::hashField>();
        result.cellSortedParticleIndexBegin_ = cells_.device<spatialGridCell::sortedParticleIndexBeginField>();
        result.cellSortedParticleIndexEnd_ = cells_.device<spatialGridCell::sortedParticleIndexEndField>();
        result.minimumBoundary_ = minimumBoundary_;
        result.levelCount_ = levelCount_;
        for (int levelIndex = 0; levelIndex < levelCount_; ++levelIndex)
        {
            result.levels_[levelIndex] = levels_[levelIndex];
        }
        return result;
    }

    /** Builds a read-only device view without allocating or copying. */
    const_device_type device() const noexcept
    {
        const_device_type result;
        result.particleHash_ = particles_.device<spatialGridParticle::hashField>();
        result.particleIndex_ = particles_.device<spatialGridParticle::indexField>();
        result.cellHash_ = cells_.device<spatialGridCell::hashField>();
        result.cellSortedParticleIndexBegin_ = cells_.device<spatialGridCell::sortedParticleIndexBeginField>();
        result.cellSortedParticleIndexEnd_ = cells_.device<spatialGridCell::sortedParticleIndexEndField>();
        result.minimumBoundary_ = minimumBoundary_;
        result.levelCount_ = levelCount_;
        for (int levelIndex = 0; levelIndex < levelCount_; ++levelIndex)
        {
            result.levels_[levelIndex] = levels_[levelIndex];
        }
        return result;
    }

    void setOccupiedCellCount(int levelIndex, int count) noexcept { levels_[levelIndex].occupiedCellEnd_ = levels_[levelIndex].occupiedCellBegin_ + count; }

private:
    /** Chooses finite-diameter bands and an optional final infinite-mass level. */
    void setLevels(Real minimumFiniteMassDiameter, Real maximumFiniteMassDiameter, Real maximumInfiniteMassDiameter, bool hasInfiniteMassParticles) noexcept
    {
        levels_ = {};
        int finiteMassLevelCount = 0;
        if (minimumFiniteMassDiameter > math::defaultTolerance && maximumFiniteMassDiameter > math::defaultTolerance)
        {
            const Real diameterRatio = maximumFiniteMassDiameter / minimumFiniteMassDiameter;
            const int requestedLevelCount = diameterRatio <= 4.0 ? 1 : (diameterRatio <= 16.0 ? 2 : 3);
            const int availableLevelCount = maximumSpatialGridLevelCount - (hasInfiniteMassParticles ? 1 : 0);
            finiteMassLevelCount = requestedLevelCount < availableLevelCount ? requestedLevelCount : availableLevelCount;
            const Real levelRatio = std::pow(diameterRatio, 1.0 / Real(finiteMassLevelCount));
            for (int levelIndex = 0; levelIndex < finiteMassLevelCount; ++levelIndex)
            {
                levels_[levelIndex].maximumDiameter_ = levelIndex + 1 == finiteMassLevelCount ? maximumFiniteMassDiameter : minimumFiniteMassDiameter * std::pow(levelRatio, Real(levelIndex + 1));
            }
        }

        levelCount_ = finiteMassLevelCount;
        if (hasInfiniteMassParticles)
        {
            spatialGridLevel& infiniteMassLevel = levels_[levelCount_++];
            infiniteMassLevel.maximumDiameter_ = maximumInfiniteMassDiameter > math::defaultTolerance ? maximumInfiniteMassDiameter : 1.0;
            infiniteMassLevel.hasInfiniteMass_ = true;
        }
        if (levelCount_ == 0)
        {
            levelCount_ = 1;
        }
    }

    /** Stores a non-degenerate domain, expanding collapsed axes when required. */
    void setBoundary(const Vec3& minimumBoundary, Vec3 maximumBoundary, Real maximumDiameter) noexcept
    {
        const Real fallbackSize = maximumDiameter > math::defaultTolerance ? maximumDiameter : 1.0;
        if (maximumBoundary.x <= minimumBoundary.x)
        {
            maximumBoundary.x = minimumBoundary.x + fallbackSize;
        }
        if (maximumBoundary.y <= minimumBoundary.y)
        {
            maximumBoundary.y = minimumBoundary.y + fallbackSize;
        }
        if (maximumBoundary.z <= minimumBoundary.z)
        {
            maximumBoundary.z = minimumBoundary.z + fallbackSize;
        }
        minimumBoundary_ = minimumBoundary;
        maximumBoundary_ = maximumBoundary;
    }

    /** Assigns stable particle indices to size levels without reordering users' objects. */
    template <class ParticleStorage> void assignParticles(const ParticleStorage& particles)
    {
        std::array<std::vector<int>, maximumSpatialGridLevelCount> particleIndices;
        const int infiniteMassLevelIndex = levels_[levelCount_ - 1].hasInfiniteMass_ ? levelCount_ - 1 : -1;
        const int finiteMassLevelCount = infiniteMassLevelIndex >= 0 ? infiniteMassLevelIndex : levelCount_;
        for (int particleIndex = 0; particleIndex < static_cast<int>(particles.hostSize()); ++particleIndex)
        {
            const particle& body = static_cast<const particle&>(particles.host()[particleIndex]);
            const Real diameter = 2.0 * body.radius();
            if (body.inverseMass() == 0.0 && infiniteMassLevelIndex >= 0)
            {
                particleIndices[infiniteMassLevelIndex].push_back(particleIndex);
                continue;
            }
            int levelIndex = 0;
            while (levelIndex + 1 < finiteMassLevelCount && diameter > levels_[levelIndex].maximumDiameter_)
            {
                ++levelIndex;
            }
            particleIndices[levelIndex].push_back(particleIndex);
        }

        particles_.host().clear();
        particles_.host().reserve(particles.hostSize());
        for (int levelIndex = 0; levelIndex < levelCount_; ++levelIndex)
        {
            levels_[levelIndex].particleBegin_ = static_cast<int>(particles_.hostSize());
            levels_[levelIndex].occupiedCellBegin_ = levels_[levelIndex].particleBegin_;
            for (const int particleIndex : particleIndices[levelIndex])
            {
                particles_.host().emplace_back(particleIndex);
            }
            levels_[levelIndex].particleEnd_ = static_cast<int>(particles_.hostSize());
            levels_[levelIndex].occupiedCellEnd_ = levels_[levelIndex].occupiedCellBegin_;
        }
    }

    /** Derives cell size, reciprocal size, and per-level grid dimensions. */
    void setLevelGridValues() noexcept
    {
        const Vec3 domainSize = maximumBoundary_ - minimumBoundary_;
        for (int levelIndex = 0; levelIndex < levelCount_; ++levelIndex)
        {
            spatialGridLevel& level = levels_[levelIndex];
            if (level.maximumDiameter_ <= math::defaultTolerance)
            {
                level.inverseCellSize_ = Vec3::zero();
                level.size_ = make_int3(1, 1, 1);
                continue;
            }
            level.size_.x = domainSize.x > level.maximumDiameter_ ? static_cast<int>(domainSize.x / level.maximumDiameter_) : 1;
            level.size_.y = domainSize.y > level.maximumDiameter_ ? static_cast<int>(domainSize.y / level.maximumDiameter_) : 1;
            level.size_.z = domainSize.z > level.maximumDiameter_ ? static_cast<int>(domainSize.z / level.maximumDiameter_) : 1;
            level.inverseCellSize_ = {Real(level.size_.x) / domainSize.x, Real(level.size_.y) / domainSize.y, Real(level.size_.z) / domainSize.z};
        }
    }

    Vec3 minimumBoundary_{Vec3::zero()};                                  ///< Lower search-domain bound.
    Vec3 maximumBoundary_{Vec3::zero()};                                  ///< Upper search-domain bound.
    std::array<spatialGridLevel, maximumSpatialGridLevelCount> levels_{}; ///< Active level descriptors.
    int levelCount_{1};                                                   ///< Number of active levels.
    bool denseCells_{false};                                              ///< Whether every uniform cell has storage.
    spatialGridParticleContainer particles_;                              ///< Packed grid-entry arrays.
    spatialGridCellContainer cells_;                                      ///< Packed occupied-cell arrays.
};

} // namespace fundem
