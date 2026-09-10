#include "sphInteraction.h"

#include "execution/sphFunctions.h"
#include "math/Indexing.h"
#include "solverFunctions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace fundem::cpu
{
namespace
{

using math::Real;
using math::Vec3;

struct gridEntry {
    int hash_{-1};
    int particleIndex_{-1};
};

class uniformSpatialGrid
{
public:
    template <class PositionGetter>
    uniformSpatialGrid(int particleCount, PositionGetter position, const Vec3& minimumBoundary, const Vec3& maximumBoundary, Real cellSize)
        : minimumBoundary_(minimumBoundary), inverseCellSize_(1.0 / cellSize)
    {
        const Vec3 extent = (maximumBoundary - minimumBoundary) * inverseCellSize_;
        size_ = make_int3(std::max(1, static_cast<int>(std::ceil(extent.x))), std::max(1, static_cast<int>(std::ceil(extent.y))), std::max(1, static_cast<int>(std::ceil(extent.z))));

        std::vector<gridEntry> entries(static_cast<std::size_t>(particleCount));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
        for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
        {
            entries[particleIndex] = {hash(position(particleIndex)), particleIndex};
        }
        std::sort(entries.begin(),
                  entries.end(),
                  [](const gridEntry& first, const gridEntry& second) { return first.hash_ < second.hash_ || (first.hash_ == second.hash_ && first.particleIndex_ < second.particleIndex_); });

        particleIndices_.resize(entries.size());
        for (int sortedIndex = 0; sortedIndex < particleCount; ++sortedIndex)
        {
            const gridEntry& entry = entries[sortedIndex];
            particleIndices_[sortedIndex] = entry.particleIndex_;
            if (sortedIndex == 0 || entry.hash_ != entries[sortedIndex - 1].hash_)
            {
                cellHashes_.push_back(entry.hash_);
                cellOffsets_.push_back(sortedIndex);
            }
        }
        cellOffsets_.push_back(particleCount);
    }

    template <class Function> void forEachCandidate(const Vec3& position, Real searchRadius, Function function) const
    {
        const Vec3 local = (position - minimumBoundary_) * inverseCellSize_;
        const Real searchExtent = searchRadius * inverseCellSize_;
        const int3 begin = make_int3(clampCoordinate(static_cast<int>(std::floor(local.x - searchExtent)), size_.x),
                                     clampCoordinate(static_cast<int>(std::floor(local.y - searchExtent)), size_.y),
                                     clampCoordinate(static_cast<int>(std::floor(local.z - searchExtent)), size_.z));
        const int3 end = make_int3(clampCoordinate(static_cast<int>(std::floor(local.x + searchExtent)), size_.x),
                                   clampCoordinate(static_cast<int>(std::floor(local.y + searchExtent)), size_.y),
                                   clampCoordinate(static_cast<int>(std::floor(local.z + searchExtent)), size_.z));
        for (int z = begin.z; z <= end.z; ++z)
        {
            for (int y = begin.y; y <= end.y; ++y)
            {
                for (int x = begin.x; x <= end.x; ++x)
                {
                    const int cellHash = math::linearIndex(x, y, z, size_.x, size_.y);
                    const auto location = std::lower_bound(cellHashes_.begin(), cellHashes_.end(), cellHash);
                    if (location == cellHashes_.end() || *location != cellHash)
                    {
                        continue;
                    }
                    const int cellIndex = static_cast<int>(location - cellHashes_.begin());
                    for (int sortedIndex = cellOffsets_[cellIndex]; sortedIndex < cellOffsets_[cellIndex + 1]; ++sortedIndex)
                    {
                        function(particleIndices_[sortedIndex]);
                    }
                }
            }
        }
    }

private:
    static int clampCoordinate(int coordinate, int size) noexcept { return coordinate < 0 ? 0 : (coordinate < size ? coordinate : size - 1); }

    int hash(const Vec3& position) const noexcept
    {
        const Vec3 local = (position - minimumBoundary_) * inverseCellSize_;
        const int x = clampCoordinate(static_cast<int>(std::floor(local.x)), size_.x);
        const int y = clampCoordinate(static_cast<int>(std::floor(local.y)), size_.y);
        const int z = clampCoordinate(static_cast<int>(std::floor(local.z)), size_.z);
        return math::linearIndex(x, y, z, size_.x, size_.y);
    }

    Vec3 minimumBoundary_{Vec3::zero()};
    Real inverseCellSize_{0.0};
    int3 size_{1, 1, 1};
    std::vector<int> cellHashes_;
    std::vector<int> cellOffsets_;
    std::vector<int> particleIndices_;
};

template <class PositionGetter>
void buildNeighborList(std::vector<int>& offsets, std::vector<int>& indices, int particleCount, PositionGetter position, const uniformSpatialGrid& neighborGrid, Real searchRadius)
{
    std::vector<int> counts(static_cast<std::size_t>(particleCount), 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        neighborGrid.forEachCandidate(position(particleIndex), searchRadius, [&](int) { ++counts[particleIndex]; });
    }

    offsets.assign(static_cast<std::size_t>(particleCount) + 1, 0);
    std::partial_sum(counts.begin(), counts.end(), offsets.begin() + 1);
    indices.resize(static_cast<std::size_t>(offsets.back()));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        int writeIndex = offsets[particleIndex];
        neighborGrid.forEachCandidate(position(particleIndex), searchRadius, [&](int neighborIndex) { indices[writeIndex++] = neighborIndex; });
    }
}

void addAtomic(Vec3& target, const Vec3& value) noexcept
{
#if defined(_OPENMP)
#pragma omp atomic update
#endif
    target.x += value.x;
#if defined(_OPENMP)
#pragma omp atomic update
#endif
    target.y += value.y;
#if defined(_OPENMP)
#pragma omp atomic update
#endif
    target.z += value.z;
}

} // namespace

void SPHInteraction::initializeParticles(SPHParticleContainer& particles, Real referenceDensity, Real soundSpeed) const
{
    auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        particle.setForce(Vec3::zero());
        particle.setPriorForce(Vec3::zero());
        particle.setDensityRate(0.0);
        particle.setPressure(execution::computePressureFromDensity(particle.density(), referenceDensity, soundSpeed));
    }
}

void SPHInteraction::buildNeighborhood(const SPHParticleContainer& particles,
                                       const virtualParticleContainer& virtualParticles,
                                       const Vec3& minimumBoundary,
                                       const Vec3& maximumBoundary,
                                       Real cellSize,
                                       Real searchRadius)
{
    const auto& particleHost = particles.host();
    const auto& virtualParticleHost = virtualParticles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    const int virtualParticleCount = static_cast<int>(virtualParticleHost.size());
    const uniformSpatialGrid particleGrid(particleCount, [&](int index) -> const Vec3& { return particleHost[index].position(); }, minimumBoundary, maximumBoundary, cellSize);
    const uniformSpatialGrid virtualParticleGrid(virtualParticleCount, [&](int index) -> const Vec3& { return virtualParticleHost[index].position(); }, minimumBoundary, maximumBoundary, cellSize);

    buildNeighborList(particleNeighbors_.offsets_, particleNeighbors_.indices_, particleCount, [&](int index) -> const Vec3& { return particleHost[index].position(); }, particleGrid, searchRadius);
    buildNeighborList(
        virtualParticleNeighbors_.offsets_,
        virtualParticleNeighbors_.indices_,
        particleCount,
        [&](int index) -> const Vec3& { return particleHost[index].position(); },
        virtualParticleGrid,
        searchRadius);
}

void SPHInteraction::updateFreeSurface(SPHParticleContainer& particles, const virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const
{
    auto& particleHost = particles.host();
    const auto& virtualParticleHost = virtualParticles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        Real positionDivergence = 0.0;
        for (const int* neighbor = particleNeighbors_.begin(particleIndex); neighbor != particleNeighbors_.end(particleIndex); ++neighbor)
        {
            if (*neighbor == particleIndex)
            {
                continue;
            }
            const SPHParticle& neighborParticle = particleHost[*neighbor];
            const Real neighborVolume = parameters.particleMass_ / execution::nonZeroDensity(neighborParticle.density());
            positionDivergence += execution::freeSurfacePositionDivergenceContribution(particle.position() - neighborParticle.position(), neighborVolume, parameters.smoothingLength_);
        }
        for (const int* neighbor = virtualParticleNeighbors_.begin(particleIndex); neighbor != virtualParticleNeighbors_.end(particleIndex); ++neighbor)
        {
            const virtualParticle& wall = virtualParticleHost[*neighbor];
            positionDivergence += execution::freeSurfacePositionDivergenceContribution(particle.position() - wall.position(), wall.volume(), parameters.smoothingLength_);
        }
        particle.setFreeSurface(execution::isFreeSurface(positionDivergence));
    }
}

void SPHInteraction::reinitializeDensity(SPHParticleContainer& particles, const virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const
{
    auto& particleHost = particles.host();
    const auto& virtualParticleHost = virtualParticles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    const Real supportRadius = 2.0 * parameters.smoothingLength_;
    const Real referenceVolume = parameters.particleMass_ / parameters.referenceDensity_;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        if (particle.isConstrained())
        {
            continue;
        }
        Real kernelSum = execution::wendlandKernel3D(0.0, parameters.smoothingLength_);
        for (const int* neighbor = particleNeighbors_.begin(particleIndex); neighbor != particleNeighbors_.end(particleIndex); ++neighbor)
        {
            if (*neighbor == particleIndex)
            {
                continue;
            }
            const Real distance = math::norm(particle.position() - particleHost[*neighbor].position());
            if (distance < supportRadius)
            {
                kernelSum += execution::wendlandKernel3D(distance, parameters.smoothingLength_);
            }
        }
        for (const int* neighbor = virtualParticleNeighbors_.begin(particleIndex); neighbor != virtualParticleNeighbors_.end(particleIndex); ++neighbor)
        {
            const virtualParticle& wall = virtualParticleHost[*neighbor];
            kernelSum += wall.volume() / referenceVolume * execution::wendlandKernel3D(math::norm(particle.position() - wall.position()), parameters.smoothingLength_);
        }
        const Real density = execution::reinitializedDensity(kernelSum, parameters.latticeKernelSum3D_, parameters.referenceDensity_);
        particle.setDensity(density);
        particle.setPressure(execution::computePressureFromDensity(density, parameters.referenceDensity_, parameters.soundSpeed_));
    }
}

void SPHInteraction::updateDensity(SPHParticleContainer& particles, const virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters, Real timeStep) const
{
    auto& particleHost = particles.host();
    const auto& virtualParticleHost = virtualParticles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        if (particle.isConstrained())
        {
            if (math::isFinite(particle.densityRate()))
                particle.setDensityRate(0.0);
            continue;
        }
        Real densityRate = 0.0;
        for (const int* neighbor = particleNeighbors_.begin(particleIndex); neighbor != particleNeighbors_.end(particleIndex); ++neighbor)
        {
            if (*neighbor == particleIndex)
            {
                continue;
            }
            const SPHParticle& neighborParticle = particleHost[*neighbor];
            densityRate += execution::calculateFluidDensityRate(particle.position(),
                                                                particle.velocity(),
                                                                particle.density(),
                                                                particle.pressure(),
                                                                neighborParticle.position(),
                                                                neighborParticle.velocity(),
                                                                neighborParticle.density(),
                                                                neighborParticle.pressure(),
                                                                parameters.particleMass_,
                                                                parameters.referenceDensity_,
                                                                parameters.smoothingLength_,
                                                                parameters.soundSpeed_);
        }
        for (const int* neighbor = virtualParticleNeighbors_.begin(particleIndex); neighbor != virtualParticleNeighbors_.end(particleIndex); ++neighbor)
        {
            const virtualParticle& wall = virtualParticleHost[*neighbor];
            densityRate += execution::calculateWallDensityRate(particle.position(),
                                                               particle.velocity(),
                                                               particle.density(),
                                                               particle.pressure(),
                                                               wall.position(),
                                                               wall.normal(),
                                                               wall.velocity(),
                                                               wall.acceleration(),
                                                               wall.volume(),
                                                               parameters.referenceDensity_,
                                                               parameters.smoothingLength_,
                                                               parameters.soundSpeed_,
                                                               parameters.gravity_);
        }
        // Keep an invalid computed rate latched until initialization so a later split stage cannot hide it.
        if (math::isFinite(particle.densityRate()))
            particle.*SPHParticle::densityRateField::member = densityRate;
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        if (particle.isConstrained())
        {
            continue;
        }
        Real density = particle.density();
        Real pressure = particle.pressure();
        execution::integrateDensityAndPressure(density, pressure, particle.densityRate(), parameters.referenceDensity_, parameters.soundSpeed_, timeStep);
        particle.setDensity(density);
        particle.setPressure(pressure);
    }
}

void SPHInteraction::updatePriorForceAndBoundaryForce(SPHParticleContainer& particles, virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const
{
    auto& particleHost = particles.host();
    auto& virtualParticleHost = virtualParticles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    std::vector<Vec3> wallForces(virtualParticleHost.size(), Vec3::zero());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        Vec3 acceleration = Vec3::zero();
        for (const int* neighbor = particleNeighbors_.begin(particleIndex); neighbor != particleNeighbors_.end(particleIndex); ++neighbor)
        {
            if (*neighbor == particleIndex)
            {
                continue;
            }
            const SPHParticle& neighborParticle = particleHost[*neighbor];
            acceleration += execution::calculateFluidViscousAcceleration(particle.position(),
                                                                         particle.velocity(),
                                                                         particle.density(),
                                                                         neighborParticle.position(),
                                                                         neighborParticle.velocity(),
                                                                         neighborParticle.density(),
                                                                         parameters.particleMass_,
                                                                         parameters.dynamicViscosity_,
                                                                         parameters.smoothingLength_);
        }
        for (const int* neighbor = virtualParticleNeighbors_.begin(particleIndex); neighbor != virtualParticleNeighbors_.end(particleIndex); ++neighbor)
        {
            const virtualParticle& wall = virtualParticleHost[*neighbor];
            const Vec3 wallAcceleration = execution::calculateWallViscousAcceleration(particle.position(),
                                                                                      particle.velocity(),
                                                                                      particle.density(),
                                                                                      wall.position(),
                                                                                      wall.velocity(),
                                                                                      wall.volume(),
                                                                                      parameters.dynamicViscosity_,
                                                                                      parameters.smoothingLength_);
            acceleration += wallAcceleration;
            addAtomic(wallForces[*neighbor], execution::calculateWallForce(parameters.particleMass_, wallAcceleration));
        }
        particle.setPriorForce(parameters.particleMass_ * acceleration);
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (static_cast<int>(virtualParticleHost.size()) >= parallelParticleThreshold)
#endif
    for (int wallIndex = 0; wallIndex < static_cast<int>(virtualParticleHost.size()); ++wallIndex)
    {
        virtualParticleHost[wallIndex].setPriorForce(wallForces[wallIndex]);
    }
}

void SPHInteraction::updatePressureForceAndBoundaryForce(SPHParticleContainer& particles, virtualParticleContainer& virtualParticles, const SPHInteractionParameters& parameters) const
{
    auto& particleHost = particles.host();
    auto& virtualParticleHost = virtualParticles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    std::vector<Vec3> wallForces(virtualParticleHost.size());
    for (int wallIndex = 0; wallIndex < static_cast<int>(virtualParticleHost.size()); ++wallIndex)
    {
        wallForces[wallIndex] = virtualParticleHost[wallIndex].priorForce();
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        SPHParticle& particle = particleHost[particleIndex];
        Vec3 acceleration = Vec3::zero();
        for (const int* neighbor = particleNeighbors_.begin(particleIndex); neighbor != particleNeighbors_.end(particleIndex); ++neighbor)
        {
            if (*neighbor == particleIndex)
            {
                continue;
            }
            const SPHParticle& neighborParticle = particleHost[*neighbor];
            acceleration += execution::calculateFluidPressureAcceleration(particle.position(),
                                                                          particle.velocity(),
                                                                          particle.density(),
                                                                          particle.pressure(),
                                                                          neighborParticle.position(),
                                                                          neighborParticle.velocity(),
                                                                          neighborParticle.density(),
                                                                          neighborParticle.pressure(),
                                                                          parameters.particleMass_,
                                                                          parameters.referenceDensity_,
                                                                          parameters.smoothingLength_,
                                                                          parameters.soundSpeed_);
        }
        for (const int* neighbor = virtualParticleNeighbors_.begin(particleIndex); neighbor != virtualParticleNeighbors_.end(particleIndex); ++neighbor)
        {
            const virtualParticle& wall = virtualParticleHost[*neighbor];
            const Vec3 wallAcceleration = execution::calculateWallPressureAcceleration(particle.position(),
                                                                                       particle.velocity(),
                                                                                       particle.density(),
                                                                                       particle.pressure(),
                                                                                       wall.position(),
                                                                                       wall.normal(),
                                                                                       wall.velocity(),
                                                                                       wall.acceleration(),
                                                                                       wall.volume(),
                                                                                       parameters.referenceDensity_,
                                                                                       parameters.smoothingLength_,
                                                                                       parameters.soundSpeed_,
                                                                                       parameters.gravity_);
            acceleration += wallAcceleration;
            addAtomic(wallForces[*neighbor], execution::calculateWallForce(parameters.particleMass_, wallAcceleration));
        }
        particle.setForce(particle.priorForce() + parameters.particleMass_ * acceleration);
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (static_cast<int>(virtualParticleHost.size()) >= parallelParticleThreshold)
#endif
    for (int wallIndex = 0; wallIndex < static_cast<int>(virtualParticleHost.size()); ++wallIndex)
    {
        virtualParticleHost[wallIndex].setForce(wallForces[wallIndex]);
    }
}

void SPHInteraction::integrateVelocity(SPHParticleContainer& particles, const Vec3& gravity, Real timeStep) const
{
    auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        particleHost[particleIndex].updateVelocity(gravity, timeStep);
    }
}

void SPHInteraction::integratePosition(SPHParticleContainer& particles, Real timeStep) const
{
    auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        particleHost[particleIndex].updatePosition(timeStep);
    }
}

SPHStateStatistics SPHInteraction::stateStatistics(const SPHParticleContainer& particles, const Vec3& gravity) const
{
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    if (particleCount == 0)
    {
        return {};
    }
    Real maximumVelocity = 0.0;
    Real maximumAcceleration = 0.0;
    Real minimumDensity = std::numeric_limits<Real>::max();
    Real maximumDensity = 0.0;
    int invalidValueCount = 0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold) reduction(max : maximumVelocity, maximumAcceleration, maximumDensity) reduction(min : minimumDensity)        \
    reduction(+ : invalidValueCount)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const SPHParticle& particle = particleHost[particleIndex];
        const Vec3 acceleration = particle.inverseMass() * particle.force() + gravity;
        const bool validVelocity = math::isFinite(particle.velocity());
        const bool validAcceleration = math::isFinite(acceleration);
        const bool validDensity = math::isFinite(particle.density());
        maximumVelocity = std::max(maximumVelocity, validVelocity ? math::norm(particle.velocity()) : 0.0);
        maximumAcceleration = std::max(maximumAcceleration, validAcceleration ? math::norm(acceleration) : 0.0);
        minimumDensity = std::min(minimumDensity, validDensity ? particle.density() : 0.0);
        maximumDensity = std::max(maximumDensity, validDensity ? particle.density() : 0.0);
        const bool validThermodynamics = math::isFinite(particle.pressure()) && math::isFinite(particle.densityRate());
        invalidValueCount += validVelocity && validAcceleration && validDensity && validThermodynamics && math::isFinite(particle.position()) ? 0 : 1;
    }
    return {maximumVelocity, maximumAcceleration, minimumDensity, maximumDensity, invalidValueCount};
}

SPHKinematicsStatistics SPHInteraction::kinematicsStatistics(const virtualParticleContainer& particles) const
{
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    Real maximumVelocity = 0.0;
    Real maximumAcceleration = 0.0;
    int invalidValueCount = 0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold) reduction(max : maximumVelocity, maximumAcceleration) reduction(+ : invalidValueCount)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const virtualParticle& particle = particleHost[particleIndex];
        const bool validVelocity = math::isFinite(particle.velocity());
        const bool validAcceleration = math::isFinite(particle.acceleration());
        maximumVelocity = std::max(maximumVelocity, validVelocity ? math::norm(particle.velocity()) : 0.0);
        maximumAcceleration = std::max(maximumAcceleration, validAcceleration ? math::norm(particle.acceleration()) : 0.0);
        invalidValueCount += validVelocity && validAcceleration ? 0 : 1;
    }
    return {maximumVelocity, maximumAcceleration, invalidValueCount};
}

} // namespace fundem::cpu
