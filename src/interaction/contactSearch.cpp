#include "contactSearch.h"

#include "execution/contactDetection.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <utility>

namespace fundem
{
namespace
{

using Real = math::Real;
using Vec3 = math::Vec3;

constexpr int parallelParticleThreshold = 1024;
constexpr int parallelInteractionThreshold = 256;
constexpr long long parallelPairTestThreshold = 4096;

template <class Value> void compactValues(std::vector<Value>& values, const std::vector<unsigned char>& keep)
{
    const int valueCount = static_cast<int>(values.size());
    std::vector<int> offsets(valueCount + 1, 0);
    for (int index = 0; index < valueCount; ++index)
    {
        offsets[index + 1] = keep[index] != 0 ? 1 : 0;
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    std::vector<Value> compacted(offsets.back());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (valueCount >= parallelInteractionThreshold)
#endif
    for (int index = 0; index < valueCount; ++index)
    {
        if (keep[index] != 0)
        {
            compacted[offsets[index]] = std::move(values[index]);
        }
    }
    values = std::move(compacted);
}

template <class Value> void mergeLocalValues(std::vector<Value>& destination, std::vector<std::vector<Value>>& localValues)
{
    const int localValueCount = static_cast<int>(localValues.size());
    std::vector<int> offsets(localValueCount + 1, 0);
    for (int index = 0; index < localValueCount; ++index)
    {
        offsets[index + 1] = static_cast<int>(localValues[index].size());
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    destination.resize(offsets.back());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (offsets.back() >= parallelInteractionThreshold)
#endif
    for (int index = 0; index < localValueCount; ++index)
    {
        std::move(localValues[index].begin(), localValues[index].end(), destination.begin() + offsets[index]);
    }
}

struct contactKey {
    int masterParticleIndex_{-1};
    int slaveParticleIndex_{-1};
    int particleSurfaceNodeIndex_{-1};

    bool operator==(const contactKey& other) const noexcept
    {
        return masterParticleIndex_ == other.masterParticleIndex_ && slaveParticleIndex_ == other.slaveParticleIndex_ && particleSurfaceNodeIndex_ == other.particleSurfaceNodeIndex_;
    }
};

struct contactKeyHash {
    std::size_t operator()(const contactKey& key) const noexcept
    {
        std::size_t seed = std::hash<int>{}(key.masterParticleIndex_);
        combine(seed, std::hash<int>{}(key.slaveParticleIndex_));
        combine(seed, std::hash<int>{}(key.particleSurfaceNodeIndex_));
        return seed;
    }

private:
    static void combine(std::size_t& seed, std::size_t value) noexcept { seed ^= value + std::size_t{0x9e3779b9U} + (seed << 6U) + (seed >> 2U); }
};

struct contactHistoryState {
    Vec3 slidingSpringDeformation_{Vec3::zero()};
    Vec3 rollingSpringDeformation_{Vec3::zero()};
    Vec3 torsionalSpringDeformation_{Vec3::zero()};
};

using contactHistoryMap = std::unordered_map<contactKey, contactHistoryState, contactKeyHash>;

contactKey makeContactKey(const contact& currentContact) noexcept { return {currentContact.masterParticleIndex(), currentContact.slaveParticleIndex(), currentContact.particleSurfaceNodeIndex()}; }

contactHistoryMap collectContactHistory(const contactContainer& contacts)
{
    contactHistoryMap history;
    history.reserve(contacts.hostSize());
    for (const contact& currentContact : contacts.host())
    {
        history.insert_or_assign(makeContactKey(currentContact),
                                 contactHistoryState{currentContact.slidingSpringDeformation(), currentContact.rollingSpringDeformation(), currentContact.torsionalSpringDeformation()});
    }
    return history;
}

void restoreContactHistory(contact& currentContact, const contactHistoryMap& history)
{
    const auto previous = history.find(makeContactKey(currentContact));
    if (previous == history.end())
    {
        return;
    }
    currentContact.setSlidingSpringDeformation(previous->second.slidingSpringDeformation_);
    currentContact.setRollingSpringDeformation(previous->second.rollingSpringDeformation_);
    currentContact.setTorsionalSpringDeformation(previous->second.torsionalSpringDeformation_);
}

void readPositionsAndRadii(std::vector<Vec3>& positions, std::vector<Real>& radii, const particleContainer& spheres)
{
    const auto& sphereHost = spheres.host();
    const int particleCount = static_cast<int>(sphereHost.size());
    positions.resize(particleCount);
    radii.resize(particleCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        positions[particleIndex] = sphereHost[particleIndex].position();
        radii[particleIndex] = sphereHost[particleIndex].radius();
    }
}

void readPositionsAndRadii(std::vector<Vec3>& positions, std::vector<Real>& radii, const LSParticleContainer& LSParticles)
{
    const auto& LSParticleHost = LSParticles.host();
    const int particleCount = static_cast<int>(LSParticleHost.size());
    positions.resize(particleCount);
    radii.resize(particleCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        positions[particleIndex] = LSParticleHost[particleIndex].position();
        radii[particleIndex] = LSParticleHost[particleIndex].boundingRadius();
    }
}

void replaceContacts(contactContainer& destination, std::vector<contact>&& source)
{
    destination.resetDevice();
    destination.host() = std::move(source);
}

template <class ParticleStorage> void enforceInfiniteMassSlave(std::vector<particlePair>& pairs, const ParticleStorage& particles)
{
    const int pairCount = static_cast<int>(pairs.size());
    std::vector<unsigned char> keep(pairCount, 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (pairCount >= parallelInteractionThreshold)
#endif
    for (int pairIndex = 0; pairIndex < pairCount; ++pairIndex)
    {
        particlePair pair = pairs[pairIndex];
        const bool masterInfiniteMass = particles.host()[pair.masterParticleIndex_].inverseMass() == 0.0;
        const bool slaveInfiniteMass = particles.host()[pair.slaveParticleIndex_].inverseMass() == 0.0;
        if (masterInfiniteMass && slaveInfiniteMass)
        {
            continue;
        }
        if (masterInfiniteMass)
        {
            std::swap(pair.masterParticleIndex_, pair.slaveParticleIndex_);
        }
        pairs[pairIndex] = pair;
        keep[pairIndex] = 1;
    }
    compactValues(pairs, keep);
}

void enforceLevelSetSlaveOrdering(std::vector<particlePair>& pairs, const LSParticleContainer& particles)
{
    const int pairCount = static_cast<int>(pairs.size());
    std::vector<unsigned char> keep(pairCount, 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (pairCount >= parallelInteractionThreshold)
#endif
    for (int pairIndex = 0; pairIndex < pairCount; ++pairIndex)
    {
        particlePair pair = pairs[pairIndex];
        const LSParticle& master = particles.host()[pair.masterParticleIndex_];
        const LSParticle& slave = particles.host()[pair.slaveParticleIndex_];
        const bool masterInfiniteMass = master.inverseMass() == 0.0;
        const bool slaveInfiniteMass = slave.inverseMass() == 0.0;
        if (masterInfiniteMass && slaveInfiniteMass)
        {
            continue;
        }
        if (masterInfiniteMass || (!slaveInfiniteMass && master.boundingRadius() > slave.boundingRadius()))
        {
            std::swap(pair.masterParticleIndex_, pair.slaveParticleIndex_);
        }
        pairs[pairIndex] = pair;
        keep[pairIndex] = 1;
    }
    compactValues(pairs, keep);
}

contact makeSphereContact(bool& valid, const particleContainer& spheres, const particlePair& pair, const contactHistoryMap& history)
{
    contact newContact;
    Vec3 point;
    Vec3 normal;
    Real overlap = 0.0;
    Real area = 0.0;
    Real effectiveRadius = 0.0;
    const particle& masterSphere = spheres.host()[pair.masterParticleIndex_];
    const particle& slaveSphere = spheres.host()[pair.slaveParticleIndex_];
    valid = execution::detectSphereContact(point, normal, overlap, area, effectiveRadius, masterSphere.position(), masterSphere.radius(), slaveSphere.position(), slaveSphere.radius());
    if (!valid)
    {
        return newContact;
    }

    if (!newContact.setMasterSlaveParticle(spheres, pair.masterParticleIndex_, pair.slaveParticleIndex_))
    {
        valid = false;
        return newContact;
    }

    newContact.setPoint(point);
    newContact.setNormal(normal);
    newContact.setOverlap(overlap);
    newContact.setArea(area);
    newContact.setEffectiveMass(execution::effectiveMass(masterSphere.inverseMass(), slaveSphere.inverseMass()));
    newContact.setEffectiveRadius(effectiveRadius);
    restoreContactHistory(newContact, history);
    return newContact;
}

contact makeSphereLevelSetContact(bool& valid, const particleContainer& masterSpheres, const LSParticleContainer& slaveLSParticles, const particlePair& pair, const contactHistoryMap& history)
{
    contact newContact;
    const particle& masterSphere = masterSpheres.host()[pair.masterParticleIndex_];
    const LSParticle& slaveLSParticle = slaveLSParticles.host()[pair.slaveParticleIndex_];
    const std::vector<LSGridNode>& gridNodes = slaveLSParticle.gridNodes();
    const Real gridNodeInverseSpacing = slaveLSParticle.gridNodeInverseSpacing();
    if (gridNodes.empty() || gridNodeInverseSpacing <= 0.0)
    {
        valid = false;
        return newContact;
    }

    Real overlap = 0.0;
    Vec3 normal;
    valid = execution::detectLevelSetContact(overlap,
                                             normal,
                                             gridNodes.data(),
                                             masterSphere.position(),
                                             slaveLSParticle.position(),
                                             slaveLSParticle.orientation(),
                                             slaveLSParticle.gridNodeOrigin(),
                                             gridNodeInverseSpacing,
                                             slaveLSParticle.gridNodeSize(),
                                             0,
                                             masterSphere.radius());
    if (!valid || !newContact.setMasterSlaveParticle(masterSpheres, pair.masterParticleIndex_, slaveLSParticles, pair.slaveParticleIndex_))
    {
        valid = false;
        return newContact;
    }

    newContact.setPoint(execution::sphereLevelSetContactPoint(masterSphere.position(), masterSphere.radius(), normal, overlap));
    newContact.setNormal(normal);
    newContact.setOverlap(overlap);
    newContact.setArea(math::pi * masterSphere.radius() * overlap);
    newContact.setEffectiveMass(execution::effectiveMass(masterSphere.inverseMass(), slaveLSParticle.inverseMass()));
    // Physical radius used by the unified spherical contact model.
    newContact.setEffectiveRadius(masterSphere.radius());
    restoreContactHistory(newContact, history);
    return newContact;
}

void appendLevelSetContacts(std::vector<contact>& contacts, const LSParticleContainer& LSParticles, const particlePair& pair, const contactHistoryMap& history)
{
    const LSParticle& masterLSParticle = LSParticles.host()[pair.masterParticleIndex_];
    const LSParticle& slaveLSParticle = LSParticles.host()[pair.slaveParticleIndex_];
    const std::vector<LSGridNode>& gridNodes = slaveLSParticle.gridNodes();
    const Real gridNodeInverseSpacing = slaveLSParticle.gridNodeInverseSpacing();
    if (gridNodes.empty() || gridNodeInverseSpacing <= 0.0)
    {
        return;
    }

    const Real slaveRadiusSquared = slaveLSParticle.boundingRadius() * slaveLSParticle.boundingRadius();
    const std::vector<LSSurfaceNode>& surfaceNodes = masterLSParticle.surfaceNodes();
    const int contactBegin = static_cast<int>(contacts.size());
    for (int surfaceNodeIndex = 0; surfaceNodeIndex < static_cast<int>(surfaceNodes.size()); ++surfaceNodeIndex)
    {
        const LSSurfaceNode& surfaceNode = surfaceNodes[surfaceNodeIndex];
        const Vec3 surfaceNodePosition = masterLSParticle.position() + math::rotateUnit(masterLSParticle.orientation(), surfaceNode.position_);
        if (math::distanceSquared(surfaceNodePosition, slaveLSParticle.position()) > slaveRadiusSquared)
        {
            continue;
        }

        Real overlap = 0.0;
        Vec3 normal;
        if (!execution::detectLevelSetContact(overlap,
                                              normal,
                                              gridNodes.data(),
                                              surfaceNodePosition,
                                              slaveLSParticle.position(),
                                              slaveLSParticle.orientation(),
                                              slaveLSParticle.gridNodeOrigin(),
                                              gridNodeInverseSpacing,
                                              slaveLSParticle.gridNodeSize(),
                                              0))
        {
            continue;
        }

        contact newContact;
        if (!newContact.setMasterSlaveParticle(LSParticles, pair.masterParticleIndex_, pair.slaveParticleIndex_))
        {
            continue;
        }
        newContact.setParticleSurfaceNodeIndex(surfaceNodeIndex);
        newContact.setPoint(execution::surfaceNodeLevelSetContactPoint(surfaceNodePosition, normal, overlap));
        newContact.setNormal(normal);
        newContact.setOverlap(overlap);
        newContact.setArea(surfaceNode.area_);
        newContact.setEffectiveRadius(0.0);
        restoreContactHistory(newContact, history);
        contacts.push_back(std::move(newContact));
    }

    Real totalContactArea = 0.0;
    for (int contactIndex = contactBegin; contactIndex < static_cast<int>(contacts.size()); ++contactIndex)
    {
        totalContactArea += contacts[contactIndex].area();
    }

    const Real pairEffectiveMass = execution::effectiveMass(masterLSParticle.inverseMass(), slaveLSParticle.inverseMass());
    const Real effectiveMassPerArea = totalContactArea > 0.0 ? pairEffectiveMass / totalContactArea : 0.0;
    for (int contactIndex = contactBegin; contactIndex < static_cast<int>(contacts.size()); ++contactIndex)
    {
        contacts[contactIndex].setEffectiveMass(contacts[contactIndex].area() * effectiveMassPerArea);
    }
}

} // namespace

bool contactSearch::setCellSize(Real value) noexcept
{
    if (!math::isFinite(value) || value < 0.0)
    {
        return false;
    }
    requestedCellSize_ = value;
    return true;
}

void contactSearch::beginSearch() noexcept
{
    cells_.clear();
    gridParticles_.clear();
    candidatePairs_.clear();
}

void contactSearch::configureGrid(Real maximumRadius) noexcept
{
    effectiveCellSize_ = requestedCellSize_ > 0.0 ? requestedCellSize_ : 2.0 * maximumRadius;
    if (!math::isFinite(effectiveCellSize_) || effectiveCellSize_ <= math::defaultTolerance)
    {
        effectiveCellSize_ = 1.0;
    }
    inverseCellSize_ = 1.0 / effectiveCellSize_;
}

void contactSearch::buildGrid(const std::vector<Vec3>& positions)
{
    cells_.clear();
    gridParticles_.clear();
    if (positions.empty())
    {
        gridOrigin_ = Vec3::zero();
        return;
    }

    const int particleCount = static_cast<int>(positions.size());
    Real minimumX = positions.front().x;
    Real minimumY = positions.front().y;
    Real minimumZ = positions.front().z;
#if defined(_OPENMP)
#pragma omp parallel for reduction(min : minimumX, minimumY, minimumZ) schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        minimumX = std::min(minimumX, positions[particleIndex].x);
        minimumY = std::min(minimumY, positions[particleIndex].y);
        minimumZ = std::min(minimumZ, positions[particleIndex].z);
    }
    gridOrigin_ = {minimumX, minimumY, minimumZ};

    gridParticles_.resize(particleCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        gridParticles_[particleIndex] = {calculateCell(positions[particleIndex]), particleIndex};
    }
    std::sort(gridParticles_.begin(),
              gridParticles_.end(),
              [](const gridParticle& first, const gridParticle& second) { return first.cell_ == second.cell_ ? first.particleIndex_ < second.particleIndex_ : first.cell_ < second.cell_; });

    cells_.reserve(particleCount);
    int sortedParticleIndexBegin = 0;
    while (sortedParticleIndexBegin < particleCount)
    {
        int sortedParticleIndexEnd = sortedParticleIndexBegin + 1;
        while (sortedParticleIndexEnd < particleCount && gridParticles_[sortedParticleIndexEnd].cell_ == gridParticles_[sortedParticleIndexBegin].cell_)
        {
            ++sortedParticleIndexEnd;
        }
        cells_.push_back({gridParticles_[sortedParticleIndexBegin].cell_, sortedParticleIndexBegin, sortedParticleIndexEnd});
        sortedParticleIndexBegin = sortedParticleIndexEnd;
    }
}

contactSearch::gridCell contactSearch::calculateCell(const Vec3& position) const noexcept
{
    return {static_cast<int>(std::floor((position.x - gridOrigin_.x) * inverseCellSize_)),
            static_cast<int>(std::floor((position.y - gridOrigin_.y) * inverseCellSize_)),
            static_cast<int>(std::floor((position.z - gridOrigin_.z) * inverseCellSize_))};
}

const contactSearch::occupiedCell* contactSearch::findCell(const gridCell& cell) const noexcept
{
    const auto value = std::lower_bound(cells_.begin(), cells_.end(), cell, [](const occupiedCell& occupied, const gridCell& target) { return occupied.cell_ < target; });
    return value != cells_.end() && value->cell_ == cell ? &*value : nullptr;
}

int contactSearch::searchCellRange(Real distance) const noexcept { return static_cast<int>(std::ceil(distance * inverseCellSize_)); }

void contactSearch::findCandidatePairs(const std::vector<Vec3>& positions, const std::vector<Real>& radii)
{
    const int particleCount = static_cast<int>(positions.size());
    Real maximumRadius = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(max : maximumRadius) schedule(static) if (particleCount >= parallelParticleThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        maximumRadius = std::max(maximumRadius, radii[particleIndex]);
    }
    configureGrid(maximumRadius);

    std::vector<std::vector<particlePair>> particlePairs(particleCount);
    if (directSearch_)
    {
        [[maybe_unused]] const long long pairTestCount = static_cast<long long>(particleCount) * (particleCount - 1) / 2;
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (pairTestCount >= parallelPairTestThreshold)
#endif
        for (int masterIndex = 0; masterIndex < particleCount; ++masterIndex)
        {
            std::vector<particlePair>& localPairs = particlePairs[masterIndex];
            for (int slaveIndex = masterIndex + 1; slaveIndex < particleCount; ++slaveIndex)
            {
                if (execution::boundingSpheresOverlap(positions[masterIndex], radii[masterIndex], positions[slaveIndex], radii[slaveIndex]))
                {
                    localPairs.push_back({masterIndex, slaveIndex});
                }
            }
        }
    }
    else
    {
        buildGrid(positions);
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (particleCount >= parallelInteractionThreshold)
#endif
        for (int masterIndex = 0; masterIndex < particleCount; ++masterIndex)
        {
            std::vector<particlePair>& localPairs = particlePairs[masterIndex];
            const gridCell center = calculateCell(positions[masterIndex]);
            const int cellRange = searchCellRange(radii[masterIndex] + maximumRadius);
            for (int zOffset = -cellRange; zOffset <= cellRange; ++zOffset)
            {
                for (int yOffset = -cellRange; yOffset <= cellRange; ++yOffset)
                {
                    for (int xOffset = -cellRange; xOffset <= cellRange; ++xOffset)
                    {
                        const occupiedCell* cell = findCell({center.x_ + xOffset, center.y_ + yOffset, center.z_ + zOffset});
                        if (cell == nullptr)
                        {
                            continue;
                        }
                        for (int sortedParticleIndex = cell->sortedParticleIndexBegin_; sortedParticleIndex < cell->sortedParticleIndexEnd_; ++sortedParticleIndex)
                        {
                            const int slaveIndex = gridParticles_[sortedParticleIndex].particleIndex_;
                            if (slaveIndex <= masterIndex)
                            {
                                continue;
                            }
                            if (execution::boundingSpheresOverlap(positions[masterIndex], radii[masterIndex], positions[slaveIndex], radii[slaveIndex]))
                            {
                                localPairs.push_back({masterIndex, slaveIndex});
                            }
                        }
                    }
                }
            }
        }
    }
    mergeLocalValues(candidatePairs_, particlePairs);
}

void contactSearch::findSphereLevelSetCandidatePairs(const std::vector<Vec3>& masterPositions,
                                                     const std::vector<Real>& masterRadii,
                                                     const std::vector<Vec3>& slavePositions,
                                                     const std::vector<Real>& slaveRadii)
{
    const int masterCount = static_cast<int>(masterPositions.size());
    Real maximumRadius = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(max : maximumRadius) schedule(static) if (masterCount >= parallelParticleThreshold)
#endif
    for (int masterIndex = 0; masterIndex < masterCount; ++masterIndex)
    {
        maximumRadius = std::max(maximumRadius, masterRadii[masterIndex]);
    }
    const int slaveCount = static_cast<int>(slaveRadii.size());
#if defined(_OPENMP)
#pragma omp parallel for reduction(max : maximumRadius) schedule(static) if (slaveCount >= parallelParticleThreshold)
#endif
    for (int slaveIndex = 0; slaveIndex < slaveCount; ++slaveIndex)
    {
        maximumRadius = std::max(maximumRadius, slaveRadii[slaveIndex]);
    }
    configureGrid(maximumRadius);

    std::vector<std::vector<particlePair>> particlePairs(masterCount);
    if (directSearch_)
    {
        [[maybe_unused]] const long long pairTestCount = static_cast<long long>(masterCount) * slaveCount;
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (pairTestCount >= parallelPairTestThreshold)
#endif
        for (int masterIndex = 0; masterIndex < masterCount; ++masterIndex)
        {
            std::vector<particlePair>& localPairs = particlePairs[masterIndex];
            for (int slaveIndex = 0; slaveIndex < static_cast<int>(slavePositions.size()); ++slaveIndex)
            {
                if (execution::boundingSpheresOverlap(masterPositions[masterIndex], masterRadii[masterIndex], slavePositions[slaveIndex], slaveRadii[slaveIndex]))
                {
                    localPairs.push_back({masterIndex, slaveIndex});
                }
            }
        }
    }
    else
    {
        Real maximumSlaveRadius = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for reduction(max : maximumSlaveRadius) schedule(static) if (slaveCount >= parallelParticleThreshold)
#endif
        for (int slaveIndex = 0; slaveIndex < slaveCount; ++slaveIndex)
        {
            maximumSlaveRadius = std::max(maximumSlaveRadius, slaveRadii[slaveIndex]);
        }
        buildGrid(slavePositions);
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (masterCount >= parallelInteractionThreshold)
#endif
        for (int masterIndex = 0; masterIndex < masterCount; ++masterIndex)
        {
            std::vector<particlePair>& localPairs = particlePairs[masterIndex];
            const gridCell center = calculateCell(masterPositions[masterIndex]);
            const int cellRange = searchCellRange(masterRadii[masterIndex] + maximumSlaveRadius);
            for (int zOffset = -cellRange; zOffset <= cellRange; ++zOffset)
            {
                for (int yOffset = -cellRange; yOffset <= cellRange; ++yOffset)
                {
                    for (int xOffset = -cellRange; xOffset <= cellRange; ++xOffset)
                    {
                        const occupiedCell* cell = findCell({center.x_ + xOffset, center.y_ + yOffset, center.z_ + zOffset});
                        if (cell == nullptr)
                        {
                            continue;
                        }
                        for (int sortedParticleIndex = cell->sortedParticleIndexBegin_; sortedParticleIndex < cell->sortedParticleIndexEnd_; ++sortedParticleIndex)
                        {
                            const int slaveIndex = gridParticles_[sortedParticleIndex].particleIndex_;
                            if (execution::boundingSpheresOverlap(masterPositions[masterIndex], masterRadii[masterIndex], slavePositions[slaveIndex], slaveRadii[slaveIndex]))
                            {
                                localPairs.push_back({masterIndex, slaveIndex});
                            }
                        }
                    }
                }
            }
        }
    }
    mergeLocalValues(candidatePairs_, particlePairs);
}

int contactSearch::findContacts(contactContainer& contacts, const particleContainer& spheres)
{
    beginSearch();
    const contactHistoryMap history = collectContactHistory(contacts);
    std::vector<Vec3> positions;
    std::vector<Real> radii;
    readPositionsAndRadii(positions, radii, spheres);
    findCandidatePairs(positions, radii);
    enforceInfiniteMassSlave(candidatePairs_, spheres);

    const int pairCount = static_cast<int>(candidatePairs_.size());
    std::vector<contact> detectedContacts(pairCount);
    std::vector<unsigned char> valid(pairCount, 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (pairCount >= parallelInteractionThreshold)
#endif
    for (int pairIndex = 0; pairIndex < pairCount; ++pairIndex)
    {
        bool contactValid = false;
        detectedContacts[pairIndex] = makeSphereContact(contactValid, spheres, candidatePairs_[pairIndex], history);
        valid[pairIndex] = contactValid ? 1 : 0;
    }

    compactValues(detectedContacts, valid);
    const int contactCount = static_cast<int>(detectedContacts.size());
    replaceContacts(contacts, std::move(detectedContacts));
    return contactCount;
}

int contactSearch::findContacts(contactContainer& contacts, const particleContainer& masterSpheres, const LSParticleContainer& slaveLSParticles)
{
    beginSearch();
    const contactHistoryMap history = collectContactHistory(contacts);
    std::vector<Vec3> masterSpherePositions;
    std::vector<Real> masterSphereRadii;
    std::vector<Vec3> slaveLSParticlePositions;
    std::vector<Real> slaveLSParticleRadii;
    readPositionsAndRadii(masterSpherePositions, masterSphereRadii, masterSpheres);
    readPositionsAndRadii(slaveLSParticlePositions, slaveLSParticleRadii, slaveLSParticles);
    findSphereLevelSetCandidatePairs(masterSpherePositions, masterSphereRadii, slaveLSParticlePositions, slaveLSParticleRadii);

    const int pairCount = static_cast<int>(candidatePairs_.size());
    std::vector<contact> detectedContacts(pairCount);
    std::vector<unsigned char> valid(pairCount, 0);
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 8) if (pairCount >= parallelInteractionThreshold)
#endif
    for (int pairIndex = 0; pairIndex < pairCount; ++pairIndex)
    {
        bool contactValid = false;
        detectedContacts[pairIndex] = makeSphereLevelSetContact(contactValid, masterSpheres, slaveLSParticles, candidatePairs_[pairIndex], history);
        valid[pairIndex] = contactValid ? 1 : 0;
    }

    compactValues(detectedContacts, valid);
    const int contactCount = static_cast<int>(detectedContacts.size());
    replaceContacts(contacts, std::move(detectedContacts));
    return contactCount;
}

int contactSearch::findContacts(contactContainer& contacts, const LSParticleContainer& LSParticles)
{
    beginSearch();
    const contactHistoryMap history = collectContactHistory(contacts);
    std::vector<Vec3> positions;
    std::vector<Real> radii;
    readPositionsAndRadii(positions, radii, LSParticles);
    findCandidatePairs(positions, radii);
    enforceLevelSetSlaveOrdering(candidatePairs_, LSParticles);

    const int pairCount = static_cast<int>(candidatePairs_.size());
    std::vector<std::vector<contact>> pairContacts(pairCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 4) if (pairCount >= parallelInteractionThreshold)
#endif
    for (int pairIndex = 0; pairIndex < pairCount; ++pairIndex)
    {
        appendLevelSetContacts(pairContacts[pairIndex], LSParticles, candidatePairs_[pairIndex], history);
    }

    std::vector<contact> activeContacts;
    mergeLocalValues(activeContacts, pairContacts);
    const int contactCount = static_cast<int>(activeContacts.size());
    replaceContacts(contacts, std::move(activeContacts));
    return contactCount;
}

void contactSearch::clear() noexcept
{
    effectiveCellSize_ = 1.0;
    inverseCellSize_ = 1.0;
    gridOrigin_ = Vec3::zero();
    gridParticles_.clear();
    cells_.clear();
    candidatePairs_.clear();
}

} // namespace fundem
