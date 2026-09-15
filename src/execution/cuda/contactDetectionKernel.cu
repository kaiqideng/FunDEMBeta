#include "contactDetectionKernel.cuh"

#include "atomicKernel.cuh"
#include "execution/contactDetection.h"

namespace fundem::cuda
{
namespace
{

constexpr int blockSize = 256;

using detail::contactHistoryDeviceView;
using detail::contactOutputDeviceView;
using detail::levelSetGeometryDeviceView;
using detail::levelSetParticleDeviceView;
using detail::particleDeviceView;
using detail::surfaceNodeMappingDeviceView;

template <class ParticleStorage> particleDeviceView makeParticleDeviceView(const ParticleStorage& particles) noexcept
{
    return {particles.template device<rigidBody::positionField>(),
            particles.template device<rigidBody::orientationField>(),
            particles.template device<rigidBody::inverseMassField>(),
            particles.template device<particle::radiusField>()};
}

levelSetParticleDeviceView makeLevelSetParticleDeviceView(const LSParticleContainer& particles) noexcept
{
    return {makeParticleDeviceView(particles), particles.device<LSParticle::geometryIndexField>()};
}

levelSetGeometryDeviceView makeLevelSetGeometryDeviceView(const LSGeometryContainer& geometries) noexcept { return {geometries.descriptors(), geometries.gridNodes(), geometries.surfaceNodes()}; }

surfaceNodeMappingDeviceView makeSurfaceNodeMappingDeviceView(const surfaceNodeMappingContainer& mappings) noexcept
{
    return {mappings.device<surfaceNodeMapping::geometrySurfaceNodeIndexField>(),
            mappings.device<surfaceNodeMapping::ownerParticleIndexField>(),
            mappings.device<surfaceNodeMapping::particleSurfaceNodeIndexField>()};
}

contactOutputDeviceView makeContactOutputDeviceView(contactContainer& contacts) noexcept
{
    return {contacts.device<contact::masterParticleIndexField>(),
            contacts.device<contact::slaveParticleIndexField>(),
            contacts.device<contact::particleSurfaceNodeIndexField>(),
            contacts.device<contact::pointField>(),
            contacts.device<contact::normalField>(),
            contacts.device<contact::overlapField>(),
            contacts.device<contact::areaField>(),
            contacts.device<contact::effectiveMassField>(),
            contacts.device<contact::effectiveRadiusField>(),
            contacts.device<contact::normalForceMagnitudeField>(),
            contacts.device<contact::slidingSpringDeformationField>(),
            contacts.device<contact::rollingSpringDeformationField>(),
            contacts.device<contact::torsionalSpringDeformationField>()};
}

contactHistoryDeviceView makeContactHistoryDeviceView(interactionContainer& interactions) noexcept
{
    const contactHistory::device_type histories = interactions.contactHistories();
    const neighborRangeHistory::device_type ranges = interactions.neighborRangeHistories();
    return {histories.slaveParticleIndex_, histories.slidingSpringDeformation_, histories.rollingSpringDeformation_, histories.torsionalSpringDeformation_, ranges.prefixSum_};
}

__device__ int rangeBegin(const int* prefixSum, int index) noexcept { return index > 0 ? prefixSum[index - 1] : 0; }

__device__ int spatialGridLevelIndex(math::Real diameter, const spatialGridDeviceView& spatialGrid) noexcept
{
    for (int levelIndex = 0; levelIndex < spatialGrid.levelCount_ - 1; ++levelIndex)
    {
        if (diameter <= spatialGrid.levels_[levelIndex].maximumDiameter_)
        {
            return levelIndex;
        }
    }
    return spatialGrid.levelCount_ - 1;
}

__device__ bool spatialGridRange(int3& cellPositionBegin, int3& cellPositionEnd, const math::Vec3& position, math::Real queryRadius, const spatialGridDeviceView& spatialGrid, int levelIndex) noexcept
{
    const spatialGridLevel& level = spatialGrid.levels_[levelIndex];
    const math::Vec3 local = math::hadamardProduct(position - spatialGrid.minimumBoundary_, level.inverseCellSize_);
    const math::Vec3 extent = (queryRadius + 0.5 * level.maximumDiameter_) * level.inverseCellSize_;
    cellPositionBegin =
        make_int3(static_cast<int>(math::detail::floor(local.x - extent.x)), static_cast<int>(math::detail::floor(local.y - extent.y)), static_cast<int>(math::detail::floor(local.z - extent.z)));
    cellPositionEnd =
        make_int3(static_cast<int>(math::detail::floor(local.x + extent.x)), static_cast<int>(math::detail::floor(local.y + extent.y)), static_cast<int>(math::detail::floor(local.z + extent.z)));
    cellPositionBegin.x = cellPositionBegin.x < 0 ? 0 : (cellPositionBegin.x < level.size_.x ? cellPositionBegin.x : level.size_.x - 1);
    cellPositionBegin.y = cellPositionBegin.y < 0 ? 0 : (cellPositionBegin.y < level.size_.y ? cellPositionBegin.y : level.size_.y - 1);
    cellPositionBegin.z = cellPositionBegin.z < 0 ? 0 : (cellPositionBegin.z < level.size_.z ? cellPositionBegin.z : level.size_.z - 1);
    cellPositionEnd.x = cellPositionEnd.x < 0 ? 0 : (cellPositionEnd.x < level.size_.x ? cellPositionEnd.x : level.size_.x - 1);
    cellPositionEnd.y = cellPositionEnd.y < 0 ? 0 : (cellPositionEnd.y < level.size_.y ? cellPositionEnd.y : level.size_.y - 1);
    cellPositionEnd.z = cellPositionEnd.z < 0 ? 0 : (cellPositionEnd.z < level.size_.z ? cellPositionEnd.z : level.size_.z - 1);
    return true;
}

__device__ int findSpatialGridCell(const spatialGridDeviceView& spatialGrid, int levelIndex, int hash) noexcept
{
    int begin = spatialGrid.levels_[levelIndex].occupiedCellBegin_;
    int end = spatialGrid.levels_[levelIndex].occupiedCellEnd_;
    while (begin < end)
    {
        const int middle = begin + (end - begin) / 2;
        if (spatialGrid.cellHash_[middle] < hash)
        {
            begin = middle + 1;
        }
        else
        {
            end = middle;
        }
    }
    return begin < spatialGrid.levels_[levelIndex].occupiedCellEnd_ && spatialGrid.cellHash_[begin] == hash ? begin : -1;
}

class spatialGridSearch
{
public:
    __device__ spatialGridSearch(const spatialGridDeviceView& spatialGrid, const math::Vec3& position, math::Real queryRadius, int firstLevel) noexcept
        : spatialGrid_(spatialGrid), position_(position), queryRadius_(queryRadius), nextLevel_(firstLevel)
    {}

    __device__ bool next(int& particleIndex, int& levelIndex) noexcept
    {
        for (;;)
        {
            if (sortedIndex_ < sortedEnd_)
            {
                particleIndex = spatialGrid_.particleIndex_[sortedIndex_++];
                levelIndex = currentLevel_;
                return true;
            }
            if (loadNextOccupiedCell())
            {
                continue;
            }
            if (!beginNextLevel())
            {
                return false;
            }
        }
    }

private:
    __device__ bool beginNextLevel() noexcept
    {
        while (nextLevel_ < spatialGrid_.levelCount_)
        {
            currentLevel_ = nextLevel_++;
            if (spatialGridRange(cellPositionBegin_, cellPositionEnd_, position_, queryRadius_, spatialGrid_, currentLevel_))
            {
                cellPosition_ = cellPositionBegin_;
                cellPositionValid_ = true;
                return true;
            }
        }
        return false;
    }

    __device__ bool loadNextOccupiedCell() noexcept
    {
        const spatialGridLevel& level = spatialGrid_.levels_[currentLevel_];
        while (cellPositionValid_)
        {
            const int hash = math::linearIndex(cellPosition_.x, cellPosition_.y, cellPosition_.z, level.size_.x, level.size_.y);
            advanceCellPosition();
            const int cellIndex = findSpatialGridCell(spatialGrid_, currentLevel_, hash);
            if (cellIndex >= 0)
            {
                sortedIndex_ = spatialGrid_.cellSortedParticleIndexBegin_[cellIndex];
                sortedEnd_ = spatialGrid_.cellSortedParticleIndexEnd_[cellIndex];
                return true;
            }
        }
        return false;
    }

    __device__ void advanceCellPosition() noexcept
    {
        if (cellPosition_.x < cellPositionEnd_.x)
        {
            ++cellPosition_.x;
            return;
        }
        cellPosition_.x = cellPositionBegin_.x;
        if (cellPosition_.y < cellPositionEnd_.y)
        {
            ++cellPosition_.y;
            return;
        }
        cellPosition_.y = cellPositionBegin_.y;
        if (cellPosition_.z < cellPositionEnd_.z)
        {
            ++cellPosition_.z;
            return;
        }
        cellPositionValid_ = false;
    }

    const spatialGridDeviceView& spatialGrid_;
    math::Vec3 position_;
    math::Real queryRadius_{0.0};
    int nextLevel_{0};
    int currentLevel_{0};
    int3 cellPositionBegin_{0, 0, 0};
    int3 cellPositionEnd_{0, 0, 0};
    int3 cellPosition_{0, 0, 0};
    int sortedIndex_{0};
    int sortedEnd_{0};
    bool cellPositionValid_{false};
};

__device__ void restoreContactHistory(math::Vec3& slidingSpringDeformation,
                                      math::Vec3& rollingSpringDeformation,
                                      math::Vec3& torsionalSpringDeformation,
                                      const contactHistoryDeviceView& histories,
                                      int masterParticleOrSurfaceNodeIndex,
                                      int slaveParticleIndex,
                                      bool restoreRotational) noexcept
{
    slidingSpringDeformation = math::Vec3::zero();
    rollingSpringDeformation = math::Vec3::zero();
    torsionalSpringDeformation = math::Vec3::zero();
    const int historyBegin = rangeBegin(histories.prefixSum_, masterParticleOrSurfaceNodeIndex);
    const int historyEnd = histories.prefixSum_[masterParticleOrSurfaceNodeIndex];
    for (int historyIndex = historyBegin; historyIndex < historyEnd; ++historyIndex)
    {
        if (histories.slaveParticleIndex_[historyIndex] != slaveParticleIndex)
        {
            continue;
        }
        slidingSpringDeformation = histories.slidingSpringDeformation_[historyIndex];
        if (restoreRotational)
        {
            rollingSpringDeformation = histories.rollingSpringDeformation_[historyIndex];
            torsionalSpringDeformation = histories.torsionalSpringDeformation_[historyIndex];
        }
        return;
    }
}

__device__ bool detectSpherePair(math::Vec3& point,
                                 math::Vec3& normal,
                                 math::Real& overlap,
                                 math::Real& area,
                                 math::Real& effectiveRadius,
                                 const particleDeviceView& spheres,
                                 int masterIndex,
                                 int slaveIndex) noexcept
{
    if (spheres.inverseMass_[masterIndex] == 0.0 && spheres.inverseMass_[slaveIndex] == 0.0)
    {
        return false;
    }
    return execution::detectSphereContact(point,
                                          normal,
                                          overlap,
                                          area,
                                          effectiveRadius,
                                          spheres.position_[masterIndex],
                                          spheres.radius_[masterIndex],
                                          spheres.position_[slaveIndex],
                                          spheres.radius_[slaveIndex]);
}

__device__ bool detectSphereLevelSetPair(math::Real& overlap,
                                         math::Vec3& normal,
                                         const particleDeviceView& sphere,
                                         int sphereIndex,
                                         const levelSetParticleDeviceView& levelSet,
                                         const levelSetGeometryDeviceView& geometries,
                                         int levelSetIndex) noexcept
{
    const int geometryIndex = levelSet.geometryIndex_[levelSetIndex];
    return execution::detectLevelSetContact(overlap,
                                            normal,
                                            geometries.gridNodes_.signedDistance_,
                                            sphere.position_[sphereIndex],
                                            levelSet.body_.position_[levelSetIndex],
                                            levelSet.body_.orientation_[levelSetIndex],
                                            geometries.descriptors_.gridNodeOrigin_[geometryIndex],
                                            geometries.descriptors_.gridNodeInverseSpacing_[geometryIndex],
                                            geometries.descriptors_.gridNodeSize_[geometryIndex],
                                            geometries.descriptors_.signedDistanceOffset_[geometryIndex],
                                            sphere.radius_[sphereIndex]);
}

__device__ bool detectSurfaceNodePair(math::Real& overlap,
                                      math::Vec3& normal,
                                      const math::Vec3& surfaceNodePosition,
                                      const levelSetParticleDeviceView& slave,
                                      const levelSetGeometryDeviceView& geometries,
                                      int slaveIndex) noexcept
{
    if (math::distanceSquared(surfaceNodePosition, slave.body_.position_[slaveIndex]) > slave.body_.radius_[slaveIndex] * slave.body_.radius_[slaveIndex])
    {
        return false;
    }
    const int geometryIndex = slave.geometryIndex_[slaveIndex];
    return execution::detectLevelSetContact(overlap,
                                            normal,
                                            geometries.gridNodes_.signedDistance_,
                                            surfaceNodePosition,
                                            slave.body_.position_[slaveIndex],
                                            slave.body_.orientation_[slaveIndex],
                                            geometries.descriptors_.gridNodeOrigin_[geometryIndex],
                                            geometries.descriptors_.gridNodeInverseSpacing_[geometryIndex],
                                            geometries.descriptors_.gridNodeSize_[geometryIndex],
                                            geometries.descriptors_.signedDistanceOffset_[geometryIndex]);
}

__global__ void countSphereContactsKernel(neighborRange::device_type ranges, particleDeviceView spheres, spatialGridDeviceView spatialGrid, int particleCount)
{
    const int masterIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (masterIndex >= particleCount)
    {
        return;
    }
    if (spheres.inverseMass_[masterIndex] == 0.0)
    {
        ranges.count_[masterIndex] = 0;
        return;
    }

    const math::Real masterRadius = spheres.radius_[masterIndex];
    const int firstLevel = spatialGridLevelIndex(2.0 * masterRadius, spatialGrid);
    spatialGridSearch search(spatialGrid, spheres.position_[masterIndex], masterRadius, firstLevel);
    int contactCount = 0;
    int slaveIndex;
    int slaveLevel;
    while (search.next(slaveIndex, slaveLevel))
    {
        if (slaveLevel == firstLevel && masterIndex >= slaveIndex)
        {
            continue;
        }
        math::Vec3 point;
        math::Vec3 normal;
        math::Real overlap;
        math::Real area;
        math::Real effectiveRadius;
        contactCount += detectSpherePair(point, normal, overlap, area, effectiveRadius, spheres, masterIndex, slaveIndex) ? 1 : 0;
    }
    ranges.count_[masterIndex] = contactCount;
}

__global__ void writeSphereContactsKernel(contactOutputDeviceView contacts,
                                          contactHistoryDeviceView histories,
                                          neighborRange::device_type ranges,
                                          particleDeviceView spheres,
                                          spatialGridDeviceView spatialGrid,
                                          int particleCount)
{
    const int masterIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (masterIndex >= particleCount)
    {
        return;
    }
    if (spheres.inverseMass_[masterIndex] == 0.0)
    {
        return;
    }

    const math::Real masterRadius = spheres.radius_[masterIndex];
    const int firstLevel = spatialGridLevelIndex(2.0 * masterRadius, spatialGrid);
    spatialGridSearch search(spatialGrid, spheres.position_[masterIndex], masterRadius, firstLevel);
    const int writeOffset = rangeBegin(ranges.prefixSum_, masterIndex);
    int localContactIndex = 0;
    int slaveIndex;
    int slaveLevel;
    while (search.next(slaveIndex, slaveLevel))
    {
        if (slaveLevel == firstLevel && masterIndex >= slaveIndex)
        {
            continue;
        }

        math::Vec3 point;
        math::Vec3 normal;
        math::Real overlap;
        math::Real area;
        math::Real effectiveRadius;
        if (!detectSpherePair(point, normal, overlap, area, effectiveRadius, spheres, masterIndex, slaveIndex))
        {
            continue;
        }

        const int writeIndex = writeOffset + localContactIndex++;
        contacts.masterParticleIndex_[writeIndex] = masterIndex;
        contacts.slaveParticleIndex_[writeIndex] = slaveIndex;
        contacts.particleSurfaceNodeIndex_[writeIndex] = -1;
        contacts.point_[writeIndex] = point;
        contacts.normal_[writeIndex] = normal;
        contacts.overlap_[writeIndex] = overlap;
        contacts.area_[writeIndex] = area;
        contacts.effectiveMass_[writeIndex] = execution::effectiveMass(spheres.inverseMass_[masterIndex], spheres.inverseMass_[slaveIndex]);
        contacts.effectiveRadius_[writeIndex] = effectiveRadius;
        contacts.normalForceMagnitude_[writeIndex] = 0.0;
        restoreContactHistory(contacts.slidingSpringDeformation_[writeIndex],
                              contacts.rollingSpringDeformation_[writeIndex],
                              contacts.torsionalSpringDeformation_[writeIndex],
                              histories,
                              masterIndex,
                              slaveIndex,
                              true);
    }
}

__global__ void countSphereLevelSetContactsKernel(neighborRange::device_type ranges,
                                                  particleDeviceView masterSpheres,
                                                  levelSetParticleDeviceView slaveLSParticles,
                                                  levelSetGeometryDeviceView geometries,
                                                  spatialGridDeviceView slaveSpatialGrid,
                                                  int masterCount)
{
    const int masterIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (masterIndex >= masterCount)
    {
        return;
    }

    const math::Real masterRadius = masterSpheres.radius_[masterIndex];
    spatialGridSearch search(slaveSpatialGrid, masterSpheres.position_[masterIndex], masterRadius, 0);
    int contactCount = 0;
    int slaveIndex;
    int slaveLevel;
    while (search.next(slaveIndex, slaveLevel))
    {
        math::Real overlap;
        math::Vec3 normal;
        contactCount += detectSphereLevelSetPair(overlap, normal, masterSpheres, masterIndex, slaveLSParticles, geometries, slaveIndex) ? 1 : 0;
    }
    ranges.count_[masterIndex] = contactCount;
}

__global__ void writeSphereLevelSetContactsKernel(contactOutputDeviceView contacts,
                                                  contactHistoryDeviceView histories,
                                                  neighborRange::device_type ranges,
                                                  particleDeviceView masterSpheres,
                                                  levelSetParticleDeviceView slaveLSParticles,
                                                  levelSetGeometryDeviceView geometries,
                                                  spatialGridDeviceView slaveSpatialGrid,
                                                  int masterCount)
{
    const int masterIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (masterIndex >= masterCount)
    {
        return;
    }

    const math::Real masterRadius = masterSpheres.radius_[masterIndex];
    spatialGridSearch search(slaveSpatialGrid, masterSpheres.position_[masterIndex], masterRadius, 0);
    const int writeOffset = rangeBegin(ranges.prefixSum_, masterIndex);
    int localContactIndex = 0;
    int slaveIndex;
    int slaveLevel;
    while (search.next(slaveIndex, slaveLevel))
    {
        math::Real overlap;
        math::Vec3 normal;
        if (!detectSphereLevelSetPair(overlap, normal, masterSpheres, masterIndex, slaveLSParticles, geometries, slaveIndex))
        {
            continue;
        }

        const int writeIndex = writeOffset + localContactIndex++;
        contacts.masterParticleIndex_[writeIndex] = masterIndex;
        contacts.slaveParticleIndex_[writeIndex] = slaveIndex;
        contacts.particleSurfaceNodeIndex_[writeIndex] = -1;
        contacts.point_[writeIndex] = execution::sphereLevelSetContactPoint(masterSpheres.position_[masterIndex], masterRadius, normal, overlap);
        contacts.normal_[writeIndex] = normal;
        contacts.overlap_[writeIndex] = overlap;
        contacts.area_[writeIndex] = math::pi * masterRadius * overlap;
        contacts.effectiveMass_[writeIndex] = execution::effectiveMass(masterSpheres.inverseMass_[masterIndex], slaveLSParticles.body_.inverseMass_[slaveIndex]);
        // Physical radius used by the unified spherical contact model.
        contacts.effectiveRadius_[writeIndex] = masterRadius;
        contacts.normalForceMagnitude_[writeIndex] = 0.0;
        restoreContactHistory(contacts.slidingSpringDeformation_[writeIndex],
                              contacts.rollingSpringDeformation_[writeIndex],
                              contacts.torsionalSpringDeformation_[writeIndex],
                              histories,
                              masterIndex,
                              slaveIndex,
                              true);
    }
}

__global__ void countLevelSetNeighborsKernel(neighborRange::device_type ranges, levelSetParticleDeviceView particles, spatialGridDeviceView spatialGrid, int particleCount)
{
    const int masterIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (masterIndex >= particleCount)
    {
        return;
    }
    if (particles.body_.inverseMass_[masterIndex] == 0.0)
    {
        ranges.count_[masterIndex] = 0;
        return;
    }

    const math::Real masterRadius = particles.body_.radius_[masterIndex];
    const int firstLevel = spatialGridLevelIndex(2.0 * masterRadius, spatialGrid);
    spatialGridSearch search(spatialGrid, particles.body_.position_[masterIndex], masterRadius, firstLevel);
    int neighborCount = 0;
    int slaveIndex;
    int slaveLevel;
    while (search.next(slaveIndex, slaveLevel))
    {
        if (slaveLevel == firstLevel && masterIndex >= slaveIndex)
        {
            continue;
        }
        neighborCount += execution::boundingSpheresOverlap(particles.body_.position_[masterIndex], masterRadius, particles.body_.position_[slaveIndex], particles.body_.radius_[slaveIndex]) ? 1 : 0;
    }
    ranges.count_[masterIndex] = neighborCount;
}

__global__ void writeLevelSetNeighborsKernel(neighborRange::device_type ranges,
                                             particleNeighbor::device_type neighbors,
                                             levelSetParticleDeviceView particles,
                                             spatialGridDeviceView spatialGrid,
                                             int particleCount)
{
    const int masterIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (masterIndex >= particleCount)
    {
        return;
    }
    if (particles.body_.inverseMass_[masterIndex] == 0.0)
    {
        return;
    }

    const math::Real masterRadius = particles.body_.radius_[masterIndex];
    const int firstLevel = spatialGridLevelIndex(2.0 * masterRadius, spatialGrid);
    spatialGridSearch search(spatialGrid, particles.body_.position_[masterIndex], masterRadius, firstLevel);
    const int writeOffset = rangeBegin(ranges.prefixSum_, masterIndex);
    int localNeighborIndex = 0;
    int slaveIndex;
    int slaveLevel;
    while (search.next(slaveIndex, slaveLevel))
    {
        if (slaveLevel == firstLevel && masterIndex >= slaveIndex)
        {
            continue;
        }
        if (!execution::boundingSpheresOverlap(particles.body_.position_[masterIndex], masterRadius, particles.body_.position_[slaveIndex], particles.body_.radius_[slaveIndex]))
        {
            continue;
        }
        const int writeIndex = writeOffset + localNeighborIndex++;
        neighbors.particleIndex_[writeIndex] = slaveIndex;
        neighbors.patchArea_[writeIndex] = 0.0;
    }
}

__global__ void countLevelSetContactsKernel(neighborRange::device_type surfaceNodeNeighborRanges,
                                            surfaceNodeMappingDeviceView mappings,
                                            levelSetParticleDeviceView particles,
                                            neighborRange::device_type particleNeighborRanges,
                                            particleNeighbor::device_type particleNeighbors,
                                            levelSetGeometryDeviceView geometries,
                                            int surfaceNodeCount)
{
    const int surfaceNodeIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (surfaceNodeIndex >= surfaceNodeCount)
    {
        return;
    }

    const int masterIndex = mappings.ownerParticleIndex_[surfaceNodeIndex];
    const int geometrySurfaceNodeIndex = mappings.geometrySurfaceNodeIndex_[surfaceNodeIndex];
    const math::Vec3 surfaceNodePosition =
        particles.body_.position_[masterIndex] + math::rotateUnit(particles.body_.orientation_[masterIndex], geometries.surfaceNodes_.position_[geometrySurfaceNodeIndex]);
    const int neighborBegin = rangeBegin(particleNeighborRanges.prefixSum_, masterIndex);
    const int neighborEnd = particleNeighborRanges.prefixSum_[masterIndex];
    int contactCount = 0;
    for (int neighborIndex = neighborBegin; neighborIndex < neighborEnd; ++neighborIndex)
    {
        math::Real overlap;
        math::Vec3 normal;
        contactCount += detectSurfaceNodePair(overlap, normal, surfaceNodePosition, particles, geometries, particleNeighbors.particleIndex_[neighborIndex]) ? 1 : 0;
    }
    surfaceNodeNeighborRanges.count_[surfaceNodeIndex] = contactCount;
}

__global__ void writeLevelSetContactsKernel(contactOutputDeviceView contacts,
                                            contactHistoryDeviceView histories,
                                            neighborRange::device_type surfaceNodeNeighborRanges,
                                            surfaceNodeMappingDeviceView mappings,
                                            levelSetParticleDeviceView particles,
                                            neighborRange::device_type particleNeighborRanges,
                                            particleNeighbor::device_type particleNeighbors,
                                            levelSetGeometryDeviceView geometries,
                                            int surfaceNodeCount)
{
    const int surfaceNodeIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (surfaceNodeIndex >= surfaceNodeCount)
    {
        return;
    }

    const int masterIndex = mappings.ownerParticleIndex_[surfaceNodeIndex];
    const int geometrySurfaceNodeIndex = mappings.geometrySurfaceNodeIndex_[surfaceNodeIndex];
    const math::Vec3 surfaceNodePosition =
        particles.body_.position_[masterIndex] + math::rotateUnit(particles.body_.orientation_[masterIndex], geometries.surfaceNodes_.position_[geometrySurfaceNodeIndex]);
    const math::Real surfaceNodeArea = geometries.surfaceNodes_.area_[geometrySurfaceNodeIndex];
    const int neighborBegin = rangeBegin(particleNeighborRanges.prefixSum_, masterIndex);
    const int neighborEnd = particleNeighborRanges.prefixSum_[masterIndex];
    const int writeOffset = rangeBegin(surfaceNodeNeighborRanges.prefixSum_, surfaceNodeIndex);
    int localContactIndex = 0;
    for (int neighborIndex = neighborBegin; neighborIndex < neighborEnd; ++neighborIndex)
    {
        const int slaveIndex = particleNeighbors.particleIndex_[neighborIndex];
        math::Real overlap;
        math::Vec3 normal;
        if (!detectSurfaceNodePair(overlap, normal, surfaceNodePosition, particles, geometries, slaveIndex))
        {
            continue;
        }

        const int writeIndex = writeOffset + localContactIndex++;
        contacts.masterParticleIndex_[writeIndex] = masterIndex;
        contacts.slaveParticleIndex_[writeIndex] = slaveIndex;
        contacts.particleSurfaceNodeIndex_[writeIndex] = mappings.particleSurfaceNodeIndex_[surfaceNodeIndex];
        contacts.point_[writeIndex] = execution::surfaceNodeLevelSetContactPoint(surfaceNodePosition, normal, overlap);
        contacts.normal_[writeIndex] = normal;
        contacts.overlap_[writeIndex] = overlap;
        contacts.area_[writeIndex] = surfaceNodeArea;
        contacts.effectiveMass_[writeIndex] = 0.0;
        contacts.effectiveRadius_[writeIndex] = 0.0;
        contacts.normalForceMagnitude_[writeIndex] = 0.0;
        restoreContactHistory(contacts.slidingSpringDeformation_[writeIndex],
                              contacts.rollingSpringDeformation_[writeIndex],
                              contacts.torsionalSpringDeformation_[writeIndex],
                              histories,
                              surfaceNodeIndex,
                              slaveIndex,
                              false);
        detail::atomicAddReal(&particleNeighbors.patchArea_[neighborIndex], surfaceNodeArea);
    }
}

__global__ void levelSetContactEffectiveMassKernel(contactOutputDeviceView contacts,
                                                   neighborRange::device_type particleNeighborRanges,
                                                   particleNeighbor::device_type particleNeighbors,
                                                   particleDeviceView particles,
                                                   int contactCount)
{
    const int contactIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (contactIndex >= contactCount)
    {
        return;
    }

    const int masterIndex = contacts.masterParticleIndex_[contactIndex];
    const int slaveIndex = contacts.slaveParticleIndex_[contactIndex];
    const int neighborBegin = rangeBegin(particleNeighborRanges.prefixSum_, masterIndex);
    const int neighborEnd = particleNeighborRanges.prefixSum_[masterIndex];
    math::Real patchArea = 0.0;
    for (int neighborIndex = neighborBegin; neighborIndex < neighborEnd; ++neighborIndex)
    {
        if (particleNeighbors.particleIndex_[neighborIndex] == slaveIndex)
        {
            patchArea = particleNeighbors.patchArea_[neighborIndex];
            break;
        }
    }

    const math::Real pairEffectiveMass = execution::effectiveMass(particles.inverseMass_[masterIndex], particles.inverseMass_[slaveIndex]);
    contacts.effectiveMass_[contactIndex] = patchArea > 0.0 ? pairEffectiveMass * contacts.area_[contactIndex] / patchArea : 0.0;
}

void launchCountSphereContacts(neighborRange::device_type ranges, particleDeviceView spheres, spatialGridDeviceView spatialGrid, int particleCount, cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }
    const int blockCount = (particleCount + blockSize - 1) / blockSize;
    countSphereContactsKernel<<<blockCount, blockSize, 0, stream>>>(ranges, spheres, spatialGrid, particleCount);
    host_device_detail::checkCuda(cudaGetLastError(), "countSphereContactsKernel launch");
}

void launchWriteSphereContacts(contactOutputDeviceView contacts,
                               contactHistoryDeviceView histories,
                               neighborRange::device_type ranges,
                               particleDeviceView spheres,
                               spatialGridDeviceView spatialGrid,
                               int particleCount,
                               cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }
    const int blockCount = (particleCount + blockSize - 1) / blockSize;
    writeSphereContactsKernel<<<blockCount, blockSize, 0, stream>>>(contacts, histories, ranges, spheres, spatialGrid, particleCount);
    host_device_detail::checkCuda(cudaGetLastError(), "writeSphereContactsKernel launch");
}

void launchCountLevelSetNeighbors(neighborRange::device_type ranges, levelSetParticleDeviceView particles, spatialGridDeviceView spatialGrid, int particleCount, cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }
    const int blockCount = (particleCount + blockSize - 1) / blockSize;
    countLevelSetNeighborsKernel<<<blockCount, blockSize, 0, stream>>>(ranges, particles, spatialGrid, particleCount);
    host_device_detail::checkCuda(cudaGetLastError(), "countLevelSetNeighborsKernel launch");
}

void launchWriteLevelSetNeighbors(neighborRange::device_type ranges,
                                  particleNeighbor::device_type neighbors,
                                  levelSetParticleDeviceView particles,
                                  spatialGridDeviceView spatialGrid,
                                  int particleCount,
                                  cudaStream_t stream)
{
    if (particleCount <= 0)
    {
        return;
    }
    const int blockCount = (particleCount + blockSize - 1) / blockSize;
    writeLevelSetNeighborsKernel<<<blockCount, blockSize, 0, stream>>>(ranges, neighbors, particles, spatialGrid, particleCount);
    host_device_detail::checkCuda(cudaGetLastError(), "writeLevelSetNeighborsKernel launch");
}

void launchCountLevelSetContacts(interactionContainer& interactions,
                                 const surfaceNodeMappingContainer& mappings,
                                 const LSParticleContainer& particles,
                                 const LSGeometryContainer& geometries,
                                 cudaStream_t stream)
{
    const int surfaceNodeCount = static_cast<int>(mappings.deviceSize());
    if (surfaceNodeCount <= 0)
    {
        return;
    }
    const int blockCount = (surfaceNodeCount + blockSize - 1) / blockSize;
    countLevelSetContactsKernel<<<blockCount, blockSize, 0, stream>>>(interactions.surfaceNodeNeighborRanges(),
                                                                      makeSurfaceNodeMappingDeviceView(mappings),
                                                                      makeLevelSetParticleDeviceView(particles),
                                                                      interactions.particleNeighborRanges(),
                                                                      interactions.particleNeighborList(),
                                                                      makeLevelSetGeometryDeviceView(geometries),
                                                                      surfaceNodeCount);
    host_device_detail::checkCuda(cudaGetLastError(), "countLevelSetContactsKernel launch");
}

void launchWriteLevelSetContacts(interactionContainer& interactions,
                                 const surfaceNodeMappingContainer& mappings,
                                 const LSParticleContainer& particles,
                                 const LSGeometryContainer& geometries,
                                 cudaStream_t stream)
{
    const int surfaceNodeCount = static_cast<int>(mappings.deviceSize());
    if (surfaceNodeCount <= 0)
    {
        return;
    }
    const int blockCount = (surfaceNodeCount + blockSize - 1) / blockSize;
    writeLevelSetContactsKernel<<<blockCount, blockSize, 0, stream>>>(makeContactOutputDeviceView(interactions.contacts()),
                                                                      makeContactHistoryDeviceView(interactions),
                                                                      interactions.surfaceNodeNeighborRanges(),
                                                                      makeSurfaceNodeMappingDeviceView(mappings),
                                                                      makeLevelSetParticleDeviceView(particles),
                                                                      interactions.particleNeighborRanges(),
                                                                      interactions.particleNeighborList(),
                                                                      makeLevelSetGeometryDeviceView(geometries),
                                                                      surfaceNodeCount);
    host_device_detail::checkCuda(cudaGetLastError(), "writeLevelSetContactsKernel launch");
}

void launchLevelSetContactEffectiveMass(interactionContainer& interactions, const LSParticleContainer& particles, cudaStream_t stream)
{
    const int contactCount = static_cast<int>(interactions.contacts().deviceSize());
    if (contactCount <= 0)
    {
        return;
    }
    const int blockCount = (contactCount + blockSize - 1) / blockSize;
    levelSetContactEffectiveMassKernel<<<blockCount, blockSize, 0, stream>>>(makeContactOutputDeviceView(interactions.contacts()),
                                                                             interactions.particleNeighborRanges(),
                                                                             interactions.particleNeighborList(),
                                                                             makeParticleDeviceView(particles),
                                                                             contactCount);
    host_device_detail::checkCuda(cudaGetLastError(), "levelSetContactEffectiveMassKernel launch");
}

void launchCountSphereLevelSetContacts(interactionContainer& interactions,
                                       const particleContainer& masterSpheres,
                                       const LSParticleContainer& slaveLSParticles,
                                       const LSGeometryContainer& geometries,
                                       const spatialGridDeviceView& slaveSpatialGrid,
                                       cudaStream_t stream)
{
    const int masterCount = static_cast<int>(masterSpheres.deviceSize());
    if (masterCount <= 0)
    {
        return;
    }
    const int blockCount = (masterCount + blockSize - 1) / blockSize;
    countSphereLevelSetContactsKernel<<<blockCount, blockSize, 0, stream>>>(interactions.particleNeighborRanges(),
                                                                            makeParticleDeviceView(masterSpheres),
                                                                            makeLevelSetParticleDeviceView(slaveLSParticles),
                                                                            makeLevelSetGeometryDeviceView(geometries),
                                                                            slaveSpatialGrid,
                                                                            masterCount);
    host_device_detail::checkCuda(cudaGetLastError(), "countSphereLevelSetContactsKernel launch");
}

void launchWriteSphereLevelSetContacts(interactionContainer& interactions,
                                       const particleContainer& masterSpheres,
                                       const LSParticleContainer& slaveLSParticles,
                                       const LSGeometryContainer& geometries,
                                       const spatialGridDeviceView& slaveSpatialGrid,
                                       cudaStream_t stream)
{
    const int masterCount = static_cast<int>(masterSpheres.deviceSize());
    if (masterCount <= 0)
    {
        return;
    }
    const int blockCount = (masterCount + blockSize - 1) / blockSize;
    writeSphereLevelSetContactsKernel<<<blockCount, blockSize, 0, stream>>>(makeContactOutputDeviceView(interactions.contacts()),
                                                                            makeContactHistoryDeviceView(interactions),
                                                                            interactions.particleNeighborRanges(),
                                                                            makeParticleDeviceView(masterSpheres),
                                                                            makeLevelSetParticleDeviceView(slaveLSParticles),
                                                                            makeLevelSetGeometryDeviceView(geometries),
                                                                            slaveSpatialGrid,
                                                                            masterCount);
    host_device_detail::checkCuda(cudaGetLastError(), "writeSphereLevelSetContactsKernel launch");
}

} // namespace

void launchSphereSphereContactDetection(interactionContainer& interactions, const particleContainer& spheres, const spatialGridDeviceView& spatialGrid, cudaStream_t stream)
{
    interactions.saveCurrentStepToHistory(stream);
    const particleDeviceView sphereView = makeParticleDeviceView(spheres);
    const int particleCount = static_cast<int>(spheres.deviceSize());
    launchCountSphereContacts(interactions.particleNeighborRanges(), sphereView, spatialGrid, particleCount, stream);
    interactions.contactCount(stream);
    launchWriteSphereContacts(makeContactOutputDeviceView(interactions.contacts()),
                              makeContactHistoryDeviceView(interactions),
                              interactions.particleNeighborRanges(),
                              sphereView,
                              spatialGrid,
                              particleCount,
                              stream);
}

void launchSphereLevelSetContactDetection(interactionContainer& interactions,
                                          const particleContainer& masterSpheres,
                                          const LSParticleContainer& slaveLSParticles,
                                          const LSGeometryContainer& geometries,
                                          const spatialGridDeviceView& slaveSpatialGrid,
                                          cudaStream_t stream)
{
    interactions.saveCurrentStepToHistory(stream);
    launchCountSphereLevelSetContacts(interactions, masterSpheres, slaveLSParticles, geometries, slaveSpatialGrid, stream);
    interactions.contactCount(stream);
    launchWriteSphereLevelSetContacts(interactions, masterSpheres, slaveLSParticles, geometries, slaveSpatialGrid, stream);
}

void launchLevelSetLevelSetContactDetection(interactionContainer& interactions,
                                            const LSParticleContainer& LSParticles,
                                            const LSGeometryContainer& geometries,
                                            const spatialGridDeviceView& spatialGrid,
                                            cudaStream_t stream)
{
    interactions.saveCurrentStepToHistory(stream);
    const surfaceNodeMappingContainer& surfaceNodeMappings = interactions.surfaceNodeMappings();

    const levelSetParticleDeviceView particleView = makeLevelSetParticleDeviceView(LSParticles);
    const int particleCount = static_cast<int>(LSParticles.deviceSize());
    launchCountLevelSetNeighbors(interactions.particleNeighborRanges(), particleView, spatialGrid, particleCount, stream);
    interactions.resizeParticleNeighborList(stream);
    launchWriteLevelSetNeighbors(interactions.particleNeighborRanges(), interactions.particleNeighborList(), particleView, spatialGrid, particleCount, stream);

    launchCountLevelSetContacts(interactions, surfaceNodeMappings, LSParticles, geometries, stream);
    interactions.contactCount(stream);
    launchWriteLevelSetContacts(interactions, surfaceNodeMappings, LSParticles, geometries, stream);
    launchLevelSetContactEffectiveMass(interactions, LSParticles, stream);
}

} // namespace fundem::cuda
