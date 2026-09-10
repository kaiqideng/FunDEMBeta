#include "sphInteractionKernel.cuh"

#include "atomicKernel.cuh"
#include "data/HostAoSDeviceSoA.h"
#include "execution/motionIntegration.h"
#include "execution/sphFunctions.h"
#include "execution/sphCouplingFunctions.h"

#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/transform_reduce.h>

#include <limits>

namespace fundem::cuda
{
namespace
{

constexpr int blockSize = 256;
constexpr math::Real invalidNeighborDisplacement = std::numeric_limits<math::Real>::infinity();

struct maximumNeighborDisplacement {
    __host__ __device__ math::Real operator()(math::Real first, math::Real second) const noexcept { return first > second ? first : second; }
};

struct readNeighborDisplacement {
    const math::Vec3* fluid_;
    const math::Vec3* boundary_;
    const math::Vec3* reference_;
    int fluidCount_;
    __host__ __device__ math::Real operator()(int index) const noexcept
    {
        const math::Vec3 position = index < fluidCount_ ? fluid_[index] : boundary_[index - fluidCount_];
        const math::Vec3 displacement = position - reference_[index];
        const math::Real squared = math::dot(displacement, displacement);
        return math::isFinite(squared) ? squared : invalidNeighborDisplacement;
    }
};

struct readSPHState {
    const math::Vec3* velocity_{nullptr};
    const math::Vec3* force_{nullptr};
    const math::Real* inverseMass_{nullptr};
    const math::Real* density_{nullptr};
    const math::Real* pressure_{nullptr};
    const math::Real* densityRate_{nullptr};
    const math::Vec3* position_{nullptr};
    math::Vec3 gravity_{math::Vec3::zero()};

    __host__ __device__ SPHStateStatistics operator()(int index) const noexcept
    {
        const math::Vec3 velocity = velocity_[index];
        const math::Vec3 acceleration = inverseMass_[index] * force_[index] + gravity_;
        const math::Real density = density_[index];
        const bool validVelocity = math::isFinite(velocity);
        const bool validAcceleration = math::isFinite(acceleration);
        const bool validDensity = math::isFinite(density);
        return {validVelocity ? math::norm(velocity) : 0.0,
                validAcceleration ? math::norm(acceleration) : 0.0,
                validDensity ? density : 0.0,
                validDensity ? density : 0.0,
                validVelocity && validAcceleration && validDensity && math::isFinite(pressure_[index]) && math::isFinite(densityRate_[index]) && math::isFinite(position_[index]) ? 0 : 1};
    }
};

struct combineSPHStateStatistics {
    __host__ __device__ SPHStateStatistics operator()(const SPHStateStatistics& first, const SPHStateStatistics& second) const noexcept
    {
        return {first.maximumVelocity_ > second.maximumVelocity_ ? first.maximumVelocity_ : second.maximumVelocity_,
                first.maximumAcceleration_ > second.maximumAcceleration_ ? first.maximumAcceleration_ : second.maximumAcceleration_,
                first.minimumDensity_ < second.minimumDensity_ ? first.minimumDensity_ : second.minimumDensity_,
                first.maximumDensity_ > second.maximumDensity_ ? first.maximumDensity_ : second.maximumDensity_,
                first.invalidValueCount_ + second.invalidValueCount_};
    }
};

struct readSPHKinematics {
    const math::Vec3* velocity_{nullptr};
    const math::Vec3* acceleration_{nullptr};

    __host__ __device__ SPHKinematicsStatistics operator()(int index) const noexcept
    {
        const math::Vec3 velocity = velocity_[index];
        const math::Vec3 acceleration = acceleration_[index];
        const bool validVelocity = math::isFinite(velocity);
        const bool validAcceleration = math::isFinite(acceleration);
        return {validVelocity ? math::norm(velocity) : 0.0, validAcceleration ? math::norm(acceleration) : 0.0, validVelocity && validAcceleration ? 0 : 1};
    }
};

struct combineSPHKinematicsStatistics {
    __host__ __device__ SPHKinematicsStatistics operator()(const SPHKinematicsStatistics& first, const SPHKinematicsStatistics& second) const noexcept
    {
        return {first.maximumVelocity_ > second.maximumVelocity_ ? first.maximumVelocity_ : second.maximumVelocity_,
                first.maximumAcceleration_ > second.maximumAcceleration_ ? first.maximumAcceleration_ : second.maximumAcceleration_,
                first.invalidValueCount_ + second.invalidValueCount_};
    }
};

struct SPHParticleDeviceView {
    math::Vec3* position_{nullptr};
    math::Vec3* velocity_{nullptr};
    math::Vec3* force_{nullptr};
    math::Vec3* priorForce_{nullptr};
    const math::Real* inverseMass_{nullptr};
    math::Real* density_{nullptr};
    math::Real* pressure_{nullptr};
    math::Real* densityRate_{nullptr};
    int* freeSurface_{nullptr};
    int* constrained_{nullptr};
};

struct VirtualParticleDeviceView {
    const math::Vec3* localPosition_{nullptr};
    const math::Vec3* localNormal_{nullptr};
    const math::Real* volume_{nullptr};
    const int* ownerLSParticleIndex_{nullptr};
    math::Vec3* position_{nullptr};
    math::Vec3* normal_{nullptr};
    math::Vec3* velocity_{nullptr};
    math::Vec3* acceleration_{nullptr};
    math::Vec3* force_{nullptr};
    math::Vec3* priorForce_{nullptr};
    math::Vec3* couplingTorque_{nullptr};
    math::Vec3* couplingForceIncrement_{nullptr};
    math::Vec3* couplingTorqueIncrement_{nullptr};
    math::Vec3* previousKinematicsVelocity_{nullptr};
};

struct ConstVirtualParticleDeviceView {
    const math::Vec3* localPosition_{nullptr};
    const math::Vec3* localNormal_{nullptr};
    const math::Real* volume_{nullptr};
    const int* ownerLSParticleIndex_{nullptr};
    const math::Vec3* position_{nullptr};
    const math::Vec3* normal_{nullptr};
    const math::Vec3* velocity_{nullptr};
    const math::Vec3* acceleration_{nullptr};
    const math::Vec3* force_{nullptr};
    const math::Vec3* priorForce_{nullptr};
    const math::Vec3* couplingTorque_{nullptr};
    const math::Vec3* couplingForceIncrement_{nullptr};
    const math::Vec3* couplingTorqueIncrement_{nullptr};
    const math::Vec3* previousKinematicsVelocity_{nullptr};
};

struct LSParticleDeviceView {
    const math::Vec3* position_{nullptr};
    const math::Quaternion* orientation_{nullptr};
    math::Vec3* velocity_{nullptr};
    math::Vec3* angularVelocity_{nullptr};
    math::Vec3* force_{nullptr};
    math::Vec3* torque_{nullptr};
    const math::Real* inverseMass_{nullptr};
    const math::Mat3* inertiaTensor_{nullptr};
    const math::Mat3* inverseInertiaTensor_{nullptr};
};

struct ConstLSParticleDeviceView {
    const math::Vec3* position_{nullptr};
    const math::Quaternion* orientation_{nullptr};
    const math::Vec3* velocity_{nullptr};
    const math::Vec3* angularVelocity_{nullptr};
    const math::Vec3* force_{nullptr};
    const math::Vec3* torque_{nullptr};
    const math::Real* inverseMass_{nullptr};
    const math::Mat3* inertiaTensor_{nullptr};
    const math::Mat3* inverseInertiaTensor_{nullptr};
};

SPHParticleDeviceView makeSPHParticleDeviceView(SPHParticleContainer& particles) noexcept
{
    return {particles.device<pointMass::positionField>(),
            particles.device<pointMass::velocityField>(),
            particles.device<pointMass::forceField>(),
            particles.device<SPHParticle::priorForceField>(),
            particles.device<pointMass::inverseMassField>(),
            particles.device<SPHParticle::densityField>(),
            particles.device<SPHParticle::pressureField>(),
            particles.device<SPHParticle::densityRateField>(),
            particles.device<SPHParticle::freeSurfaceField>(),
            particles.device<SPHParticle::constrainedField>()};
}

VirtualParticleDeviceView makeVirtualParticleDeviceView(virtualParticleContainer& particles) noexcept
{
    return {particles.device<virtualParticle::localPositionField>(),
            particles.device<virtualParticle::localNormalField>(),
            particles.device<virtualParticle::volumeField>(),
            particles.device<virtualParticle::ownerLSParticleIndexField>(),
            particles.device<pointMass::positionField>(),
            particles.device<virtualParticle::normalField>(),
            particles.device<pointMass::velocityField>(),
            particles.device<virtualParticle::accelerationField>(),
            particles.device<pointMass::forceField>(),
            particles.device<virtualParticle::priorForceField>(),
            particles.device<virtualParticle::couplingTorqueField>(),
            particles.device<virtualParticle::couplingForceIncrementField>(),
            particles.device<virtualParticle::couplingTorqueIncrementField>(),
            particles.device<virtualParticle::previousKinematicsVelocityField>()};
}

ConstVirtualParticleDeviceView makeVirtualParticleDeviceView(const virtualParticleContainer& particles) noexcept
{
    return {particles.device<virtualParticle::localPositionField>(),
            particles.device<virtualParticle::localNormalField>(),
            particles.device<virtualParticle::volumeField>(),
            particles.device<virtualParticle::ownerLSParticleIndexField>(),
            particles.device<pointMass::positionField>(),
            particles.device<virtualParticle::normalField>(),
            particles.device<pointMass::velocityField>(),
            particles.device<virtualParticle::accelerationField>(),
            particles.device<pointMass::forceField>(),
            particles.device<virtualParticle::priorForceField>(),
            particles.device<virtualParticle::couplingTorqueField>(),
            particles.device<virtualParticle::couplingForceIncrementField>(),
            particles.device<virtualParticle::couplingTorqueIncrementField>(),
            particles.device<virtualParticle::previousKinematicsVelocityField>()};
}

LSParticleDeviceView makeLSParticleDeviceView(LSParticleContainer& particles) noexcept
{
    return {particles.device<rigidBody::positionField>(),
            particles.device<rigidBody::orientationField>(),
            particles.device<rigidBody::velocityField>(),
            particles.device<rigidBody::angularVelocityField>(),
            particles.device<rigidBody::forceField>(),
            particles.device<rigidBody::torqueField>(),
            particles.device<rigidBody::inverseMassField>(),
            particles.device<rigidBody::inertiaTensorField>(),
            particles.device<rigidBody::inverseInertiaTensorField>()};
}

ConstLSParticleDeviceView makeLSParticleDeviceView(const LSParticleContainer& particles) noexcept
{
    return {particles.device<rigidBody::positionField>(),
            particles.device<rigidBody::orientationField>(),
            particles.device<rigidBody::velocityField>(),
            particles.device<rigidBody::angularVelocityField>(),
            particles.device<rigidBody::forceField>(),
            particles.device<rigidBody::torqueField>(),
            particles.device<rigidBody::inverseMassField>(),
            particles.device<rigidBody::inertiaTensorField>(),
            particles.device<rigidBody::inverseInertiaTensorField>()};
}

__device__ bool spatialGridRange(int3& begin, int3& end, const math::Vec3& position, math::Real queryRadius, const spatialGridContainer::const_device_type& spatialGrid) noexcept
{
    const spatialGridLevel& level = spatialGrid.levels_[0];
    const math::Vec3 local = math::hadamardProduct(position - spatialGrid.minimumBoundary_, level.inverseCellSize_);
    const math::Vec3 extent = queryRadius * level.inverseCellSize_;
    begin = make_int3(static_cast<int>(math::detail::floor(local.x - extent.x)), static_cast<int>(math::detail::floor(local.y - extent.y)), static_cast<int>(math::detail::floor(local.z - extent.z)));
    end = make_int3(static_cast<int>(math::detail::floor(local.x + extent.x)), static_cast<int>(math::detail::floor(local.y + extent.y)), static_cast<int>(math::detail::floor(local.z + extent.z)));
    begin.x = begin.x < 0 ? 0 : (begin.x < level.size_.x ? begin.x : level.size_.x - 1);
    begin.y = begin.y < 0 ? 0 : (begin.y < level.size_.y ? begin.y : level.size_.y - 1);
    begin.z = begin.z < 0 ? 0 : (begin.z < level.size_.z ? begin.z : level.size_.z - 1);
    end.x = end.x < 0 ? 0 : (end.x < level.size_.x ? end.x : level.size_.x - 1);
    end.y = end.y < 0 ? 0 : (end.y < level.size_.y ? end.y : level.size_.y - 1);
    end.z = end.z < 0 ? 0 : (end.z < level.size_.z ? end.z : level.size_.z - 1);
    return true;
}

class uniformSpatialGridSearch
{
public:
    __device__ uniformSpatialGridSearch(const spatialGridContainer::const_device_type& spatialGrid, const math::Vec3& position, math::Real queryRadius) noexcept : spatialGrid_(spatialGrid)
    {
        cellValid_ = spatialGridRange(cellBegin_, cellEnd_, position, queryRadius, spatialGrid_);
        cell_ = cellBegin_;
    }

    __device__ bool next(int& particleIndex) noexcept
    {
        for (;;)
        {
            if (sortedIndex_ < sortedEnd_)
            {
                particleIndex = spatialGrid_.particleIndex_[sortedIndex_++];
                return true;
            }
            if (!loadNextCell())
            {
                return false;
            }
        }
    }

private:
    __device__ bool loadNextCell() noexcept
    {
        const spatialGridLevel& level = spatialGrid_.levels_[0];
        while (cellValid_)
        {
            const int hash = math::linearIndex(cell_.x, cell_.y, cell_.z, level.size_.x, level.size_.y);
            advanceCell();
            const int cellIndex = level.occupiedCellBegin_ + hash;
            const int sortedIndex = spatialGrid_.cellSortedParticleIndexBegin_[cellIndex];
            const int sortedEnd = spatialGrid_.cellSortedParticleIndexEnd_[cellIndex];
            if (sortedIndex >= 0 && sortedEnd > sortedIndex)
            {
                sortedIndex_ = sortedIndex;
                sortedEnd_ = sortedEnd;
                return true;
            }
        }
        return false;
    }

    __device__ void advanceCell() noexcept
    {
        if (cell_.x < cellEnd_.x)
        {
            ++cell_.x;
            return;
        }
        cell_.x = cellBegin_.x;
        if (cell_.y < cellEnd_.y)
        {
            ++cell_.y;
            return;
        }
        cell_.y = cellBegin_.y;
        if (cell_.z < cellEnd_.z)
        {
            ++cell_.z;
            return;
        }
        cellValid_ = false;
    }

    const spatialGridContainer::const_device_type& spatialGrid_;
    int3 cellBegin_{0, 0, 0};
    int3 cellEnd_{0, 0, 0};
    int3 cell_{0, 0, 0};
    int sortedIndex_{0};
    int sortedEnd_{0};
    bool cellValid_{false};
};

__global__ void initializeSPHParticlesKernel(SPHParticleDeviceView particles, math::Real referenceDensity, math::Real soundSpeed, int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }
    particles.force_[particleIndex] = math::Vec3::zero();
    particles.priorForce_[particleIndex] = math::Vec3::zero();
    particles.densityRate_[particleIndex] = 0.0;
    particles.pressure_[particleIndex] = execution::computePressureFromDensity(particles.density_[particleIndex], referenceDensity, soundSpeed);
}

__global__ void updateVirtualParticlePositionAndNormalKernel(VirtualParticleDeviceView virtualParticles, ConstLSParticleDeviceView LSParticles, int virtualParticleCount)
{
    const int virtualParticleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (virtualParticleIndex >= virtualParticleCount)
    {
        return;
    }
    const int ownerIndex = virtualParticles.ownerLSParticleIndex_[virtualParticleIndex];
    execution::calculateVirtualParticlePositionAndNormal(virtualParticles.position_[virtualParticleIndex],
                                                         virtualParticles.normal_[virtualParticleIndex],
                                                         virtualParticles.localPosition_[virtualParticleIndex],
                                                         virtualParticles.localNormal_[virtualParticleIndex],
                                                         LSParticles.position_[ownerIndex],
                                                         LSParticles.orientation_[ownerIndex]);
}

__global__ void accumulateVirtualParticleKinematicsKernel(VirtualParticleDeviceView virtualParticles, ConstLSParticleDeviceView LSParticles, math::Vec3 gravity, math::Real timeStep, int virtualParticleCount)
{
    const int virtualParticleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (virtualParticleIndex >= virtualParticleCount)
    {
        return;
    }
    const int ownerIndex = virtualParticles.ownerLSParticleIndex_[virtualParticleIndex];
    if (LSParticles.inverseMass_[ownerIndex] <= 0.0 && timeStep > 0.0)
    {
        const math::Vec3 velocity = LSParticles.velocity_[ownerIndex] + math::cross(LSParticles.angularVelocity_[ownerIndex], math::rotateUnit(LSParticles.orientation_[ownerIndex], virtualParticles.localPosition_[virtualParticleIndex]));
        execution::accumulatePrescribedSPHKinematics(virtualParticles.velocity_[virtualParticleIndex], virtualParticles.acceleration_[virtualParticleIndex],
                                                    virtualParticles.previousKinematicsVelocity_[virtualParticleIndex], velocity, timeStep);
        return;
    }
    execution::accumulateVirtualParticleVelocityAndAcceleration(virtualParticles.velocity_[virtualParticleIndex],
                                                                virtualParticles.acceleration_[virtualParticleIndex],
                                                                virtualParticles.localPosition_[virtualParticleIndex],
                                                                LSParticles.orientation_[ownerIndex],
                                                                LSParticles.velocity_[ownerIndex],
                                                                LSParticles.angularVelocity_[ownerIndex],
                                                                LSParticles.force_[ownerIndex],
                                                                LSParticles.torque_[ownerIndex],
                                                                LSParticles.inverseMass_[ownerIndex],
                                                                LSParticles.inertiaTensor_[ownerIndex],
                                                                LSParticles.inverseInertiaTensor_[ownerIndex],
                                                                gravity);
}

__global__ void resetVirtualParticleKinematicsReferenceKernel(VirtualParticleDeviceView virtualParticles, ConstLSParticleDeviceView LSParticles, int virtualParticleCount)
{
    const int index = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= virtualParticleCount)
        return;
    const int owner = virtualParticles.ownerLSParticleIndex_[index];
    virtualParticles.previousKinematicsVelocity_[index] = LSParticles.velocity_[owner] + math::cross(LSParticles.angularVelocity_[owner], math::rotateUnit(LSParticles.orientation_[owner], virtualParticles.localPosition_[index]));
}

__global__ void averageVirtualParticleKinematicsKernel(VirtualParticleDeviceView virtualParticles, ConstLSParticleDeviceView LSParticles, math::Real inverseSampleCount, int virtualParticleCount)
{
    const int virtualParticleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (virtualParticleIndex >= virtualParticleCount)
    {
        return;
    }
    const int ownerIndex = virtualParticles.ownerLSParticleIndex_[virtualParticleIndex];
    execution::calculateVirtualParticlePositionAndNormal(virtualParticles.position_[virtualParticleIndex],
                                                         virtualParticles.normal_[virtualParticleIndex],
                                                         virtualParticles.localPosition_[virtualParticleIndex],
                                                         virtualParticles.localNormal_[virtualParticleIndex],
                                                         LSParticles.position_[ownerIndex],
                                                         LSParticles.orientation_[ownerIndex]);
    virtualParticles.velocity_[virtualParticleIndex] *= inverseSampleCount;
    virtualParticles.acceleration_[virtualParticleIndex] *= inverseSampleCount;
}

__global__ void addVirtualParticleForceAndTorqueKernel(LSParticleDeviceView LSParticles, ConstVirtualParticleDeviceView virtualParticles, int virtualParticleCount)
{
    const int virtualParticleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (virtualParticleIndex >= virtualParticleCount)
    {
        return;
    }
    const int ownerIndex = virtualParticles.ownerLSParticleIndex_[virtualParticleIndex];
    // Both backends hold world force AND world torque between fluid updates.
    detail::atomicAddVec3(LSParticles.force_, ownerIndex, virtualParticles.force_[virtualParticleIndex]);
    detail::atomicAddVec3(LSParticles.torque_, ownerIndex, virtualParticles.couplingTorque_[virtualParticleIndex]);
}

__global__ void refreshVirtualParticleCouplingLoadsKernel(VirtualParticleDeviceView virtualParticles, ConstLSParticleDeviceView LSParticles, int virtualParticleCount)
{
    const int index = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= virtualParticleCount)
        return;
    const int owner = virtualParticles.ownerLSParticleIndex_[index];
    const math::Vec3 torque = math::cross(math::rotateUnit(LSParticles.orientation_[owner], virtualParticles.localPosition_[index]), virtualParticles.force_[index]);
    virtualParticles.couplingForceIncrement_[index] = virtualParticles.force_[index] - virtualParticles.couplingForceIncrement_[index];
    virtualParticles.couplingTorqueIncrement_[index] = torque - virtualParticles.couplingTorque_[index];
    virtualParticles.couplingTorque_[index] = torque;
}

__global__ void applyVirtualParticleImpulseCorrectionKernel(LSParticleDeviceView LSParticles, ConstVirtualParticleDeviceView virtualParticles, math::Real correctionTime, int virtualParticleCount)
{
    const int index = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (index >= virtualParticleCount)
        return;
    const int owner = virtualParticles.ownerLSParticleIndex_[index];
    math::Vec3 velocityCorrection;
    math::Vec3 angularVelocityCorrection;
    execution::calculateSPHCouplingVelocityCorrection(velocityCorrection, angularVelocityCorrection,
                                                     virtualParticles.couplingForceIncrement_[index], virtualParticles.couplingTorqueIncrement_[index],
                                                     LSParticles.orientation_[owner], LSParticles.inverseMass_[owner], LSParticles.inverseInertiaTensor_[owner], correctionTime);
    detail::atomicAddVec3(LSParticles.velocity_, owner, velocityCorrection);
    detail::atomicAddVec3(LSParticles.angularVelocity_, owner, angularVelocityCorrection);
}

__global__ void updateSPHFreeSurfaceKernel(SPHParticleDeviceView particles,
                                           ConstVirtualParticleDeviceView virtualParticles,
                                           spatialGridContainer::const_device_type particleSpatialGrid,
                                           spatialGridContainer::const_device_type virtualParticleSpatialGrid,
                                           math::Real smoothingLength,
                                           math::Real neighborSearchRadius,
                                           math::Real particleMass,
                                           int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }
    const math::Vec3 position = particles.position_[particleIndex];
    math::Real positionDivergence = 0.0;

    uniformSpatialGridSearch particleSearch(particleSpatialGrid, position, neighborSearchRadius);
    int neighborIndex;
    while (particleSearch.next(neighborIndex))
    {
        if (neighborIndex == particleIndex)
        {
            continue;
        }
        const math::Real neighborVolume = particleMass / execution::nonZeroDensity(particles.density_[neighborIndex]);
        positionDivergence += execution::freeSurfacePositionDivergenceContribution(position - particles.position_[neighborIndex], neighborVolume, smoothingLength);
    }

    uniformSpatialGridSearch wallSearch(virtualParticleSpatialGrid, position, neighborSearchRadius);
    while (wallSearch.next(neighborIndex))
    {
        positionDivergence += execution::freeSurfacePositionDivergenceContribution(position - virtualParticles.position_[neighborIndex], virtualParticles.volume_[neighborIndex], smoothingLength);
    }
    particles.freeSurface_[particleIndex] = execution::isFreeSurface(positionDivergence) ? 1 : 0;
}

__global__ void reinitializeSPHDensityKernel(SPHParticleDeviceView particles,
                                             ConstVirtualParticleDeviceView virtualParticles,
                                             spatialGridContainer::const_device_type particleSpatialGrid,
                                             spatialGridContainer::const_device_type virtualParticleSpatialGrid,
                                             math::Real smoothingLength,
                                             math::Real neighborSearchRadius,
                                             math::Real particleMass,
                                             math::Real referenceDensity,
                                             math::Real latticeKernelSum3D,
                                             math::Real soundSpeed,
                                             int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }
    if (particles.constrained_[particleIndex] != 0)
    {
        return;
    }
    const math::Vec3 position = particles.position_[particleIndex];
    const math::Real supportRadius = 2.0 * smoothingLength;
    const math::Real referenceVolume = particleMass / referenceDensity;
    math::Real kernelSum = execution::wendlandKernel3D(0.0, smoothingLength);
    uniformSpatialGridSearch particleSearch(particleSpatialGrid, position, neighborSearchRadius);
    int neighborIndex;
    while (particleSearch.next(neighborIndex))
    {
        if (neighborIndex == particleIndex)
        {
            continue;
        }
        const math::Real neighborDistance = math::norm(position - particles.position_[neighborIndex]);
        if (neighborDistance >= supportRadius)
        {
            continue;
        }
        kernelSum += execution::wendlandKernel3D(neighborDistance, smoothingLength);
    }

    uniformSpatialGridSearch wallSearch(virtualParticleSpatialGrid, position, neighborSearchRadius);
    while (wallSearch.next(neighborIndex))
    {
        kernelSum += virtualParticles.volume_[neighborIndex] / referenceVolume * execution::wendlandKernel3D(math::norm(position - virtualParticles.position_[neighborIndex]), smoothingLength);
    }

    const math::Real density = execution::reinitializedDensity(kernelSum, latticeKernelSum3D, referenceDensity);
    particles.density_[particleIndex] = density;
    particles.pressure_[particleIndex] = execution::computePressureFromDensity(density, referenceDensity, soundSpeed);
}

__global__ void computeSPHDensityRateKernel(SPHParticleDeviceView particles,
                                            ConstVirtualParticleDeviceView virtualParticles,
                                            spatialGridContainer::const_device_type particleSpatialGrid,
                                            spatialGridContainer::const_device_type virtualParticleSpatialGrid,
                                            math::Real smoothingLength,
                                            math::Real neighborSearchRadius,
                                            math::Real particleMass,
                                            math::Real referenceDensity,
                                            math::Real soundSpeed,
                                            math::Vec3 gravity,
                                            int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }
    if (particles.constrained_[particleIndex] != 0)
    {
        if (math::isFinite(particles.densityRate_[particleIndex]))
            particles.densityRate_[particleIndex] = 0.0;
        return;
    }

    const math::Vec3 position = particles.position_[particleIndex];
    const math::Vec3 velocity = particles.velocity_[particleIndex];
    const math::Real density = particles.density_[particleIndex];
    const math::Real pressure = particles.pressure_[particleIndex];
    math::Real densityRate = 0.0;

    uniformSpatialGridSearch particleSearch(particleSpatialGrid, position, neighborSearchRadius);
    int neighborIndex;
    while (particleSearch.next(neighborIndex))
    {
        if (neighborIndex == particleIndex)
        {
            continue;
        }
        densityRate += execution::calculateFluidDensityRate(position,
                                                            velocity,
                                                            density,
                                                            pressure,
                                                            particles.position_[neighborIndex],
                                                            particles.velocity_[neighborIndex],
                                                            particles.density_[neighborIndex],
                                                            particles.pressure_[neighborIndex],
                                                            particleMass,
                                                            referenceDensity,
                                                            smoothingLength,
                                                            soundSpeed);
    }

    uniformSpatialGridSearch wallSearch(virtualParticleSpatialGrid, position, neighborSearchRadius);
    while (wallSearch.next(neighborIndex))
    {
        densityRate += execution::calculateWallDensityRate(position,
                                                           velocity,
                                                           density,
                                                           pressure,
                                                           virtualParticles.position_[neighborIndex],
                                                           virtualParticles.normal_[neighborIndex],
                                                           virtualParticles.velocity_[neighborIndex],
                                                           virtualParticles.acceleration_[neighborIndex],
                                                           virtualParticles.volume_[neighborIndex],
                                                           referenceDensity,
                                                           smoothingLength,
                                                           soundSpeed,
                                                           gravity);
    }
    // Retain a failure across both split stages until the common state diagnostic observes it.
    if (math::isFinite(particles.densityRate_[particleIndex]))
        particles.densityRate_[particleIndex] = densityRate;
}

__global__ void integrateSPHDensityKernel(SPHParticleDeviceView particles, math::Real referenceDensity, math::Real soundSpeed, math::Real timeStep, int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }
    if (particles.constrained_[particleIndex] != 0)
    {
        return;
    }
    execution::integrateDensityAndPressure(particles.density_[particleIndex], particles.pressure_[particleIndex], particles.densityRate_[particleIndex], referenceDensity, soundSpeed, timeStep);
}

__global__ void updateSPHPriorForceAndBoundaryForceKernel(SPHParticleDeviceView particles,
                                                          VirtualParticleDeviceView virtualParticles,
                                                          spatialGridContainer::const_device_type particleSpatialGrid,
                                                          spatialGridContainer::const_device_type virtualParticleSpatialGrid,
                                                          math::Real smoothingLength,
                                                          math::Real neighborSearchRadius,
                                                          math::Real particleMass,
                                                          math::Real dynamicViscosity,
                                                          int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }

    const math::Vec3 position = particles.position_[particleIndex];
    const math::Vec3 velocity = particles.velocity_[particleIndex];
    const math::Real density = particles.density_[particleIndex];
    math::Vec3 acceleration = math::Vec3::zero();

    uniformSpatialGridSearch particleSearch(particleSpatialGrid, position, neighborSearchRadius);
    int neighborIndex;
    while (particleSearch.next(neighborIndex))
    {
        if (neighborIndex == particleIndex)
        {
            continue;
        }
        acceleration += execution::calculateFluidViscousAcceleration(position,
                                                                     velocity,
                                                                     density,
                                                                     particles.position_[neighborIndex],
                                                                     particles.velocity_[neighborIndex],
                                                                     particles.density_[neighborIndex],
                                                                     particleMass,
                                                                     dynamicViscosity,
                                                                     smoothingLength);
    }

    uniformSpatialGridSearch wallSearch(virtualParticleSpatialGrid, position, neighborSearchRadius);
    while (wallSearch.next(neighborIndex))
    {
        const math::Vec3 wallAcceleration = execution::calculateWallViscousAcceleration(position,
                                                                                        velocity,
                                                                                        density,
                                                                                        virtualParticles.position_[neighborIndex],
                                                                                        virtualParticles.velocity_[neighborIndex],
                                                                                        virtualParticles.volume_[neighborIndex],
                                                                                        dynamicViscosity,
                                                                                        smoothingLength);
        acceleration += wallAcceleration;
        detail::atomicAddVec3(virtualParticles.priorForce_, neighborIndex, execution::calculateWallForce(particleMass, wallAcceleration));
    }
    particles.priorForce_[particleIndex] = particleMass * acceleration;
}

__global__ void updateSPHPressureForceAndBoundaryForceKernel(SPHParticleDeviceView particles,
                                                             VirtualParticleDeviceView virtualParticles,
                                                             spatialGridContainer::const_device_type particleSpatialGrid,
                                                             spatialGridContainer::const_device_type virtualParticleSpatialGrid,
                                                             math::Real smoothingLength,
                                                             math::Real neighborSearchRadius,
                                                             math::Real particleMass,
                                                             math::Real referenceDensity,
                                                             math::Real soundSpeed,
                                                             math::Vec3 gravity,
                                                             int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex >= particleCount)
    {
        return;
    }

    const math::Vec3 position = particles.position_[particleIndex];
    const math::Vec3 velocity = particles.velocity_[particleIndex];
    const math::Real density = particles.density_[particleIndex];
    const math::Real pressure = particles.pressure_[particleIndex];
    math::Vec3 acceleration = math::Vec3::zero();

    uniformSpatialGridSearch particleSearch(particleSpatialGrid, position, neighborSearchRadius);
    int neighborIndex;
    while (particleSearch.next(neighborIndex))
    {
        if (neighborIndex == particleIndex)
        {
            continue;
        }
        acceleration += execution::calculateFluidPressureAcceleration(position,
                                                                      velocity,
                                                                      density,
                                                                      pressure,
                                                                      particles.position_[neighborIndex],
                                                                      particles.velocity_[neighborIndex],
                                                                      particles.density_[neighborIndex],
                                                                      particles.pressure_[neighborIndex],
                                                                      particleMass,
                                                                      referenceDensity,
                                                                      smoothingLength,
                                                                      soundSpeed);
    }

    uniformSpatialGridSearch wallSearch(virtualParticleSpatialGrid, position, neighborSearchRadius);
    while (wallSearch.next(neighborIndex))
    {
        const math::Vec3 wallAcceleration = execution::calculateWallPressureAcceleration(position,
                                                                                         velocity,
                                                                                         density,
                                                                                         pressure,
                                                                                         virtualParticles.position_[neighborIndex],
                                                                                         virtualParticles.normal_[neighborIndex],
                                                                                         virtualParticles.velocity_[neighborIndex],
                                                                                         virtualParticles.acceleration_[neighborIndex],
                                                                                         virtualParticles.volume_[neighborIndex],
                                                                                         referenceDensity,
                                                                                         smoothingLength,
                                                                                         soundSpeed,
                                                                                         gravity);
        acceleration += wallAcceleration;
        detail::atomicAddVec3(virtualParticles.force_, neighborIndex, execution::calculateWallForce(particleMass, wallAcceleration));
    }
    particles.force_[particleIndex] = particles.priorForce_[particleIndex] + particleMass * acceleration;
}

__global__ void integrateSPHVelocityKernel(SPHParticleDeviceView particles, math::Vec3 gravity, math::Real timeStep, int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex < particleCount)
    {
        execution::integrateVelocity(particles.velocity_[particleIndex], particles.force_[particleIndex], particles.inverseMass_[particleIndex], gravity, timeStep);
    }
}

__global__ void integrateSPHPositionKernel(SPHParticleDeviceView particles, math::Real timeStep, int particleCount)
{
    const int particleIndex = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    if (particleIndex < particleCount)
    {
        execution::integratePosition(particles.position_[particleIndex], particles.velocity_[particleIndex], timeStep);
    }
}

int blockCount(int count) noexcept { return (count + blockSize - 1) / blockSize; }

} // namespace

void launchInitializeSPHParticles(SPHParticleContainer& particles, math::Real referenceDensity, math::Real soundSpeed, cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    initializeSPHParticlesKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles), referenceDensity, soundSpeed, count);
    host_device_detail::checkCuda(cudaGetLastError(), "initializeSPHParticlesKernel launch");
}

void launchUpdateVirtualParticlePositionAndNormal(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    updateVirtualParticlePositionAndNormalKernel<<<blockCount(count), blockSize, 0, stream>>>(makeVirtualParticleDeviceView(virtualParticles), makeLSParticleDeviceView(LSParticles), count);
    host_device_detail::checkCuda(cudaGetLastError(), "updateVirtualParticlePositionAndNormalKernel launch");
}

void launchClearVirtualParticleKinematics(virtualParticleContainer& virtualParticles, cudaStream_t stream)
{
    const std::size_t byteCount = virtualParticles.deviceSize() * sizeof(math::Vec3);
    if (byteCount == 0)
    {
        return;
    }
    host_device_detail::checkCuda(cudaMemsetAsync(virtualParticles.device<pointMass::velocityField>(), 0, byteCount, stream), "cudaMemsetAsync virtual-particle velocity");
    host_device_detail::checkCuda(cudaMemsetAsync(virtualParticles.device<virtualParticle::accelerationField>(), 0, byteCount, stream), "cudaMemsetAsync virtual-particle acceleration");
}

void launchAccumulateVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, const math::Vec3& gravity, math::Real timeStep, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    accumulateVirtualParticleKinematicsKernel<<<blockCount(count), blockSize, 0, stream>>>(makeVirtualParticleDeviceView(virtualParticles), makeLSParticleDeviceView(LSParticles), gravity, timeStep, count);
    host_device_detail::checkCuda(cudaGetLastError(), "accumulateVirtualParticleKinematicsKernel launch");
}

void launchResetVirtualParticleKinematicsReference(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0)
        return;
    resetVirtualParticleKinematicsReferenceKernel<<<blockCount(count), blockSize, 0, stream>>>(makeVirtualParticleDeviceView(virtualParticles), makeLSParticleDeviceView(LSParticles), count);
    host_device_detail::checkCuda(cudaGetLastError(), "resetVirtualParticleKinematicsReferenceKernel launch");
}

void launchAverageVirtualParticleKinematics(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, math::Real inverseSampleCount, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    averageVirtualParticleKinematicsKernel<<<blockCount(count), blockSize, 0, stream>>>(makeVirtualParticleDeviceView(virtualParticles),
                                                                                        makeLSParticleDeviceView(LSParticles),
                                                                                        inverseSampleCount,
                                                                                        count);
    host_device_detail::checkCuda(cudaGetLastError(), "averageVirtualParticleKinematicsKernel launch");
}

void launchAddVirtualParticleForceAndTorque(LSParticleContainer& LSParticles, const virtualParticleContainer& virtualParticles, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    addVirtualParticleForceAndTorqueKernel<<<blockCount(count), blockSize, 0, stream>>>(makeLSParticleDeviceView(LSParticles), makeVirtualParticleDeviceView(virtualParticles), count);
    host_device_detail::checkCuda(cudaGetLastError(), "addVirtualParticleForceAndTorqueKernel launch");
}

void launchCacheVirtualParticleCouplingForce(virtualParticleContainer& virtualParticles, cudaStream_t stream)
{
    const auto count = virtualParticles.deviceSize();
    if (count == 0)
        return;
    host_device_detail::checkCuda(cudaMemcpyAsync(virtualParticles.device<virtualParticle::couplingForceIncrementField>(), virtualParticles.device<pointMass::forceField>(),
                                                count * sizeof(math::Vec3), cudaMemcpyDeviceToDevice, stream), "cache virtual-particle coupling force");
}

void launchRefreshVirtualParticleCouplingLoads(virtualParticleContainer& virtualParticles, const LSParticleContainer& LSParticles, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0)
        return;
    refreshVirtualParticleCouplingLoadsKernel<<<blockCount(count), blockSize, 0, stream>>>(makeVirtualParticleDeviceView(virtualParticles), makeLSParticleDeviceView(LSParticles), count);
    host_device_detail::checkCuda(cudaGetLastError(), "refreshVirtualParticleCouplingLoadsKernel launch");
}

void launchApplyVirtualParticleImpulseCorrection(LSParticleContainer& LSParticles, const virtualParticleContainer& virtualParticles, math::Real correctionTime, cudaStream_t stream)
{
    const int count = static_cast<int>(virtualParticles.deviceSize());
    if (count <= 0 || correctionTime <= 0.0)
        return;
    applyVirtualParticleImpulseCorrectionKernel<<<blockCount(count), blockSize, 0, stream>>>(makeLSParticleDeviceView(LSParticles), makeVirtualParticleDeviceView(virtualParticles), correctionTime, count);
    host_device_detail::checkCuda(cudaGetLastError(), "applyVirtualParticleImpulseCorrectionKernel launch");
}

void launchUpdateSPHFreeSurface(SPHParticleContainer& particles,
                                const virtualParticleContainer& virtualParticles,
                                const spatialGridContainer::const_device_type& particleSpatialGrid,
                                const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                math::Real smoothingLength,
                                math::Real neighborSearchRadius,
                                math::Real particleMass,
                                cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    updateSPHFreeSurfaceKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles),
                                                                            makeVirtualParticleDeviceView(virtualParticles),
                                                                            particleSpatialGrid,
                                                                            virtualParticleSpatialGrid,
                                                                            smoothingLength,
                                                                            neighborSearchRadius,
                                                                            particleMass,
                                                                            count);
    host_device_detail::checkCuda(cudaGetLastError(), "updateSPHFreeSurfaceKernel launch");
}

void launchReinitializeSPHDensity(SPHParticleContainer& particles,
                                  const virtualParticleContainer& virtualParticles,
                                  const spatialGridContainer::const_device_type& particleSpatialGrid,
                                  const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                  math::Real smoothingLength,
                                  math::Real neighborSearchRadius,
                                  math::Real particleMass,
                                  math::Real referenceDensity,
                                  math::Real latticeKernelSum3D,
                                  math::Real soundSpeed,
                                  cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    reinitializeSPHDensityKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles),
                                                                              makeVirtualParticleDeviceView(virtualParticles),
                                                                              particleSpatialGrid,
                                                                              virtualParticleSpatialGrid,
                                                                              smoothingLength,
                                                                              neighborSearchRadius,
                                                                              particleMass,
                                                                              referenceDensity,
                                                                              latticeKernelSum3D,
                                                                              soundSpeed,
                                                                              count);
    host_device_detail::checkCuda(cudaGetLastError(), "reinitializeSPHDensityKernel launch");
}

void launchUpdateSPHDensity(SPHParticleContainer& particles,
                            const virtualParticleContainer& virtualParticles,
                            const spatialGridContainer::const_device_type& particleSpatialGrid,
                            const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                            math::Real smoothingLength,
                            math::Real neighborSearchRadius,
                            math::Real particleMass,
                            math::Real referenceDensity,
                            math::Real soundSpeed,
                            math::Real timeStep,
                            const math::Vec3& gravity,
                            cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    computeSPHDensityRateKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles),
                                                                             makeVirtualParticleDeviceView(virtualParticles),
                                                                             particleSpatialGrid,
                                                                             virtualParticleSpatialGrid,
                                                                             smoothingLength,
                                                                             neighborSearchRadius,
                                                                             particleMass,
                                                                             referenceDensity,
                                                                             soundSpeed,
                                                                             gravity,
                                                                             count);
    host_device_detail::checkCuda(cudaGetLastError(), "computeSPHDensityRateKernel launch");
    integrateSPHDensityKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles), referenceDensity, soundSpeed, timeStep, count);
    host_device_detail::checkCuda(cudaGetLastError(), "integrateSPHDensityKernel launch");
}

void launchUpdateSPHPriorForceAndBoundaryForce(SPHParticleContainer& particles,
                                               virtualParticleContainer& virtualParticles,
                                               const spatialGridContainer::const_device_type& particleSpatialGrid,
                                               const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                               math::Real smoothingLength,
                                               math::Real neighborSearchRadius,
                                               math::Real particleMass,
                                               math::Real dynamicViscosity,
                                               cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    const std::size_t virtualForceBytes = virtualParticles.deviceSize() * sizeof(math::Vec3);
    if (virtualForceBytes > 0)
    {
        host_device_detail::checkCuda(cudaMemsetAsync(virtualParticles.device<virtualParticle::priorForceField>(), 0, virtualForceBytes, stream), "cudaMemsetAsync virtual-particle prior force");
    }
    updateSPHPriorForceAndBoundaryForceKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles),
                                                                                           makeVirtualParticleDeviceView(virtualParticles),
                                                                                           particleSpatialGrid,
                                                                                           virtualParticleSpatialGrid,
                                                                                           smoothingLength,
                                                                                           neighborSearchRadius,
                                                                                           particleMass,
                                                                                           dynamicViscosity,
                                                                                           count);
    host_device_detail::checkCuda(cudaGetLastError(), "updateSPHPriorForceAndBoundaryForceKernel launch");
}

void launchUpdateSPHPressureForceAndBoundaryForce(SPHParticleContainer& particles,
                                                  virtualParticleContainer& virtualParticles,
                                                  const spatialGridContainer::const_device_type& particleSpatialGrid,
                                                  const spatialGridContainer::const_device_type& virtualParticleSpatialGrid,
                                                  math::Real smoothingLength,
                                                  math::Real neighborSearchRadius,
                                                  math::Real particleMass,
                                                  math::Real referenceDensity,
                                                  math::Real soundSpeed,
                                                  const math::Vec3& gravity,
                                                  cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    const std::size_t virtualForceBytes = virtualParticles.deviceSize() * sizeof(math::Vec3);
    if (virtualForceBytes > 0)
    {
        host_device_detail::checkCuda(
            cudaMemcpyAsync(virtualParticles.device<pointMass::forceField>(), virtualParticles.device<virtualParticle::priorForceField>(), virtualForceBytes, cudaMemcpyDeviceToDevice, stream),
            "cudaMemcpyAsync virtual-particle prior force");
    }
    updateSPHPressureForceAndBoundaryForceKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles),
                                                                                              makeVirtualParticleDeviceView(virtualParticles),
                                                                                              particleSpatialGrid,
                                                                                              virtualParticleSpatialGrid,
                                                                                              smoothingLength,
                                                                                              neighborSearchRadius,
                                                                                              particleMass,
                                                                                              referenceDensity,
                                                                                              soundSpeed,
                                                                                              gravity,
                                                                                              count);
    host_device_detail::checkCuda(cudaGetLastError(), "updateSPHPressureForceAndBoundaryForceKernel launch");
}

void launchIntegrateSPHVelocity(SPHParticleContainer& particles, const math::Vec3& gravity, math::Real timeStep, cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    integrateSPHVelocityKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles), gravity, timeStep, count);
    host_device_detail::checkCuda(cudaGetLastError(), "integrateSPHVelocityKernel launch");
}

void launchIntegrateSPHPosition(SPHParticleContainer& particles, math::Real timeStep, cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return;
    }
    integrateSPHPositionKernel<<<blockCount(count), blockSize, 0, stream>>>(makeSPHParticleDeviceView(particles), timeStep, count);
    host_device_detail::checkCuda(cudaGetLastError(), "integrateSPHPositionKernel launch");
}

math::Real maximumSPHNeighborDisplacementSquared(const SPHParticleContainer& particles,
                                                 const virtualParticleContainer& boundaries,
                                                 const SPHNeighborPositionContainer& reference,
                                                 cudaStream_t stream)
{
    const int fluidCount = static_cast<int>(particles.deviceSize());
    const int count = fluidCount + static_cast<int>(boundaries.deviceSize());
    if (reference.deviceSize() != static_cast<std::size_t>(count))
        return std::numeric_limits<math::Real>::infinity();
    const readNeighborDisplacement read{particles.device<pointMass::positionField>(), boundaries.device<pointMass::positionField>(), reference.device<SPHNeighborPosition::positionField>(), fluidCount};
    const auto begin = thrust::make_counting_iterator(0);
    const math::Real result = thrust::transform_reduce(thrust::cuda::par.on(stream), begin, begin + count, read, math::Real{0.0}, maximumNeighborDisplacement{});
    host_device_detail::checkCuda(cudaGetLastError(), "SPH neighbor-displacement reduction");
    return result;
}

SPHStateStatistics calculateSPHStateStatistics(const SPHParticleContainer& particles, const math::Vec3& gravity, cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return {};
    }
    const auto begin = thrust::make_counting_iterator(0);
    const SPHStateStatistics initial{0.0, 0.0, std::numeric_limits<math::Real>::max(), 0.0, 0};
    const readSPHState read{particles.device<pointMass::velocityField>(),
                            particles.device<pointMass::forceField>(),
                            particles.device<pointMass::inverseMassField>(),
                            particles.device<SPHParticle::densityField>(),
                            particles.device<SPHParticle::pressureField>(),
                            particles.device<SPHParticle::densityRateField>(),
                            particles.device<pointMass::positionField>(),
                            gravity};
    const SPHStateStatistics result = thrust::transform_reduce(thrust::cuda::par.on(stream), begin, begin + count, read, initial, combineSPHStateStatistics{});
    host_device_detail::checkCuda(cudaGetLastError(), "SPH state-statistics reduction");
    host_device_detail::synchronize(stream);
    return result;
}

SPHKinematicsStatistics calculateSPHKinematicsStatistics(const virtualParticleContainer& particles, cudaStream_t stream)
{
    const int count = static_cast<int>(particles.deviceSize());
    if (count <= 0)
    {
        return {};
    }
    const auto begin = thrust::make_counting_iterator(0);
    const readSPHKinematics read{particles.device<pointMass::velocityField>(), particles.device<virtualParticle::accelerationField>()};
    const SPHKinematicsStatistics result = thrust::transform_reduce(thrust::cuda::par.on(stream), begin, begin + count, read, SPHKinematicsStatistics{}, combineSPHKinematicsStatistics{});
    host_device_detail::checkCuda(cudaGetLastError(), "SPH kinematics-statistics reduction");
    host_device_detail::synchronize(stream);
    return result;
}

} // namespace fundem::cuda
