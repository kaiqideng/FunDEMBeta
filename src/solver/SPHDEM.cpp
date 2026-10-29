#include "SPHDEM.h"

#include "execution/contactDetection.h"
#include "execution/sphFunctions.h"
#include "execution/sphJetFunctions.h"
#include "execution/sphSourceDiscretization.h"
#include "execution/cpu/particleFunctions.h"
#include "interaction/sphInteraction.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "data/HostAoSDeviceSoA.h"
#include "execution/cuda/integrationKernel.cuh"
#include "execution/cuda/sphInteractionKernel.cuh"
#include "execution/cuda/sphJetKernel.cuh"
#include "execution/cuda/spatialGridKernel.cuh"
#endif

namespace fundem
{
namespace
{

void validateSPHPeriodicDomain(const math::periodicDomain& value, math::Real smoothingLength)
{
    if (!value.valid())
        throw std::invalid_argument("The SPH periodic domain needs finite coordinates and positive enabled-axis lengths.");
    for (int axis = 0; axis < 3; ++axis)
        if (value.periodicAxes_[axis] && smoothingLength > 0.0 && !(4.0 * smoothingLength < value.length_[axis]))
            throw std::invalid_argument("The SPH kernel-support diameter must be smaller than each enabled periodic length (4 h < L).");
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
/** Only periodic axes use the exact cell bounds; other axes keep their search padding. */
template <class Particles>
void setSPHSpatialGrid(spatialGridContainer& grid,
                       const Particles& particles,
                       math::Real support,
                       math::Vec3 minimum,
                       math::Vec3 maximum,
                       const math::periodicDomain& periodicDomain,
                       cudaStream_t stream)
{
    minimum -= math::Vec3{support};
    maximum += math::Vec3{support};
    for (int axis = 0; axis < 3; ++axis)
        if (periodicDomain.periodicAxes_[axis])
        {
            minimum[axis] = periodicDomain.minimum_[axis];
            maximum[axis] = minimum[axis] + periodicDomain.length_[axis];
        }
    grid.setUniform(particles, support, minimum, maximum, stream);
}

/** Keeps an acoustic endpoint entirely on the solver stream, without a host round trip. */
template <class... Fields>
void copySPHDeviceState(hostAoSDeviceSoA<SPHParticle, deviceLayout<Fields...>>& destination, const SPHParticleContainer& source, deviceLayout<Fields...>, cudaStream_t stream)
{
    const auto count = source.deviceSize();
    destination.resizeDevice(count);
    if (count == 0)
        return;
    (host_device_detail::checkCuda(
         cudaMemcpyAsync(destination.device<Fields>(), source.device<Fields>(), count * sizeof(typename SPHParticleContainer::template device_value_type<Fields>), cudaMemcpyDeviceToDevice, stream),
         "cudaMemcpyAsync SPH acoustic endpoint"),
     ...);
}
#endif

void orthogonalBasis(const math::Vec3& unitAxis, math::Vec3& tangent1, math::Vec3& tangent2) noexcept
{
    tangent1 = std::abs(unitAxis.x) > std::abs(unitAxis.z) ? math::Vec3{-unitAxis.y, unitAxis.x, 0.0} : math::Vec3{0.0, -unitAxis.z, unitAxis.y};
    tangent1 = math::normalizedOrZero(tangent1);
    tangent2 = math::cross(unitAxis, tangent1);
}

/** Creates a static downstream cylinder, without an inlet constraint. */
std::vector<SPHParticle> makeSPHJetParticles(const math::Vec3& center,
                                             const math::Vec3& direction,
                                             const math::Vec3& velocity,
                                             math::Real spacing,
                                             const execution::SPHJetDiscretization& discretization)
{
    math::Vec3 tangent1;
    math::Vec3 tangent2;
    orthogonalBasis(direction, tangent1, tangent2);
    std::vector<SPHParticle> particles;
    particles.reserve(static_cast<std::size_t>(discretization.particleCount_));
    for (int axialIndex = 0; axialIndex < discretization.axialCount_; ++axialIndex)
    {
        const math::Vec3 offset = (math::Real(axialIndex) + 0.5) * spacing * direction;
        const math::Vec3 axialPosition = center + offset;
        for (int y = -discretization.radialIndex_; y <= discretization.radialIndex_; ++y)
        {
            for (int x = -discretization.radialIndex_; x <= discretization.radialIndex_; ++x)
            {
                if (discretization.contains(x, y))
                    particles.emplace_back(axialPosition + math::Real(x) * spacing * tangent1 + math::Real(y) * spacing * tangent2, velocity);
            }
        }
    }
    return particles;
}

void appendVirtualParticles(std::vector<virtualParticle>& virtualParticles,
                            const SPHParticleContainer& SPHParticles,
                            const LSParticleContainer& LSParticles,
                            math::Real samplingSpacing,
                            math::Real boundaryBand)
{
    const auto& SPHParticleHost = SPHParticles.host();
    if (SPHParticleHost.empty() || samplingSpacing <= 0.0)
    {
        return;
    }

    const math::Vec3 samplingOrigin = SPHParticleHost.front().position();
    const math::Real samplingVolume = samplingSpacing * samplingSpacing * samplingSpacing;
    const auto& LSParticleHost = LSParticles.host();
    for (int particleIndex = 0; particleIndex < static_cast<int>(LSParticleHost.size()); ++particleIndex)
    {
        const LSParticle& particle = LSParticleHost[particleIndex];
        const int3 size = particle.gridNodeSize();
        const math::Real inverseGridSpacing = particle.gridNodeInverseSpacing();
        if (inverseGridSpacing <= 0.0)
        {
            continue;
        }

        const math::Real gridSpacing = 1.0 / inverseGridSpacing;
        const math::Vec3 localMinimum = particle.gridNodeOrigin();
        const math::Vec3 localMaximum = localMinimum + gridSpacing * math::Vec3{math::Real(size.x - 1), math::Real(size.y - 1), math::Real(size.z - 1)};
        math::Vec3 worldMinimum = particle.position() + math::rotateUnit(particle.orientation(), localMinimum);
        math::Vec3 worldMaximum = worldMinimum;
        for (int z = 0; z < 2; ++z)
        {
            for (int y = 0; y < 2; ++y)
            {
                for (int x = 0; x < 2; ++x)
                {
                    const math::Vec3 localCorner{x == 0 ? localMinimum.x : localMaximum.x, y == 0 ? localMinimum.y : localMaximum.y, z == 0 ? localMinimum.z : localMaximum.z};
                    const math::Vec3 worldCorner = particle.position() + math::rotateUnit(particle.orientation(), localCorner);
                    worldMinimum = math::componentMin(worldMinimum, worldCorner);
                    worldMaximum = math::componentMax(worldMaximum, worldCorner);
                }
            }
        }

        const int minimumX = static_cast<int>(std::ceil((worldMinimum.x - samplingOrigin.x) / samplingSpacing));
        const int minimumY = static_cast<int>(std::ceil((worldMinimum.y - samplingOrigin.y) / samplingSpacing));
        const int minimumZ = static_cast<int>(std::ceil((worldMinimum.z - samplingOrigin.z) / samplingSpacing));
        const int maximumX = static_cast<int>(std::floor((worldMaximum.x - samplingOrigin.x) / samplingSpacing));
        const int maximumY = static_cast<int>(std::floor((worldMaximum.y - samplingOrigin.y) / samplingSpacing));
        const int maximumZ = static_cast<int>(std::floor((worldMaximum.z - samplingOrigin.z) / samplingSpacing));
        for (int z = minimumZ; z <= maximumZ; ++z)
        {
            for (int y = minimumY; y <= maximumY; ++y)
            {
                for (int x = minimumX; x <= maximumX; ++x)
                {
                    const math::Vec3 worldPosition = samplingOrigin + samplingSpacing * math::Vec3{math::Real(x), math::Real(y), math::Real(z)};
                    math::Real overlap = 0.0;
                    math::Vec3 worldNormal = math::Vec3::zero();
                    if (!execution::detectLevelSetContact(overlap,
                                                          worldNormal,
                                                          particle.gridNodes().data(),
                                                          worldPosition,
                                                          particle.position(),
                                                          particle.orientation(),
                                                          particle.gridNodeOrigin(),
                                                          inverseGridSpacing,
                                                          size,
                                                          0) ||
                        overlap <= math::defaultTolerance || overlap > boundaryBand)
                    {
                        continue;
                    }

                    const math::Vec3 localPosition = math::inverseRotateUnit(particle.orientation(), worldPosition - particle.position());
                    const math::Vec3 localNormal = math::inverseRotateUnit(particle.orientation(), worldNormal);
                    virtualParticle value;
                    if (value.setLocalValues(localPosition, localNormal, samplingVolume, particleIndex))
                    {
                        value.setPosition(worldPosition);
                        value.setNormal(worldNormal);
                        virtualParticles.push_back(std::move(value));
                    }
                }
            }
        }
    }
}

} // namespace

SPHDEM::SPHDEM(int device) : LSDEM(device) {}

SPHDEM::SPHDEM(executionMode mode, int device) : LSDEM(device) { setExecutionMode(mode); }

SPHDEM::~SPHDEM() = default;

int SPHDEM::addSPHParticle(const SPHParticle& value) { return addSPHParticles({value}); }

int SPHDEM::addSPHParticles(std::vector<SPHParticle> values)
{
    if (initialized())
    {
        throw std::logic_error("SPH particles cannot be added after solver initialization.");
    }
    if (!SPHProperties_.set_)
    {
        throw std::logic_error("Set the SPH properties before adding SPH particles.");
    }
    const Real particleMass = SPHProperties_.particleMass();
    for (const SPHParticle& value : values)
    {
        if (!value.isValid())
        {
            throw std::invalid_argument("Cannot add an invalid SPH particle.");
        }
    }

    const int firstParticleIndex = static_cast<int>(SPHParticles_.hostSize());
    auto& particleHost = SPHParticles_.host();
    particleHost.reserve(particleHost.size() + values.size());
    for (SPHParticle& value : values)
    {
        value.setMass(particleMass);
        value.setDensity(SPHProperties_.referenceDensity_);
        particleHost.push_back(std::move(value));
    }
    invalidateContinuation();
    return firstParticleIndex;
}

int SPHDEM::addSPHBlock(const Vec3& blockMinimum, const int3& particleCount, const Vec3& initialVelocity)
{
    if (!SPHProperties_.set_)
    {
        throw std::logic_error("Set the SPH properties before adding an SPH block.");
    }
    SPHParticle prototype{blockMinimum, initialVelocity};
    if (particleCount.x <= 0 || particleCount.y <= 0 || particleCount.z <= 0 || !prototype.isValid())
    {
        throw std::invalid_argument("Invalid SPH block definition.");
    }

    std::vector<SPHParticle> particles;
    particles.reserve(particleCount.x * particleCount.y * particleCount.z);
    for (int z = 0; z < particleCount.z; ++z)
    {
        for (int y = 0; y < particleCount.y; ++y)
        {
            for (int x = 0; x < particleCount.x; ++x)
            {
                const Vec3 position = blockMinimum + SPHProperties_.spacing_ * Vec3{Real(x) + 0.5, Real(y) + 0.5, Real(z) + 0.5};
                particles.emplace_back(position, initialVelocity);
            }
        }
    }
    return addSPHParticles(std::move(particles));
}

int SPHDEM::addSPHJet(const Vec3& inletCenter, const Vec3& velocity, Real radius, Real length)
{
    if (!SPHProperties_.set_)
    {
        throw std::logic_error("Set the SPH properties before adding an SPH jet.");
    }
    Vec3 axis = velocity;
    SPHParticle prototype{inletCenter, velocity};
    execution::SPHJetDiscretization discretization;
    if (!math::tryNormalize(axis) || !execution::trySPHJetDiscretization(SPHProperties_.spacing_, radius, length, discretization) || !prototype.isValid())
    {
        throw std::invalid_argument("Invalid SPH jet definition.");
    }

    return addSPHParticles(makeSPHJetParticles(inletCenter, axis, velocity, SPHProperties_.spacing_, discretization));
}

int SPHDEM::addSPHJet(const SPHJet& value)
{
    if (initialized())
    {
        throw std::logic_error("SPH jets cannot be added after solver initialization.");
    }
    if (!SPHProperties_.set_)
    {
        throw std::logic_error("Set the SPH properties before adding an SPH jet.");
    }
    if (!value.isValid())
    {
        throw std::invalid_argument("Invalid finite-duration SPH jet definition.");
    }

    const Real pipeLength = value.pipeLength();
    execution::SPHJetDiscretization discretization;
    if (!execution::trySPHJetDiscretization(SPHProperties_.spacing_, value.radius(), pipeLength, discretization))
    {
        throw std::invalid_argument("The SPH jet must contain at least one particle layer.");
    }

    SPHJetState::constraint source;
    source.value_ = value;
    source.startTime_ = time();
    source.axialCount_ = discretization.axialCount_;
    source.bufferLayerCount_ = static_cast<int>(std::min(Real(discretization.axialCount_), std::ceil(2.0 * SPHProperties_.smoothingLength_ / SPHProperties_.spacing_) + 1.0));
    const int layerCount = discretization.particleCount_ / discretization.axialCount_;
    source.offsets_.reserve(layerCount);
    Vec3 tangent1;
    Vec3 tangent2;
    orthogonalBasis(value.direction(), tangent1, tangent2);
    for (int y = -discretization.radialIndex_; y <= discretization.radialIndex_; ++y)
        for (int x = -discretization.radialIndex_; x <= discretization.radialIndex_; ++x)
            if (discretization.contains(x, y))
                source.offsets_.push_back(SPHProperties_.spacing_ * (Real(x) * tangent1 + Real(y) * tangent2));
    SPHJets_.constraints_.reserve(SPHJets_.constraints_.size() + 1);
    source.particleIndices_.reserve(layerCount * source.bufferLayerCount_);
    SPHParticles_.host().reserve(SPHParticles_.hostSize() + layerCount * source.bufferLayerCount_);
    const int firstParticleIndex = static_cast<int>(SPHParticles_.hostSize());
    for (int layer = 0; layer < source.bufferLayerCount_; ++layer)
        appendSPHJetLayer(source);
    SPHJets_.latestEndTime_ = std::max(SPHJets_.latestEndTime_, source.endTime());
    SPHJets_.constraints_.push_back(std::move(source));
    SPHJets_.completed_ = false;
    invalidateContinuation();
    return firstParticleIndex;
}

const std::vector<int>& SPHDEM::SPHJetParticleIndices(int firstParticleIndex) const
{
    const auto found =
        std::find_if(SPHJets_.constraints_.begin(), SPHJets_.constraints_.end(), [firstParticleIndex](const auto& source) { return source.particleIndices_.front() == firstParticleIndex; });
    if (found == SPHJets_.constraints_.end())
        throw std::out_of_range("The SPH jet handle does not identify a registered source.");
    return found->particleIndices_;
}

void SPHDEM::setSPHMaximumVelocity(Real value)
{
    if (initialized())
    {
        throw std::logic_error("The SPH maximum velocity cannot be changed while the solver is initialized.");
    }
    if (!math::isFinite(value) || value <= 0.0)
    {
        throw std::invalid_argument("The SPH maximum velocity must be finite and positive.");
    }
    if (value == SPHProperties_.maximumVelocity_ && SPHProperties_.soundSpeed_ == 10.0 * value)
    {
        return;
    }
    SPHProperties_.maximumVelocity_ = value;
    SPHProperties_.soundSpeed_ = 10.0 * value;
    invalidateContinuation();
}

void SPHDEM::setSPHProperties(Real spacing, Real smoothingLength, Real referenceDensity, Real dynamicViscosity)
{
    if (initialized() || SPHParticles_.hostSize() > 0)
    {
        throw std::logic_error("Set the SPH properties before adding particles or initializing the solver.");
    }
    if (!math::isFinite(spacing) || spacing <= 0.0 || !math::isFinite(smoothingLength) || smoothingLength <= 0.0 || !math::isFinite(referenceDensity) || referenceDensity <= 0.0 ||
        !math::isFinite(dynamicViscosity) || dynamicViscosity < 0.0)
    {
        throw std::invalid_argument("SPH spacing, smoothing length and reference density must be positive; dynamic viscosity must be non-negative.");
    }
    validateSPHPeriodicDomain(SPHPeriodicDomain_, smoothingLength);
    const Real smoothingLengthRatio = smoothingLength / spacing;
    if (smoothingLengthRatio < execution::minimumSPHSmoothingLengthToSpacingRatio || smoothingLengthRatio > execution::maximumSPHSmoothingLengthToSpacingRatio)
    {
        throw std::invalid_argument("The SPH smoothing length must be between one and two particle spacings for the current three-dimensional Wendland discretization.");
    }
    const Real latticeKernelSum3D = execution::calculateLatticeKernelSum3D(spacing, smoothingLength);
    if (!math::isFinite(latticeKernelSum3D) || latticeKernelSum3D <= math::defaultTolerance)
    {
        throw std::invalid_argument("The SPH reference kernel sum is invalid.");
    }
    SPHProperties_.spacing_ = spacing;
    SPHProperties_.smoothingLength_ = smoothingLength;
    SPHProperties_.referenceDensity_ = referenceDensity;
    SPHProperties_.dynamicViscosity_ = dynamicViscosity;
    SPHProperties_.latticeKernelSum3D_ = latticeKernelSum3D;
    SPHProperties_.viscousTimeScale_ = dynamicViscosity > 0.0 ? referenceDensity * smoothingLength * smoothingLength / dynamicViscosity : 0.0;
    SPHProperties_.set_ = true;
    invalidateContinuation();
}

void SPHDEM::setSPHPeriodicDomain(const math::periodicDomain& value)
{
    validateSPHPeriodicDomain(value, SPHProperties_.smoothingLength_);
    bool changed = value.minimum_ != SPHPeriodicDomain_.minimum_ || value.length_ != SPHPeriodicDomain_.length_;
    for (int axis = 0; axis < 3; ++axis)
        changed = changed || value.periodicAxes_[axis] != SPHPeriodicDomain_.periodicAxes_[axis];
    if (!changed)
        return;
    if (initialized())
        throw std::logic_error("The SPH periodic domain cannot change while the solver is initialized.");
    SPHPeriodicDomain_ = value;
    invalidateContinuation();
}

void SPHDEM::setSPHSoundSpeed(Real value)
{
    if (initialized())
    {
        throw std::logic_error("The SPH sound speed cannot be changed while the solver is initialized.");
    }
    if (!math::isFinite(value) || value <= 0.0)
    {
        throw std::invalid_argument("The SPH sound speed must be finite and positive.");
    }
    if (value < execution::minimumSPHSoundSpeedToVelocityRatio * SPHProperties_.maximumVelocity_)
    {
        throw std::invalid_argument("The SPH sound speed must be at least ten times the configured maximum fluid velocity.");
    }
    if (value == SPHProperties_.soundSpeed_)
    {
        return;
    }
    SPHProperties_.soundSpeed_ = value;
    invalidateContinuation();
}

void SPHDEM::invalidateDeferredState() noexcept
{
    invalidateSPHNeighborhood();
    outputSnapshot_ = {};
    outputObservationSnapshot_ = false;
    SPHStepState_ = {};
}

void SPHDEM::resumeContinuation() noexcept { restoreSPHState(); }

void SPHDEM::setSPHParticleVTUFields(std::vector<SPHParticleVTUField> fields) { SPHParticleVTUFields_ = std::move(fields); }

void SPHDEM::validateExecutionMode(executionMode value) const
{
    switch (value)
    {
    case executionMode::CPU:
    case executionMode::GPU:
    case executionMode::Hybrid:
        return;
    }
    throw std::invalid_argument("Invalid execution mode.");
}

void SPHDEM::executionModeChanged(executionMode value)
{
    SPHParticles_.setUseDevice(value != executionMode::CPU);
    mutableLSParticles().setUseDevice(value == executionMode::GPU);
}

void SPHDEM::validateSystemConfiguration() const
{
#if !defined(FUNDEM_HAS_CUDA) || !FUNDEM_HAS_CUDA
    if (mode() != executionMode::CPU)
    {
        throw std::runtime_error("The selected SPHDEM execution mode requires a CUDA-enabled FunDEM build.");
    }
#endif
    if (!math::isFinite(SPHProperties_.maximumVelocity_) || SPHProperties_.maximumVelocity_ <= 0.0 || !math::isFinite(SPHProperties_.soundSpeed_) || SPHProperties_.soundSpeed_ <= 0.0)
    {
        throw std::invalid_argument("Invalid SPH maximum velocity or sound speed.");
    }
    if (SPHProperties_.soundSpeed_ < execution::minimumSPHSoundSpeedToVelocityRatio * SPHProperties_.maximumVelocity_)
    {
        throw std::invalid_argument("The SPH sound speed must be at least ten times the configured maximum fluid velocity.");
    }
    if (SPHParticles_.hostSize() > 0 && !SPHProperties_.set_)
    {
        throw std::logic_error("SPH properties have not been set.");
    }
}

void SPHDEM::initialize(cpu::mode)
{
    validateContactMaterials();
    initializeSPH(cpu::mode{});
    assembleForceAndTorque(0.0, cpu::mode{});
    virtualParticleCoupling_.resetKinematicsReference(virtualParticles_, LSParticles());
}

void SPHDEM::initialize(gpu::mode, cudaStream_t stream)
{
    initializeLSParticleDevice(stream);
    initializeSPHDevice(stream);
    assembleForceAndTorque(0.0, gpu::mode{}, stream);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    cuda::launchResetVirtualParticleKinematicsReference(virtualParticles_, LSParticles(), stream);
#endif
}

void SPHDEM::initialize(hybrid::mode, cudaStream_t stream)
{
    validateContactMaterials();
    initializeSPHDevice(stream);
    assembleForceAndTorque(0.0, cpu::mode{});
    virtualParticleCoupling_.resetKinematicsReference(virtualParticles_, LSParticles());
}

void SPHDEM::advanceCPU(Real timeStep) { advanceHostParticles(timeStep); }

void SPHDEM::advanceHybrid(Real timeStep) { advanceHostParticles(timeStep); }

void SPHDEM::advanceHostParticles(Real timeStep)
{
    restoreSPHState();
    beforeLSParticleDrift(cpu::mode{});
    beginSPHInterval(gpuStream());
    // Sample the complete left endpoint with the NEW interval load installed.
    if (!SPHBoundaryStatic_)
        virtualParticleCoupling_.accumulateKinematics(virtualParticles_, LSParticles(), gravity(), timeStep);
    const Real halfTimeStep = 0.5 * timeStep;
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cpu::integratePosition(mutableLSParticles(), timeStep);
    assembleForceAndTorque(timeStep, cpu::mode{});
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    ++SPHStepState_.completedDEMSteps_;
}

void SPHDEM::advanceGPU(Real timeStep)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const cudaStream_t stream = gpuStream();
    const Real halfTimeStep = 0.5 * timeStep;
    restoreSPHState();
    beforeLSParticleDrift(gpu::mode{}, stream);
    beginSPHInterval(stream);
    if (!SPHBoundaryStatic_)
        cuda::launchAccumulateVirtualParticleKinematics(virtualParticles_, LSParticles(), gravity(), timeStep, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(mutableLSParticles(), gravity(), halfTimeStep, stream);
    cuda::launchPositionAndOrientationIntegration(mutableLSParticles(), timeStep, stream);
    assembleForceAndTorque(timeStep, gpu::mode{}, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(mutableLSParticles(), gravity(), halfTimeStep, stream);
    ++SPHStepState_.completedDEMSteps_;
#else
    (void)timeStep;
#endif
}

void SPHDEM::assembleForceAndTorque(Real historyTimeStep, cpu::mode)
{
    cpu::clearForceAndTorque(mutableLSParticles());
    calculateLSParticleForceAndTorque(historyTimeStep, cpu::mode{});
    virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles());
    addLSParticleExternalForceAndTorque(mutableLSParticles().host());
}

void SPHDEM::assembleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const int particleCount = static_cast<int>(LSParticles().deviceSize());
    if (particleCount > 0)
    {
        host_device_detail::checkCuda(cudaMemsetAsync(mutableLSParticles().device<rigidBody::forceField>(), 0, particleCount * sizeof(Vec3), stream), "cudaMemsetAsync LSParticle force");
        host_device_detail::checkCuda(cudaMemsetAsync(mutableLSParticles().device<rigidBody::torqueField>(), 0, particleCount * sizeof(Vec3), stream), "cudaMemsetAsync LSParticle torque");
    }
    calculateLSParticleForceAndTorque(historyTimeStep, gpu::mode{}, stream);
    cuda::launchAddVirtualParticleForceAndTorque(mutableLSParticles(), virtualParticles_, 1.0, stream);
    addLSParticleExternalForceAndTorque(deviceFields(mutableLSParticles()), stream);
#else
    (void)historyTimeStep;
    (void)stream;
#endif
}

void SPHDEM::initializeSPHState()
{
    const int SPHParticleCount = static_cast<int>(SPHParticles_.hostSize());
    const int LSParticleCount = static_cast<int>(LSParticles().hostSize());
    SPHBoundaryStatic_ =
        !SPHBoundaryMotionExpected_ && std::all_of(LSParticles().host().begin(),
                                                   LSParticles().host().end(),
                                                   [](const LSParticle& value) { return value.inverseMass() == 0.0 && value.velocity() == Vec3::zero() && value.angularVelocity() == Vec3::zero(); });
    outputSnapshot_ = {};
    outputObservationSnapshot_ = false;
    virtualParticles_.resetDevice();
    virtualParticles_.host().clear();
    SPHStepState_ = {};
    SPHStepState_.neighborSearchRadius_ = 2.0 * SPHProperties_.smoothingLength_;
    SPHJets_.completed_ = SPHJets_.completedAt(time());
    auto& particleHost = SPHParticles_.host();
    for (SPHParticle& value : particleHost)
    {
        value.setConstrained(false);
    }
    applySPHJetVelocity(time(), cpu::mode{});
    for (SPHParticle& value : particleHost)
    {
        value.setPressure(execution::computePressureFromDensity(value.density(), SPHProperties_.referenceDensity_, SPHProperties_.soundSpeed_));
    }

    if (SPHParticleCount > 0)
    {
        for (const LSParticle& boundary : LSParticles().host())
        {
            const Real boundarySpacing = 1.0 / boundary.gridNodeInverseSpacing();
            if (boundarySpacing > SPHProperties_.smoothingLength_ && !math::nearlyEqual(boundarySpacing, SPHProperties_.smoothingLength_))
            {
                throw std::invalid_argument("An SPH boundary level-set grid is coarser than the smoothing length. Refine the boundary geometry grid.");
            }
        }
    }

    std::vector<virtualParticle> newVirtualParticles;
    appendVirtualParticles(newVirtualParticles, SPHParticles_, LSParticles(), SPHProperties_.spacing_, 2.0 * SPHProperties_.smoothingLength_);
    if (SPHParticleCount > 0 && LSParticleCount > 0)
    {
        std::vector<int> virtualParticleCountByOwner(static_cast<std::size_t>(LSParticleCount), 0);
        for (const virtualParticle& value : newVirtualParticles)
        {
            ++virtualParticleCountByOwner[static_cast<std::size_t>(value.ownerLSParticleIndex())];
        }
        if (std::any_of(virtualParticleCountByOwner.begin(), virtualParticleCountByOwner.end(), [](int count) { return count == 0; }))
        {
            throw std::invalid_argument("Every SPH boundary LSParticle must generate at least one virtual particle inside the two-smoothing-length boundary band.");
        }
    }
    virtualParticles_.host() = std::move(newVirtualParticles);
    auto& openings = SPHJetOpenings_.host();
    openings.clear();
    for (const auto& source : SPHJets_.constraints_)
        openings.push_back({{source.value_.outletCenter(), source.value_.direction(), source.value_.radius(), source.bufferLayerCount_ * SPHProperties_.spacing_, SPHPeriodicDomain_}});
    if (mode() == executionMode::GPU)
    {
        virtualParticleCoupling_.reset();
    }
    else
    {
        virtualParticleCoupling_.initialize(virtualParticles_, LSParticles(), gravity());
    }
    updateSPHJetOpenings(cpu::mode{});
}

void SPHDEM::initializeSPH(cpu::mode)
{
    initializeSPHState();
    if (!SPHInteractions_)
    {
        SPHInteractions_ = std::make_unique<cpu::SPHInteraction>();
    }
    const Real supportRadius = SPHProperties_.smoothingLength_ > math::defaultTolerance ? 2.0 * SPHProperties_.smoothingLength_ : math::norm(maximumBoundary() - minimumBoundary());
    const Vec3 boundaryPadding{supportRadius, supportRadius, supportRadius};
    const auto parameters = SPHProperties_.interactionParameters(gravity(), SPHPeriodicDomain_);

    SPHStepState_.neighborSearchRadius_ = supportRadius;
    SPHInteractions_->initializeParticles(SPHParticles_, SPHProperties_.referenceDensity_, SPHProperties_.soundSpeed_);
    SPHInteractions_->buildNeighborhood(SPHParticles_, virtualParticles_, minimumBoundary() - boundaryPadding, maximumBoundary() + boundaryPadding, supportRadius, supportRadius, SPHPeriodicDomain_);
    SPHInteractions_->updatePriorForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->updatePressureForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
    updateSPHAcousticTimeStep(cpu::mode{});
    updateSPHAdvectionTimeStep();
    virtualParticleCoupling_.collectForceAndTorque(virtualParticles_, LSParticles());
}

void SPHDEM::beginSPHInterval(cudaStream_t stream)
{
    if (SPHStepState_.couplingActive())
        return;

    emitSPHJetParticles(stream);
    // Remove ONLY the old hydraulic component from the cached left-endpoint
    // load. Contact/Bond histories and external callbacks must not run twice.
    SPHStepState_.couplingDEMSteps_ = SPHJetInterval(SPHStepState_.acoustic_.interval_);
    const Real interval = SPHStepState_.couplingDEMSteps_ * timeStep();
    switch (mode())
    {
    case executionMode::CPU:
        advanceCouplingInterval(interval, cpu::mode{});
        virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles());
        break;
    case executionMode::GPU:
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        advanceCouplingInterval(interval, gpu::mode{}, stream);
        cuda::launchAddVirtualParticleForceAndTorque(mutableLSParticles(), virtualParticles_, 1.0, stream);
#endif
        break;
    case executionMode::Hybrid:
        advanceCouplingInterval(interval, hybrid::mode{}, stream);
        virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles());
        break;
    }
    SPHStepState_.completedDEMSteps_ = 0;
}

void SPHDEM::emitSPHJetParticles(cudaStream_t stream)
{
    if (SPHJets_.completed_)
        return;
    auto& particles = SPHParticles_.host();
    const auto oldSize = particles.size();
    for (auto& source : SPHJets_.constraints_)
    {
        while (source.firstInletLayer_ < source.layerCount() && (time() >= source.endTime() || time() >= source.releaseTime(source.firstInletLayer_, SPHProperties_.spacing_)))
            ++source.firstInletLayer_;
        while (source.layerCount() < source.axialCount_ && (source.nextEmissionTime(SPHProperties_.spacing_) - time()) / timeStep() <= math::defaultTolerance)
            appendSPHJetLayer(source);
    }
    if (particles.size() == oldSize)
        return;
    invalidateSPHAdvection();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (mode() != executionMode::CPU)
    {
        SPHParticles_.appendHostToDevice(stream);
        const Real support = 2.0 * SPHProperties_.smoothingLength_;
        setSPHSpatialGrid(SPHSpatialGrid_, SPHParticles_, support, minimumBoundary(), maximumBoundary(), SPHPeriodicDomain_, stream);
    }
#else
    (void)stream;
#endif
}

void SPHDEM::appendSPHJetLayer(SPHJetState::constraint& source)
{
    const Vec3 center = source.value_.outletCenter() + (source.value_.speed() * (time() - source.startTime_) - (source.layerCount() + 0.5) * SPHProperties_.spacing_) * source.value_.direction();
    for (const Vec3& offset : source.offsets_)
    {
        SPHParticle particle{center + offset, source.value_.velocity()};
        particle.setMass(SPHProperties_.particleMass());
        particle.setDensity(SPHProperties_.referenceDensity_);
        particle.setConstrained(true);
        source.particleIndices_.push_back(static_cast<int>(SPHParticles_.hostSize()));
        SPHParticles_.host().push_back(std::move(particle));
    }
}

int SPHDEM::SPHJetInterval(int plannedSteps) const noexcept
{
    if (SPHJets_.completed_)
        return plannedSteps;
    for (const auto& source : SPHJets_.constraints_)
    {
        Real next = source.endTime();
        if (source.layerCount() < source.axialCount_)
            next = std::min(next, source.nextEmissionTime(SPHProperties_.spacing_));
        const Real ratio = (next - time()) / timeStep();
        if (ratio > math::defaultTolerance && ratio < plannedSteps)
            plannedSteps = std::max(1, static_cast<int>(std::ceil(ratio - math::defaultTolerance * std::max(Real{1.0}, ratio))));
    }
    return plannedSteps;
}

bool SPHDEM::updateSPHJetOpenings(cpu::mode)
{
    if (SPHJetOpenings_.hostSize() == 0)
        return false;
    const Real nominalVolume = SPHProperties_.spacing_ * SPHProperties_.spacing_ * SPHProperties_.spacing_;
    const Real clearance = 0.5 * SPHProperties_.spacing_;
    auto& boundaries = virtualParticles_.host();
    const int count = static_cast<int>(boundaries.size());
    int changed = 0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (count >= cpu::parallelParticleThreshold) reduction(| : changed)
#endif
    for (int index = 0; index < count; ++index)
    {
        auto& wall = boundaries[index];
        const bool cut = std::any_of(SPHJetOpenings_.host().begin(),
                                     SPHJetOpenings_.host().end(),
                                     [&](const auto& opening) { return execution::isInsideSPHJetOpening(wall.position(), opening.value_, clearance); });
        const Real volume = cut ? 0.0 : nominalVolume;
        if (volume != wall.volume())
        {
            wall.*virtualParticle::volumeField::member = volume;
            wall.setForce(Vec3::zero());
            wall.setPriorForce(Vec3::zero());
            wall.setCouplingTorque(Vec3::zero());
            changed = 1;
        }
    }
    return changed != 0;
}

void SPHDEM::beginOutputSnapshot()
{
    // solve() may already retain the synchronized host view for continuation.
    if (outputSnapshot_.active_)
        return;
    outputObservationSnapshot_ = true;
    captureSPHState(gpuStream());
}

void SPHDEM::endOutputSnapshot()
{
    if (outputObservationSnapshot_)
    {
        outputObservationSnapshot_ = false;
        restoreSPHState();
    }
}

void SPHDEM::beginSolveCompletionSnapshot()
{
    if (!outputSnapshot_.active_)
        captureSPHState(gpuStream());
}

void SPHDEM::captureSPHState(cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (mode() != executionMode::CPU)
    {
        LSDEM::copySystemDeviceToHost(stream);
        SPHParticles_.copyDeviceToHostAsync(stream);
        virtualParticles_.copyDeviceToHostAsync(stream);
        if (SPHStepState_.couplingActive())
            SPHIntervalStart_.copyDeviceToHostAsync(stream);
        host_device_detail::synchronize(stream);
    }
#else
    (void)stream;
#endif
    outputSnapshot_.SPHParticles_ = SPHParticles_.host();
    outputSnapshot_.virtualParticles_ = virtualParticles_.host();
    outputSnapshot_.active_ = true;

    if (SPHStepState_.couplingActive())
    {
        const Real duration = SPHStepState_.couplingDEMSteps_ * timeStep();
        const Real fraction = Real(SPHStepState_.completedDEMSteps_) / SPHStepState_.couplingDEMSteps_;
        const auto& start = SPHIntervalStart_.host();
        const auto& end = outputSnapshot_.SPHParticles_;
        auto& visible = SPHParticles_.host();
        // DKD continuous output: linear endpoint blend minus the constant-
        // acceleration curvature. No fluid force or integration is repeated.
        const int particleCount = static_cast<int>(visible.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= cpu::parallelParticleThreshold)
#endif
        for (int index = 0; index < particleCount; ++index)
        {
            SPHParticle& value = visible[index];
            const SPHParticle& a = start[index];
            const SPHParticle& b = end[index];
            value.setPosition((1.0 - fraction) * a.position() + fraction * b.position() - 0.5 * duration * fraction * (1.0 - fraction) * (b.velocity() - a.velocity()));
            value.setVelocity((1.0 - fraction) * a.velocity() + fraction * b.velocity());
            value.setDensity((1.0 - fraction) * a.density() + fraction * b.density());
            value.setPressure(execution::computePressureFromDensity(value.density(), SPHProperties_.referenceDensity_, SPHProperties_.soundSpeed_));
            value.setDensityRate((b.density() - a.density()) / duration);
        }
    }
    // Inlet flags and virtual geometry belong to visible time(), not the
    // leading endpoint. They never alter the cached force or device state.
    for (const SPHJetState::constraint& source : SPHJets_.constraints_)
    {
        const SPHJet& jet = source.value_;
        for (int layer = source.firstInletLayer_; layer < source.layerCount(); ++layer)
            for (int local = 0; local < source.layerSize(); ++local)
            {
                SPHParticle& value = SPHParticles_.host()[source.layerBegin(layer) + local];
                const bool constrained = time() < source.endTime() && time() < source.releaseTime(layer, SPHProperties_.spacing_) &&
                                         execution::isInsideSPHJetPipe(value.position(), jet.outletCenter(), jet.direction(), jet.radius(), source.bufferLayerCount_ * SPHProperties_.spacing_);
                value.setConstrained(constrained);
                if (constrained)
                {
                    value.setVelocity(jet.velocity());
                    value.setDensity(SPHProperties_.referenceDensity_);
                    value.setPressure(0.0);
                    value.setDensityRate(0.0);
                }
            }
    }
    if (!SPHBoundaryStatic_ && mode() != executionMode::GPU)
        virtualParticleCoupling_.prepareKinematics(virtualParticles_, LSParticles(), gravity(), SPHStepState_.completedDEMSteps_);
    for (virtualParticle& value : virtualParticles_.host())
    {
        const LSParticle& owner = LSParticles().host()[value.ownerLSParticleIndex()];
        Vec3 position;
        Vec3 normal;
        execution::calculateVirtualParticlePositionAndNormal(position, normal, value.localPosition(), value.localNormal(), owner.position(), owner.orientation());
        value.setPosition(position);
        value.setNormal(normal);
        value.setVelocity(owner.velocity() + math::cross(owner.angularVelocity(), position - owner.position()));
        if (!SPHBoundaryStatic_ && mode() == executionMode::GPU)
        {
            if (SPHStepState_.completedDEMSteps_ > 0)
                value.setAcceleration(value.acceleration() / Real(SPHStepState_.completedDEMSteps_));
            else
            {
                Vec3 velocity = Vec3::zero();
                Vec3 acceleration = Vec3::zero();
                execution::accumulateVirtualParticleVelocityAndAcceleration(velocity,
                                                                            acceleration,
                                                                            value.localPosition(),
                                                                            owner.orientation(),
                                                                            owner.velocity(),
                                                                            owner.angularVelocity(),
                                                                            owner.force(),
                                                                            owner.torque(),
                                                                            owner.inverseMass(),
                                                                            owner.inertiaTensor(),
                                                                            owner.inverseInertiaTensor(),
                                                                            gravity());
                value.setAcceleration(acceleration);
            }
        }
    }
    // Only this reversible Host view follows the visible wall pose. Never
    // invalidate the real search phase or alter device loads from observation.
    updateSPHJetOpenings(cpu::mode{});
}

void SPHDEM::restoreSPHState() noexcept
{
    if (!outputSnapshot_.active_)
    {
        return;
    }

    SPHParticles_.host().swap(outputSnapshot_.SPHParticles_);
    virtualParticles_.host().swap(outputSnapshot_.virtualParticles_);
    outputSnapshot_.active_ = false;
}

void SPHDEM::writeSystemVTU(int frameIndex)
{
    LSDEM::writeSystemVTU(frameIndex);
    writeSPHVTU(frameIndex);
}

void SPHDEM::writeSPHVTU(int frameIndex)
{
    if (SPHParticles_.hostSize() > 0)
    {
        makeSPHParticleVTUWriter(makeVTUFileName("particle", "SPHParticle", frameIndex), SPHParticles_, gravity(), SPHParticleVTUFields_).write(vtuOutputFormat());
    }
}

double SPHDEM::SPHDeviceMemoryGB() const noexcept
{
    return SPHParticles_.deviceMemoryGB() + SPHIntervalStart_.deviceMemoryGB() + virtualParticles_.deviceMemoryGB() + SPHSpatialGrid_.deviceMemoryGB() + virtualParticleSpatialGrid_.deviceMemoryGB() +
           SPHNeighborhood_.deviceMemoryGB() + SPHJetOpenings_.deviceMemoryGB();
}

double SPHDEM::systemDeviceMemoryGB() const noexcept { return LSDEM::systemDeviceMemoryGB() + SPHDeviceMemoryGB(); }

void SPHDEM::prepareSPHAdvectionStep(cpu::mode)
{
    updateSPHNeighborSearchRadius(SPHInteractions_->kinematicsStatistics(virtualParticles_));
    const auto parameters = SPHProperties_.interactionParameters(gravity(), SPHPeriodicDomain_);
    invalidateSPHNeighborhood();
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updateFreeSurface(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->reinitializeDensity(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->updatePriorForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
}

void SPHDEM::updateSPHNeighborSearchRadius(const SPHKinematicsStatistics& boundaryStatistics)
{
    if (boundaryStatistics.invalidValueCount_ > 0)
    {
        throw std::runtime_error("The SPH boundary contains a non-finite velocity or acceleration.");
    }

    const Real fluidVelocityScale = std::max(SPHStepState_.statistics_.maximumVelocity_, SPHProperties_.maximumVelocity_);
    const Real timeHorizon = SPHStepState_.advection_.timeStep_ > 0.0 ? SPHStepState_.advection_.timeStep_ : SPHStepState_.acoustic_.timeStep_;
    const Real fluidSearchBuffer = 2.0 * fluidVelocityScale * timeHorizon + SPHStepState_.statistics_.maximumAcceleration_ * timeHorizon * timeHorizon;
    const Real boundarySearchBuffer = (fluidVelocityScale + boundaryStatistics.maximumVelocity_) * timeHorizon +
                                      0.5 * (SPHStepState_.statistics_.maximumAcceleration_ + boundaryStatistics.maximumAcceleration_) * timeHorizon * timeHorizon;
    const Real searchBuffer = std::max(fluidSearchBuffer, boundarySearchBuffer);
    SPHStepState_.neighborSearchRadius_ = 2.0 * SPHProperties_.smoothingLength_ + searchBuffer;
}

void SPHDEM::applySPHJetVelocity(Real representedTime, cpu::mode)
{
    if (SPHJets_.completed_)
    {
        return;
    }

    auto& particleHost = SPHParticles_.host();
    for (const SPHJetState::constraint& source : SPHJets_.constraints_)
    {
        const SPHJet& jet = source.value_;
        const Vec3 velocity = jet.velocity();
        const int layerCount = source.layerSize();
        for (int layer = source.firstInletLayer_; layer < source.layerCount(); ++layer)
        {
            const bool active = representedTime < source.endTime() && representedTime < source.releaseTime(layer, SPHProperties_.spacing_);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (layerCount >= cpu::parallelParticleThreshold)
#endif
            for (int localParticleIndex = 0; localParticleIndex < layerCount; ++localParticleIndex)
            {
                SPHParticle& particle = particleHost[source.layerBegin(layer) + localParticleIndex];
                const bool constrained =
                    active && execution::isInsideSPHJetPipe(particle.position(), jet.outletCenter(), jet.direction(), jet.radius(), source.bufferLayerCount_ * SPHProperties_.spacing_);
                particle.setConstrained(constrained);
                if (constrained)
                {
                    particle.setVelocity(velocity);
                    particle.setDensity(SPHProperties_.referenceDensity_);
                    particle.setPressure(0.0);
                    particle.setDensityRate(0.0);
                }
            }
        }
    }
}

void SPHDEM::ensureSPHNeighborhood(cpu::mode)
{
    const Real supportRadius = 2.0 * SPHProperties_.smoothingLength_;
    const Real skin = std::max(Real{0.0}, SPHStepState_.neighborSearchRadius_ - supportRadius);
    if (!SPHNeighborhood_.needsRebuild(SPHParticles_, virtualParticles_, skin))
        return;
    const Vec3 padding{supportRadius, supportRadius, supportRadius};
    SPHInteractions_
        ->buildNeighborhood(SPHParticles_, virtualParticles_, minimumBoundary() - padding, maximumBoundary() + padding, supportRadius, SPHStepState_.neighborSearchRadius_, SPHPeriodicDomain_);
    SPHNeighborhood_.captureHost(SPHParticles_, virtualParticles_);
}

void SPHDEM::advanceSPH(Real timeStep, cpu::mode)
{
    // The preceding interval's hydraulic wrench has already been removed.
    if (!SPHBoundaryStatic_ && updateSPHJetOpenings(cpu::mode{}))
        invalidateSPHAdvection();
    const Real jetConstraintTime = time();
    applySPHJetVelocity(jetConstraintTime, cpu::mode{});
    if (SPHStepState_.timeInAdvectionStep_ <= math::defaultTolerance)
    {
        prepareSPHAdvectionStep(cpu::mode{});
    }

    SPHIntervalStart_.host() = SPHParticles_.host();
    const Real halfTimeStep = 0.5 * timeStep;
    const auto parameters = SPHProperties_.interactionParameters(gravity(), SPHPeriodicDomain_);
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updateDensity(SPHParticles_, virtualParticles_, parameters, halfTimeStep);
    SPHInteractions_->integratePosition(SPHParticles_, halfTimeStep);
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updatePressureForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->integrateVelocity(SPHParticles_, gravity(), timeStep);
    applySPHJetVelocity(jetConstraintTime + halfTimeStep, cpu::mode{});
    SPHInteractions_->integratePosition(SPHParticles_, halfTimeStep);
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updateDensity(SPHParticles_, virtualParticles_, parameters, halfTimeStep);

    finishSPHInterval(timeStep, nullptr);
}

void SPHDEM::advanceCouplingInterval(Real interval, cpu::mode)
{
    if (!SPHBoundaryStatic_)
    {
        virtualParticleCoupling_.prepareKinematics(virtualParticles_, LSParticles(), gravity(), SPHStepState_.completedDEMSteps_);
    }
    virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles(), -1.0);
    advanceSPH(interval, cpu::mode{});
    virtualParticleCoupling_.collectForceAndTorque(virtualParticles_, LSParticles());
    virtualParticleCoupling_.clearKinematics();
}

void SPHDEM::finishSPHInterval(Real interval, cudaStream_t stream)
{
    // Called only from the common acoustic START, before the DEM clock moves.
    const Real endpoint = time() + interval;
    const bool hostFluid = mode() == executionMode::CPU;
    if (!SPHJets_.completed_)
    {
        if (hostFluid)
            applySPHJetVelocity(endpoint, cpu::mode{});
        else
            applySPHJetVelocity(endpoint, gpu::mode{}, stream);
        SPHJets_.completed_ = SPHJets_.completedAt(endpoint);
    }
    const bool advectionCompleted = SPHStepState_.advanceAdvectionPhase(interval);
    if (hostFluid)
        updateSPHAcousticTimeStep(cpu::mode{});
    else
        updateSPHAcousticTimeStep(gpu::mode{}, stream);
    if (advectionCompleted)
        updateSPHAdvectionTimeStep();
}

void SPHDEM::updateSPHAcousticTimeStep(cpu::mode)
{
    const auto statistics =
        SPHParticles_.hostSize() > 0 && SPHProperties_.smoothingLength_ > math::defaultTolerance ? SPHInteractions_->stateStatistics(SPHParticles_, gravity()) : SPHStateStatistics{};
    updateSPHAcousticTimeStep(statistics);
}

void SPHDEM::updateSPHAcousticTimeStep(const SPHStateStatistics& statistics)
{
    if (SPHParticles_.hostSize() == 0 || SPHProperties_.smoothingLength_ <= math::defaultTolerance)
    {
        SPHStepState_.acoustic_ = {timeStep(), timeStep(), 1};
        SPHStepState_.statistics_ = {};
        return;
    }

    if (statistics.invalidValueCount_ > 0)
    {
        throw std::runtime_error("The SPH state contains a non-finite velocity, acceleration or density.");
    }
    const Real minimumDensity = execution::minimumSPHDensityRatio * SPHProperties_.referenceDensity_;
    const Real maximumDensity = execution::maximumSPHDensityRatio * SPHProperties_.referenceDensity_;
    if (statistics.minimumDensity_ < minimumDensity || statistics.maximumDensity_ > maximumDensity)
    {
        throw std::runtime_error("The SPH density left the admissible weakly-compressible range. Reduce the time step or increase the sound speed.");
    }

    SPHStepState_.statistics_ = statistics;
    SPHStepState_.acoustic_.limit_ = execution::calculateSPHAcousticTimeStep(SPHProperties_.smoothingLength_, SPHProperties_.soundSpeed_, SPHStepState_.statistics_.maximumVelocity_);
    if (!math::isFinite(SPHStepState_.acoustic_.limit_) || SPHStepState_.acoustic_.limit_ <= 0.0)
    {
        throw std::runtime_error("Failed to calculate a finite positive SPH acoustic time-step limit.");
    }

    Real maximumSPHTimeStep = SPHStepState_.acoustic_.limit_;
    if (SPHStepState_.timeInAdvectionStep_ > math::defaultTolerance)
    {
        const Real remainingAdvectionTime = SPHStepState_.advection_.timeStep_ - SPHStepState_.timeInAdvectionStep_;
        maximumSPHTimeStep = remainingAdvectionTime < maximumSPHTimeStep ? remainingAdvectionTime : maximumSPHTimeStep;
    }
    setSPHAcousticTimeStepFromLimit(maximumSPHTimeStep);
}

void SPHDEM::updateSPHAdvectionTimeStep()
{
    if (SPHParticles_.hostSize() == 0 || SPHProperties_.smoothingLength_ <= math::defaultTolerance)
    {
        SPHStepState_.advection_.limit_ = timeStep();
        SPHStepState_.advection_.timeStep_ = timeStep();
        SPHStepState_.timeInAdvectionStep_ = 0.0;
        SPHStepState_.advection_.interval_ = 1;
        return;
    }

    SPHStepState_.advection_.limit_ = execution::calculateSPHAdvectionTimeStepFromLimits(SPHProperties_.smoothingLength_,
                                                                                         SPHProperties_.viscousTimeScale_,
                                                                                         SPHStepState_.statistics_.maximumVelocity_,
                                                                                         SPHStepState_.statistics_.maximumAcceleration_,
                                                                                         SPHProperties_.maximumVelocity_);
    if (!math::isFinite(SPHStepState_.advection_.limit_) || SPHStepState_.advection_.limit_ <= 0.0)
    {
        throw std::runtime_error("Failed to calculate a finite positive SPH advection time-step limit.");
    }

    const Real maximumSPHTimeStep = SPHStepState_.advection_.limit_ < SPHStepState_.acoustic_.limit_ ? SPHStepState_.advection_.limit_ : SPHStepState_.acoustic_.limit_;
    setSPHAcousticTimeStepFromLimit(maximumSPHTimeStep);
    const Real maximumAdvectionStepInterval = std::floor(SPHStepState_.advection_.limit_ / SPHStepState_.acoustic_.timeStep_ + math::defaultTolerance);
    SPHStepState_.advection_.interval_ = static_cast<int>(std::clamp(maximumAdvectionStepInterval, Real{1.0}, Real{std::numeric_limits<int>::max()}));
    SPHStepState_.advection_.timeStep_ = SPHStepState_.advection_.interval_ * SPHStepState_.acoustic_.timeStep_;
}

void SPHDEM::setSPHAcousticTimeStepFromLimit(Real maximumTimeStep)
{
    const Real DEMTimeStep = timeStep();
    // Use relative slack: an absolute tolerance in seconds can exceed a small stability limit.
    if (DEMTimeStep > maximumTimeStep && DEMTimeStep - maximumTimeStep > math::defaultTolerance * maximumTimeStep)
    {
        throw std::runtime_error("The DEM time step exceeds the SPH acoustic or advection stability limit. Reduce the DEM time step.");
    }
    const Real maximumStepInterval = std::floor(maximumTimeStep / DEMTimeStep + math::defaultTolerance);
    SPHStepState_.acoustic_.interval_ = static_cast<int>(std::clamp(maximumStepInterval, Real{1.0}, Real{std::numeric_limits<int>::max()}));
    SPHStepState_.acoustic_.timeStep_ = SPHStepState_.acoustic_.interval_ * DEMTimeStep;
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

bool SPHDEM::updateSPHJetOpenings(gpu::mode, cudaStream_t stream)
{
    if (SPHJetOpenings_.deviceSize() == 0)
        return false;
    const Real spacing = SPHProperties_.spacing_;
    return cuda::launchUpdateSPHJetOpenings(virtualParticles_, SPHJetOpenings_, spacing * spacing * spacing, 0.5 * spacing, stream);
}

void SPHDEM::applySPHJetVelocity(Real representedTime, gpu::mode, cudaStream_t stream)
{
    if (SPHJets_.completed_)
    {
        return;
    }

    for (const SPHJetState::constraint& source : SPHJets_.constraints_)
    {
        const SPHJet& jet = source.value_;
        for (int layer = source.firstInletLayer_; layer < source.layerCount(); ++layer)
        {
            const bool active = representedTime < source.endTime() && representedTime < source.releaseTime(layer, SPHProperties_.spacing_);
            cuda::launchApplySPHJetVelocity(SPHParticles_,
                                            source.layerBegin(layer),
                                            source.layerSize(),
                                            active,
                                            jet.outletCenter(),
                                            jet.direction(),
                                            jet.radius(),
                                            source.bufferLayerCount_ * SPHProperties_.spacing_,
                                            jet.velocity(),
                                            SPHProperties_.referenceDensity_,
                                            stream);
        }
    }
}

void SPHDEM::initializeSPHDevice(cudaStream_t stream)
{
    initializeSPHState();
    SPHParticles_.copyHostToDeviceAsync(stream);
    virtualParticles_.copyHostToDeviceAsync(stream);
    SPHJetOpenings_.copyHostToDeviceAsync(stream);

    const Real supportRadius = SPHProperties_.smoothingLength_ > math::defaultTolerance ? 2.0 * SPHProperties_.smoothingLength_ : math::norm(maximumBoundary() - minimumBoundary());
    const auto parameters = SPHProperties_.interactionParameters(gravity(), SPHPeriodicDomain_);
    setSPHSpatialGrid(SPHSpatialGrid_, SPHParticles_, supportRadius, minimumBoundary(), maximumBoundary(), SPHPeriodicDomain_, stream);
    setSPHSpatialGrid(virtualParticleSpatialGrid_, virtualParticles_, supportRadius, minimumBoundary(), maximumBoundary(), SPHPeriodicDomain_, stream);

    SPHStepState_.neighborSearchRadius_ = supportRadius;
    cuda::launchInitializeSPHParticles(SPHParticles_, SPHProperties_.referenceDensity_, SPHProperties_.soundSpeed_, stream);
    if (LSParticles().usesDevice())
    {
        cuda::launchUpdateVirtualParticlePositionAndNormal(virtualParticles_, LSParticles(), stream);
        cuda::launchClearVirtualParticleKinematics(virtualParticles_, stream);
        cuda::launchAccumulateVirtualParticleKinematics(virtualParticles_, LSParticles(), gravity(), stream);
    }
    cuda::launchBuildSPHSpatialGrid(SPHSpatialGrid_, SPHParticles_, SPHPeriodicDomain_, stream);
    cuda::launchBuildVirtualParticleSpatialGrid(virtualParticleSpatialGrid_, virtualParticles_, SPHPeriodicDomain_, stream);
    cuda::launchUpdateSPHPriorForceAndBoundaryForce(SPHParticles_,
                                                    virtualParticles_,
                                                    SPHSpatialGrid_.device(),
                                                    virtualParticleSpatialGrid_.device(),
                                                    parameters,
                                                    SPHStepState_.neighborSearchRadius_,
                                                    stream);
    cuda::launchUpdateSPHPressureForceAndBoundaryForce(SPHParticles_,
                                                       virtualParticles_,
                                                       SPHSpatialGrid_.device(),
                                                       virtualParticleSpatialGrid_.device(),
                                                       parameters,
                                                       SPHStepState_.neighborSearchRadius_,
                                                       stream);
    cuda::launchClearVirtualParticleKinematics(virtualParticles_, stream);
    updateSPHAcousticTimeStep(gpu::mode{}, stream);
    updateSPHAdvectionTimeStep();
    if (LSParticles().usesDevice())
    {
        cuda::launchRefreshVirtualParticleCouplingLoads(virtualParticles_, LSParticles(), stream);
    }
    else
    {
        virtualParticleCoupling_.copyForceAndTorqueFromDevice(virtualParticles_, LSParticles(), stream);
    }
}

void SPHDEM::prepareSPHAdvectionStep(gpu::mode, cudaStream_t stream)
{
    const auto parameters = SPHProperties_.interactionParameters(gravity(), SPHPeriodicDomain_);
    updateSPHNeighborSearchRadius(cuda::calculateSPHKinematicsStatistics(virtualParticles_, stream));

    invalidateSPHNeighborhood();
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHFreeSurface(SPHParticles_, virtualParticles_, SPHSpatialGrid_.device(), virtualParticleSpatialGrid_.device(), parameters, SPHStepState_.neighborSearchRadius_, stream);
    cuda::launchReinitializeSPHDensity(SPHParticles_, virtualParticles_, SPHSpatialGrid_.device(), virtualParticleSpatialGrid_.device(), parameters, SPHStepState_.neighborSearchRadius_, stream);
    cuda::launchUpdateSPHPriorForceAndBoundaryForce(SPHParticles_,
                                                    virtualParticles_,
                                                    SPHSpatialGrid_.device(),
                                                    virtualParticleSpatialGrid_.device(),
                                                    parameters,
                                                    SPHStepState_.neighborSearchRadius_,
                                                    stream);
}

void SPHDEM::ensureSPHNeighborhood(gpu::mode, cudaStream_t stream)
{
    const Real skin = std::max(Real{0.0}, SPHStepState_.neighborSearchRadius_ - 2.0 * SPHProperties_.smoothingLength_);
    if (SPHNeighborhood_.valid())
    {
        const Real displacement = cuda::maximumSPHNeighborDisplacementSquared(SPHParticles_, virtualParticles_, SPHNeighborhood_.positions(), stream);
        if (!math::isFinite(displacement))
            throw std::runtime_error("Non-finite SPH position while checking the neighborhood displacement.");
        if (!SPHNeighborhood_.exceedsSkin(displacement, skin))
            return;
    }
    cuda::launchBuildSPHSpatialGrid(SPHSpatialGrid_, SPHParticles_, SPHPeriodicDomain_, stream);
    cuda::launchBuildVirtualParticleSpatialGrid(virtualParticleSpatialGrid_, virtualParticles_, SPHPeriodicDomain_, stream);
    SPHNeighborhood_.captureDevice(SPHParticles_, virtualParticles_, stream);
}

void SPHDEM::advanceSPH(Real timeStep, gpu::mode, cudaStream_t stream)
{
    if (!SPHBoundaryStatic_ && updateSPHJetOpenings(gpu::mode{}, stream))
        invalidateSPHAdvection();
    const Real jetConstraintTime = time();
    applySPHJetVelocity(jetConstraintTime, gpu::mode{}, stream);
    if (SPHStepState_.timeInAdvectionStep_ <= math::defaultTolerance)
    {
        prepareSPHAdvectionStep(gpu::mode{}, stream);
    }

    copySPHDeviceState(SPHIntervalStart_, SPHParticles_, SPHIntervalStorage::layout_type{}, stream);
    const Real halfTimeStep = 0.5 * timeStep;
    const auto parameters = SPHProperties_.interactionParameters(gravity(), SPHPeriodicDomain_);
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHDensity(SPHParticles_,
                                 virtualParticles_,
                                 SPHSpatialGrid_.device(),
                                 virtualParticleSpatialGrid_.device(),
                                 parameters,
                                 SPHStepState_.neighborSearchRadius_,
                                 halfTimeStep,
                                 stream);
    cuda::launchIntegrateSPHPosition(SPHParticles_, halfTimeStep, stream);
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHPressureForceAndBoundaryForce(SPHParticles_,
                                                       virtualParticles_,
                                                       SPHSpatialGrid_.device(),
                                                       virtualParticleSpatialGrid_.device(),
                                                       parameters,
                                                       SPHStepState_.neighborSearchRadius_,
                                                       stream);
    cuda::launchIntegrateSPHVelocity(SPHParticles_, gravity(), timeStep, stream);
    applySPHJetVelocity(jetConstraintTime + halfTimeStep, gpu::mode{}, stream);
    cuda::launchIntegrateSPHPosition(SPHParticles_, halfTimeStep, stream);
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHDensity(SPHParticles_,
                                 virtualParticles_,
                                 SPHSpatialGrid_.device(),
                                 virtualParticleSpatialGrid_.device(),
                                 parameters,
                                 SPHStepState_.neighborSearchRadius_,
                                 halfTimeStep,
                                 stream);

    finishSPHInterval(timeStep, stream);
}

void SPHDEM::advanceCouplingInterval(Real interval, gpu::mode, cudaStream_t stream)
{
    if (!SPHBoundaryStatic_)
    {
        cuda::launchPrepareVirtualParticleKinematics(virtualParticles_, LSParticles(), gravity(), SPHStepState_.completedDEMSteps_, stream);
    }
    cuda::launchAddVirtualParticleForceAndTorque(mutableLSParticles(), virtualParticles_, -1.0, stream);
    advanceSPH(interval, gpu::mode{}, stream);
    cuda::launchRefreshVirtualParticleCouplingLoads(virtualParticles_, LSParticles(), stream);
    if (!SPHBoundaryStatic_)
    {
        cuda::launchClearVirtualParticleKinematics(virtualParticles_, stream);
    }
}

void SPHDEM::advanceCouplingInterval(Real interval, hybrid::mode, cudaStream_t stream)
{
    if (!SPHBoundaryStatic_)
    {
        virtualParticleCoupling_.prepareDeviceKinematics(virtualParticles_, LSParticles(), gravity(), SPHStepState_.completedDEMSteps_, stream);
    }
    virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles(), -1.0);
    advanceSPH(interval, gpu::mode{}, stream);
    virtualParticleCoupling_.copyForceAndTorqueFromDevice(virtualParticles_, LSParticles(), stream);
    virtualParticleCoupling_.clearKinematics();
}

void SPHDEM::updateSPHAcousticTimeStep(gpu::mode, cudaStream_t stream)
{
    const auto statistics =
        SPHParticles_.hostSize() > 0 && SPHProperties_.smoothingLength_ > math::defaultTolerance ? cuda::calculateSPHStateStatistics(SPHParticles_, gravity(), stream) : SPHStateStatistics{};
    updateSPHAcousticTimeStep(statistics);
}

void SPHDEM::copySystemDeviceToHost(cudaStream_t stream)
{
    if (!outputSnapshot_.active_)
        captureSPHState(stream);
    else
        LSDEM::copySystemDeviceToHost(stream);
}

#else

void SPHDEM::initializeSPHDevice(cudaStream_t) { throw std::runtime_error("The SPH GPU backend is unavailable in this FunDEM build."); }
void SPHDEM::prepareSPHAdvectionStep(gpu::mode, cudaStream_t) {}
void SPHDEM::applySPHJetVelocity(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::advanceSPH(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::advanceCouplingInterval(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::advanceCouplingInterval(Real, hybrid::mode, cudaStream_t) {}
void SPHDEM::updateSPHAcousticTimeStep(gpu::mode, cudaStream_t) {}
void SPHDEM::copySystemDeviceToHost(cudaStream_t) {}

#endif

} // namespace fundem
