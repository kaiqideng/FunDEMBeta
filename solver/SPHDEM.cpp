#include "SPHDEM.h"

#include "data/energyOutput.h"
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

void orthogonalBasis(const math::Vec3& unitAxis, math::Vec3& tangent1, math::Vec3& tangent2) noexcept
{
    tangent1 = std::abs(unitAxis.x) > std::abs(unitAxis.z) ? math::Vec3{-unitAxis.y, unitAxis.x, 0.0} : math::Vec3{0.0, -unitAxis.z, unitAxis.y};
    tangent1 = math::normalizedOrZero(tangent1);
    tangent2 = math::cross(unitAxis, tangent1);
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

int SPHDEM::addSPHParticle(const SPHParticle& value)
{
    if (initialized())
    {
        throw std::logic_error("SPH particles cannot be added after solver initialization.");
    }
    if (!value.isValid())
    {
        throw std::invalid_argument("Cannot add an invalid SPH particle.");
    }
    if (!SPHProperties_.set_)
    {
        throw std::logic_error("Set the SPH properties before adding SPH particles.");
    }

    const int particleIndex = static_cast<int>(SPHParticles_.hostSize());
    SPHParticle newParticle = value;
    newParticle.setMass(SPHProperties_.particleMass());
    newParticle.setDensity(SPHProperties_.referenceDensity_);
    SPHParticles_.host().push_back(std::move(newParticle));
    invalidateContinuation();
    return particleIndex;
}

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

    Vec3 tangent1;
    Vec3 tangent2;
    orthogonalBasis(axis, tangent1, tangent2);
    const int radialIndex = discretization.radialIndex_;
    const int axialCount = discretization.axialCount_;

    std::vector<SPHParticle> particles;
    particles.reserve(static_cast<std::size_t>(discretization.particleCount_));
    for (int axialIndex = 0; axialIndex < axialCount; ++axialIndex)
    {
        const Vec3 axialPosition = inletCenter + (Real(axialIndex) + 0.5) * SPHProperties_.spacing_ * axis;
        for (int y = -radialIndex; y <= radialIndex; ++y)
        {
            for (int x = -radialIndex; x <= radialIndex; ++x)
            {
                const Real offset1 = x * SPHProperties_.spacing_;
                const Real offset2 = y * SPHProperties_.spacing_;
                if (discretization.contains(x, y))
                {
                    particles.emplace_back(axialPosition + offset1 * tangent1 + offset2 * tangent2, velocity);
                }
            }
        }
    }
    return addSPHParticles(std::move(particles));
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
    const Real endTime = time() + value.duration();
    execution::SPHJetDiscretization discretization;
    if (!execution::trySPHJetDiscretization(SPHProperties_.spacing_, value.radius(), pipeLength, discretization))
    {
        throw std::invalid_argument("The SPH jet must contain at least one particle layer.");
    }

    Vec3 tangent1;
    Vec3 tangent2;
    orthogonalBasis(value.direction(), tangent1, tangent2);
    const Vec3 velocity = value.velocity();
    std::vector<SPHParticle> particles;
    particles.reserve(static_cast<std::size_t>(discretization.particleCount_));
    for (int axialIndex = 0; axialIndex < discretization.axialCount_; ++axialIndex)
    {
        const Vec3 axialPosition = value.outletCenter() - (Real(axialIndex) + 0.5) * SPHProperties_.spacing_ * value.direction();
        for (int y = -discretization.radialIndex_; y <= discretization.radialIndex_; ++y)
        {
            for (int x = -discretization.radialIndex_; x <= discretization.radialIndex_; ++x)
            {
                if (discretization.contains(x, y))
                {
                    particles.emplace_back(axialPosition + Real(x) * SPHProperties_.spacing_ * tangent1 + Real(y) * SPHProperties_.spacing_ * tangent2, velocity);
                }
            }
        }
    }

    SPHJets_.reserve(SPHJets_.size() + 1);
    const int particleCount = static_cast<int>(particles.size());
    const int firstParticleIndex = addSPHParticles(std::move(particles));
    SPHJets_.push_back({value, pipeLength, endTime, firstParticleIndex, particleCount});
    SPHStepState_.latestSPHJetEndTime_ = endTime > SPHStepState_.latestSPHJetEndTime_ ? endTime : SPHStepState_.latestSPHJetEndTime_;
    SPHStepState_.jetsCompleted_ = false;
    return firstParticleIndex;
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
    SPHStepState_.resume_ = false;
    outputSnapshot_ = {};
    SPHInitializationState_.fluidParticleCount_ = -1;
    SPHInitializationState_.boundaryParticleCount_ = -1;
}

void SPHDEM::resumeContinuation() noexcept
{
    if (mode() == executionMode::CPU && outputSnapshot_.active_)
    {
        outputSnapshot_.preserveCurrentHost_ = false;
        restoreSPHState(nullptr);
    }
    SPHStepState_.resume_ = false;
}

void SPHDEM::setSPHParticleVTUFields(std::vector<SPHParticleVTUField> fields) { SPHParticleVTUFields_ = std::move(fields); }

energyRecord SPHDEM::systemEnergy() const { return LSDEM::systemEnergy(); }

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
    initializeSPH(cpu::mode{});
    assembleForceAndTorque(0.0, cpu::mode{});
    virtualParticleCoupling_.resetKinematicsReference(virtualParticles_, LSParticles());
}

void SPHDEM::initialize(gpu::mode, cudaStream_t stream)
{
    initializeLSParticleDevice(stream);
    initializeSPH(gpu::mode{}, stream);
    assembleForceAndTorque(0.0, gpu::mode{}, stream);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    cuda::launchResetVirtualParticleKinematicsReference(virtualParticles_, LSParticles(), stream);
#endif
}

void SPHDEM::initialize(hybrid::mode, cudaStream_t stream)
{
    initializeSPH(hybrid::mode{}, stream);
    assembleForceAndTorque(0.0, hybrid::mode{}, stream);
    virtualParticleCoupling_.resetKinematicsReference(virtualParticles_, LSParticles());
}

void SPHDEM::advanceCPU(Real timeStep)
{
    // Sample the complete previous endpoint loads, before either half kick or
    // force clearing. This includes retained fluid and external-module loads.
    if (!SPHInitializationState_.boundaryStatic_)
        virtualParticleCoupling_.accumulateKinematics(virtualParticles_, LSParticles(), gravity(), timeStep);
    const Real halfTimeStep = 0.5 * timeStep;
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cpu::integratePosition(mutableLSParticles(), timeStep);
    assembleForceAndTorque(timeStep, cpu::mode{});
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    if (SPHStepState_.couplingImpulseTime_ > 0.0)
    {
        virtualParticleCoupling_.applyImpulseCorrection(mutableLSParticles(), SPHStepState_.couplingImpulseTime_);
        SPHStepState_.couplingImpulseTime_ = 0.0;
    }
}

void SPHDEM::advanceGPU(Real timeStep)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const cudaStream_t stream = gpuStream();
    const Real halfTimeStep = 0.5 * timeStep;
    if (!SPHInitializationState_.boundaryStatic_)
        cuda::launchAccumulateVirtualParticleKinematics(virtualParticles_, LSParticles(), gravity(), timeStep, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(mutableLSParticles(), gravity(), halfTimeStep, stream);
    cuda::launchPositionAndOrientationIntegration(mutableLSParticles(), timeStep, stream);
    assembleForceAndTorque(timeStep, gpu::mode{}, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(mutableLSParticles(), gravity(), halfTimeStep, stream);
    if (SPHStepState_.couplingImpulseTime_ > 0.0)
    {
        cuda::launchApplyVirtualParticleImpulseCorrection(mutableLSParticles(), virtualParticles_, SPHStepState_.couplingImpulseTime_, stream);
        SPHStepState_.couplingImpulseTime_ = 0.0;
    }
#else
    (void)timeStep;
#endif
}

void SPHDEM::advanceHybrid(Real timeStep)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const cudaStream_t stream = gpuStream();
    const Real halfTimeStep = 0.5 * timeStep;
    if (!SPHInitializationState_.boundaryStatic_)
        virtualParticleCoupling_.accumulateKinematics(virtualParticles_, LSParticles(), gravity(), timeStep);
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cpu::integratePosition(mutableLSParticles(), timeStep);
    assembleForceAndTorque(timeStep, hybrid::mode{}, stream);
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    if (SPHStepState_.couplingImpulseTime_ > 0.0)
    {
        virtualParticleCoupling_.applyImpulseCorrection(mutableLSParticles(), SPHStepState_.couplingImpulseTime_);
        SPHStepState_.couplingImpulseTime_ = 0.0;
    }
#else
    (void)timeStep;
#endif
}

void SPHDEM::assembleForceAndTorque(Real historyTimeStep, cpu::mode)
{
    cpu::clearForceAndTorque(mutableLSParticles());
    calculateLSParticleForceAndTorque(historyTimeStep, cpu::mode{});
    calculateSPHForceAndTorque(historyTimeStep, cpu::mode{});
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
    cuda::launchBuildSpatialGrid(mutableLSParticleSpatialGrid(), LSParticles(), stream);
    calculateLSParticleForceAndTorque(historyTimeStep, gpu::mode{}, stream);
    calculateSPHForceAndTorque(historyTimeStep, gpu::mode{}, stream);
    addLSParticleExternalForceAndTorque(deviceFields(mutableLSParticles()), stream);
#else
    (void)historyTimeStep;
    (void)stream;
#endif
}

void SPHDEM::assembleForceAndTorque(Real historyTimeStep, hybrid::mode, cudaStream_t stream)
{
    cpu::clearForceAndTorque(mutableLSParticles());
    calculateLSParticleForceAndTorque(historyTimeStep, cpu::mode{});
    calculateSPHForceAndTorque(historyTimeStep, hybrid::mode{}, stream);
    addLSParticleExternalForceAndTorque(mutableLSParticles().host());
}

void SPHDEM::initializeSPHState()
{
    const int SPHParticleCount = static_cast<int>(SPHParticles_.hostSize());
    const int LSParticleCount = static_cast<int>(LSParticles().hostSize());
    SPHInitializationState_.boundaryStatic_ =
        !SPHBoundaryMotionExpected_ && std::all_of(LSParticles().host().begin(),
                                                   LSParticles().host().end(),
                                                   [](const LSParticle& value) { return value.inverseMass() == 0.0 && value.velocity() == Vec3::zero() && value.angularVelocity() == Vec3::zero(); });
    const bool particleConfigurationChanged = SPHParticleCount != SPHInitializationState_.fluidParticleCount_ || LSParticleCount != SPHInitializationState_.boundaryParticleCount_;
    SPHStepState_.resume_ = false;
    virtualParticles_.resetDevice();
    virtualParticles_.host().clear();
    SPHStepState_.pendingDEMSteps_ = 0;
    SPHStepState_.couplingImpulseTime_ = 0.0;
    SPHStepState_.acousticTimeStepLimit_ = 0.0;
    SPHStepState_.acousticTimeStep_ = 0.0;
    SPHStepState_.advectionTimeStepLimit_ = 0.0;
    SPHStepState_.advectionTimeStep_ = 0.0;
    SPHStepState_.observedMaximumVelocity_ = 0.0;
    SPHStepState_.observedMaximumAcceleration_ = 0.0;
    SPHStepState_.neighborSearchRadius_ = 2.0 * SPHProperties_.smoothingLength_;
    SPHStepState_.representedTime_ = time();
    SPHStepState_.jetsCompleted_ = SPHJets_.empty() || SPHStepState_.representedTime_ >= SPHStepState_.latestSPHJetEndTime_;
    auto& particleHost = SPHParticles_.host();
    for (SPHParticle& value : particleHost)
    {
        value.setConstrained(false);
    }
    applySPHJetVelocity(SPHStepState_.representedTime_, cpu::mode{});
    if (particleConfigurationChanged)
    {
        SPHStepState_.timeInAdvectionStep_ = 0.0;
    }
    SPHStepState_.acousticStepInterval_ = 1;
    SPHStepState_.advectionStepInterval_ = 1;
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
    if (mode() == executionMode::GPU)
    {
        virtualParticleCoupling_.reset();
    }
    else
    {
        virtualParticleCoupling_.initialize(virtualParticles_, LSParticles(), gravity());
    }
    SPHInitializationState_.fluidParticleCount_ = SPHParticleCount;
    SPHInitializationState_.boundaryParticleCount_ = LSParticleCount;
}

void SPHDEM::initializeSPH(cpu::mode)
{
    initializeSPHState();
    if (!SPHInteractions_)
    {
        SPHInteractions_ = std::make_unique<cpu::SPHInteraction>();
    }
    const Real supportRadius = SPHProperties_.smoothingLength_ > math::defaultTolerance ? 2.0 * SPHProperties_.smoothingLength_ : math::norm(maximumBoundary() - minimumBoundary());
    const Real particleMass = SPHProperties_.particleMass();
    const Vec3 boundaryPadding{supportRadius, supportRadius, supportRadius};
    const cpu::SPHInteractionParameters parameters{SPHProperties_.smoothingLength_,
                                                   particleMass,
                                                   SPHProperties_.referenceDensity_,
                                                   SPHProperties_.latticeKernelSum3D_,
                                                   SPHProperties_.soundSpeed_,
                                                   SPHProperties_.dynamicViscosity_,
                                                   gravity()};

    SPHStepState_.neighborSearchRadius_ = supportRadius;
    SPHStepState_.timeInAdvectionStep_ = 0.0;
    SPHInteractions_->initializeParticles(SPHParticles_, SPHProperties_.referenceDensity_, SPHProperties_.soundSpeed_);
    SPHInteractions_->buildNeighborhood(SPHParticles_, virtualParticles_, minimumBoundary() - boundaryPadding, maximumBoundary() + boundaryPadding, supportRadius, supportRadius);
    SPHInteractions_->updatePriorForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->updatePressureForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
    updateSPHAcousticTimeStep(cpu::mode{});
    updateSPHAdvectionTimeStep();
    virtualParticleCoupling_.collectForceAndTorque(virtualParticles_, LSParticles());
}

void SPHDEM::initializeSPH(gpu::mode, cudaStream_t stream)
{
    initializeSPHState();
    initializeSPHDevice(stream);
}

void SPHDEM::initializeSPH(hybrid::mode, cudaStream_t stream)
{
    initializeSPHState();
    initializeSPHDevice(stream);
}

void SPHDEM::flushSPHToCurrentTime(cudaStream_t stream, bool commitCouplingImpulse)
{
    if (SPHStepState_.pendingDEMSteps_ == 0)
    {
        return;
    }

    const Real representedTime = SPHStepState_.pendingDEMSteps_ * timeStep();
    switch (mode())
    {
    case executionMode::CPU:
        flushSPH(representedTime, cpu::mode{});
        break;
    case executionMode::GPU:
        flushSPH(representedTime, gpu::mode{}, stream);
        break;
    case executionMode::Hybrid:
        flushSPH(representedTime, hybrid::mode{}, stream);
        break;
    }
    SPHStepState_.pendingDEMSteps_ = 0;
    // During this acoustic interval the solid used the old load for all but
    // the final DEM half kick. Match that impulse only after the kick finishes.
    // Output/completion projections must never modify the rigid-body state.
    if (commitCouplingImpulse)
        SPHStepState_.couplingImpulseTime_ = representedTime - 0.5 * timeStep();
}

void SPHDEM::beginOutputSnapshot()
{
    if (SPHStepState_.pendingDEMSteps_ == 0 || (SPHStepState_.resume_ && mode() == executionMode::CPU))
    {
        return;
    }
    const cudaStream_t stream = gpuStream();
    // After a device solve(), host results are current but device state is
    // deferred. Preserve those current host results when restoring the device.
    captureSPHState(SPHStepState_.resume_, stream);
    outputObservationSnapshot_ = true;
    flushSPHToCurrentTime(stream);
}

void SPHDEM::endOutputSnapshot()
{
    if (outputObservationSnapshot_)
    {
        outputObservationSnapshot_ = false;
        restoreSPHState(gpuStream());
    }
}

void SPHDEM::beginSolveCompletionSnapshot()
{
    if (SPHStepState_.pendingDEMSteps_ == 0)
    {
        return;
    }
    const cudaStream_t stream = gpuStream();
    captureSPHState(true, stream);
    flushSPHToCurrentTime(stream);
}

void SPHDEM::endSolveCompletionSnapshot()
{
    const bool hadSnapshot = outputSnapshot_.active_;
    if (mode() != executionMode::CPU)
    {
        restoreSPHState(gpuStream());
    }
    SPHStepState_.resume_ = hadSnapshot;
}

void SPHDEM::captureSPHState(bool preserveCurrentHost, cudaStream_t stream)
{
    if (outputSnapshot_.active_)
    {
        throw std::logic_error("An SPH state snapshot is already active.");
    }
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (mode() != executionMode::CPU)
    {
        SPHParticles_.copyDeviceToHostAsync(stream);
        virtualParticles_.copyDeviceToHostAsync(stream);
        host_device_detail::synchronize(stream);
    }
#else
    (void)stream;
#endif
    outputSnapshot_.SPHParticles_ = SPHParticles_.host();
    outputSnapshot_.virtualParticles_ = virtualParticles_.host();
    outputSnapshot_.stepState_ = SPHStepState_;
    outputSnapshot_.couplingState_ = virtualParticleCoupling_.captureState();
    outputSnapshot_.preserveCurrentHost_ = preserveCurrentHost;
    outputSnapshot_.active_ = true;
}

void SPHDEM::restoreSPHState(cudaStream_t stream)
{
    if (!outputSnapshot_.active_)
    {
        return;
    }

    SPHParticleContainer::host_container_type currentSPHParticles;
    virtualParticleContainer::host_container_type currentVirtualParticles;
    if (outputSnapshot_.preserveCurrentHost_)
    {
        currentSPHParticles = std::move(SPHParticles_.host());
        currentVirtualParticles = std::move(virtualParticles_.host());
    }
    SPHParticles_.host() = std::move(outputSnapshot_.SPHParticles_);
    virtualParticles_.host() = std::move(outputSnapshot_.virtualParticles_);
    SPHStepState_ = outputSnapshot_.stepState_;
    virtualParticleCoupling_.restoreState(std::move(outputSnapshot_.couplingState_));
    // A projected observation may have rebuilt neighbor lists or device grids.
    // Rebuild search structures against restored integration state on demand,
    // without resetting density, forces, or the saved advection phase.
    invalidateSPHNeighborhood();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (mode() != executionMode::CPU)
    {
        SPHParticles_.copyHostToDeviceAsync(stream);
        virtualParticles_.copyHostToDeviceAsync(stream);
        host_device_detail::synchronize(stream);
    }
#else
    (void)stream;
#endif
    if (outputSnapshot_.preserveCurrentHost_)
    {
        SPHParticles_.host() = std::move(currentSPHParticles);
        virtualParticles_.host() = std::move(currentVirtualParticles);
    }
    outputSnapshot_.active_ = false;
    outputSnapshot_.preserveCurrentHost_ = false;
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
    return SPHParticles_.deviceMemoryGB() + virtualParticles_.deviceMemoryGB() + SPHSpatialGrid_.deviceMemoryGB() + virtualParticleSpatialGrid_.deviceMemoryGB() + SPHNeighborhood_.deviceMemoryGB();
}

double SPHDEM::systemDeviceMemoryGB() const noexcept { return LSDEM::systemDeviceMemoryGB() + SPHDeviceMemoryGB(); }

void SPHDEM::calculateSPHForceAndTorque(Real DEMTimeStep, cpu::mode)
{
    if (DEMTimeStep > 0.0)
    {
        ++SPHStepState_.pendingDEMSteps_;
        if (SPHStepState_.pendingDEMSteps_ >= SPHStepState_.acousticStepInterval_)
        {
            flushSPHToCurrentTime(nullptr, true);
        }
    }
    virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles());
}

void SPHDEM::prepareSPHAdvectionStep(cpu::mode)
{
    const Real particleMass = SPHProperties_.particleMass();
    const cpu::SPHKinematicsStatistics boundaryStatistics = SPHInteractions_->kinematicsStatistics(virtualParticles_);
    if (boundaryStatistics.invalidValueCount_ > 0)
    {
        throw std::runtime_error("The SPH boundary contains a non-finite velocity or acceleration.");
    }

    const Real fluidVelocityScale = SPHStepState_.observedMaximumVelocity_ > SPHProperties_.maximumVelocity_ ? SPHStepState_.observedMaximumVelocity_ : SPHProperties_.maximumVelocity_;
    const Real timeHorizon = SPHStepState_.advectionTimeStep_ > 0.0 ? SPHStepState_.advectionTimeStep_ : SPHStepState_.acousticTimeStep_;
    const Real fluidSearchBuffer = 2.0 * fluidVelocityScale * timeHorizon + SPHStepState_.observedMaximumAcceleration_ * timeHorizon * timeHorizon;
    const Real boundarySearchBuffer = (fluidVelocityScale + boundaryStatistics.maximumVelocity_) * timeHorizon +
                                      0.5 * (SPHStepState_.observedMaximumAcceleration_ + boundaryStatistics.maximumAcceleration_) * timeHorizon * timeHorizon;
    const Real searchBuffer = fluidSearchBuffer > boundarySearchBuffer ? fluidSearchBuffer : boundarySearchBuffer;
    SPHStepState_.neighborSearchRadius_ = 2.0 * SPHProperties_.smoothingLength_ + searchBuffer;

    const cpu::SPHInteractionParameters parameters{SPHProperties_.smoothingLength_,
                                                   particleMass,
                                                   SPHProperties_.referenceDensity_,
                                                   SPHProperties_.latticeKernelSum3D_,
                                                   SPHProperties_.soundSpeed_,
                                                   SPHProperties_.dynamicViscosity_,
                                                   gravity()};
    invalidateSPHNeighborhood();
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updateFreeSurface(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->reinitializeDensity(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->updatePriorForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
}

void SPHDEM::applySPHJetVelocity(Real representedTime, cpu::mode)
{
    if (SPHStepState_.jetsCompleted_)
    {
        return;
    }

    auto& particleHost = SPHParticles_.host();
    for (const SPHJetConstraint& constraint : SPHJets_)
    {
        const bool active = representedTime < constraint.endTime_;
        const SPHJet& jet = constraint.value_;
        const Vec3 velocity = jet.velocity();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (constraint.particleCount_ >= cpu::parallelParticleThreshold)
#endif
        for (int localParticleIndex = 0; localParticleIndex < constraint.particleCount_; ++localParticleIndex)
        {
            SPHParticle& particle = particleHost[constraint.particleIndexBegin_ + localParticleIndex];
            const bool constrained = active && execution::isInsideSPHJetPipe(particle.position(), jet.outletCenter(), jet.direction(), jet.radius(), constraint.pipeLength_);
            particle.setConstrained(constrained);
            if (constrained)
            {
                particle.setVelocity(velocity);
            }
        }
    }
}

void SPHDEM::updateSPHJetCompletion() noexcept
{
    if (!SPHStepState_.jetsCompleted_ && SPHStepState_.representedTime_ >= SPHStepState_.latestSPHJetEndTime_)
    {
        SPHStepState_.jetsCompleted_ = true;
    }
}

void SPHDEM::ensureSPHNeighborhood(cpu::mode)
{
    const Real supportRadius = 2.0 * SPHProperties_.smoothingLength_;
    const Real skin = std::max(Real{0.0}, SPHStepState_.neighborSearchRadius_ - supportRadius);
    if (!SPHNeighborhood_.needsRebuild(SPHParticles_, virtualParticles_, skin))
        return;
    const Vec3 padding{supportRadius, supportRadius, supportRadius};
    SPHInteractions_->buildNeighborhood(SPHParticles_, virtualParticles_, minimumBoundary() - padding, maximumBoundary() + padding, supportRadius, SPHStepState_.neighborSearchRadius_);
    SPHNeighborhood_.captureHost(SPHParticles_, virtualParticles_);
}

void SPHDEM::advanceSPH(Real timeStep, cpu::mode)
{
    const Real jetConstraintTime = SPHStepState_.representedTime_;
    applySPHJetVelocity(jetConstraintTime, cpu::mode{});
    if (SPHStepState_.timeInAdvectionStep_ <= math::defaultTolerance)
    {
        prepareSPHAdvectionStep(cpu::mode{});
    }

    const Real halfTimeStep = 0.5 * timeStep;
    const Real particleMass = SPHProperties_.particleMass();
    const cpu::SPHInteractionParameters parameters{SPHProperties_.smoothingLength_,
                                                   particleMass,
                                                   SPHProperties_.referenceDensity_,
                                                   SPHProperties_.latticeKernelSum3D_,
                                                   SPHProperties_.soundSpeed_,
                                                   SPHProperties_.dynamicViscosity_,
                                                   gravity()};
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updateDensity(SPHParticles_, virtualParticles_, parameters, halfTimeStep);
    SPHInteractions_->integratePosition(SPHParticles_, halfTimeStep);
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updatePressureForceAndBoundaryForce(SPHParticles_, virtualParticles_, parameters);
    SPHInteractions_->integrateVelocity(SPHParticles_, gravity(), timeStep);
    // A jet active at the substep start remains active for the full substep, so
    // its end time is quantized upward by at most one acoustic substep.
    applySPHJetVelocity(jetConstraintTime, cpu::mode{});
    SPHInteractions_->integratePosition(SPHParticles_, halfTimeStep);
    ensureSPHNeighborhood(cpu::mode{});
    SPHInteractions_->updateDensity(SPHParticles_, virtualParticles_, parameters, halfTimeStep);

    SPHStepState_.representedTime_ += timeStep;
    if (!SPHStepState_.jetsCompleted_ && SPHStepState_.representedTime_ >= SPHStepState_.latestSPHJetEndTime_)
    {
        applySPHJetVelocity(SPHStepState_.representedTime_, cpu::mode{});
    }
    updateSPHJetCompletion();
    SPHStepState_.timeInAdvectionStep_ += timeStep;
    const bool advectionStepCompleted = SPHStepState_.timeInAdvectionStep_ + math::defaultTolerance >= SPHStepState_.advectionTimeStep_;
    if (advectionStepCompleted)
    {
        SPHStepState_.timeInAdvectionStep_ = 0.0;
    }
    updateSPHAcousticTimeStep(cpu::mode{});
    if (advectionStepCompleted)
    {
        updateSPHAdvectionTimeStep();
    }
}

void SPHDEM::flushSPH(Real representedTime, cpu::mode)
{
    if (!SPHInitializationState_.boundaryStatic_)
    {
        virtualParticleCoupling_.averageKinematics(virtualParticles_, LSParticles(), SPHStepState_.pendingDEMSteps_);
    }
    advanceSPH(representedTime, cpu::mode{});
    virtualParticleCoupling_.collectForceAndTorque(virtualParticles_, LSParticles());
    virtualParticleCoupling_.clearKinematics();
}

void SPHDEM::updateSPHAcousticTimeStep(cpu::mode)
{
    if (SPHParticles_.hostSize() == 0 || SPHProperties_.smoothingLength_ <= math::defaultTolerance)
    {
        SPHStepState_.acousticTimeStepLimit_ = timeStep();
        SPHStepState_.acousticTimeStep_ = timeStep();
        SPHStepState_.observedMaximumVelocity_ = 0.0;
        SPHStepState_.observedMaximumAcceleration_ = 0.0;
        SPHStepState_.acousticStepInterval_ = 1;
        return;
    }

    const cpu::SPHStateStatistics statistics = SPHInteractions_->stateStatistics(SPHParticles_, gravity());
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

    SPHStepState_.observedMaximumVelocity_ = statistics.maximumVelocity_;
    SPHStepState_.observedMaximumAcceleration_ = statistics.maximumAcceleration_;
    SPHStepState_.acousticTimeStepLimit_ = execution::calculateSPHAcousticTimeStep(SPHProperties_.smoothingLength_, SPHProperties_.soundSpeed_, SPHStepState_.observedMaximumVelocity_);
    if (!math::isFinite(SPHStepState_.acousticTimeStepLimit_) || SPHStepState_.acousticTimeStepLimit_ <= 0.0)
    {
        throw std::runtime_error("Failed to calculate a finite positive SPH acoustic time-step limit.");
    }

    Real maximumSPHTimeStep = SPHStepState_.acousticTimeStepLimit_;
    if (SPHStepState_.timeInAdvectionStep_ > math::defaultTolerance)
    {
        const Real remainingAdvectionTime = SPHStepState_.advectionTimeStep_ - SPHStepState_.timeInAdvectionStep_;
        maximumSPHTimeStep = remainingAdvectionTime < maximumSPHTimeStep ? remainingAdvectionTime : maximumSPHTimeStep;
    }
    setSPHAcousticTimeStepFromLimit(maximumSPHTimeStep);
}

void SPHDEM::updateSPHAdvectionTimeStep()
{
    if (SPHParticles_.hostSize() == 0 || SPHProperties_.smoothingLength_ <= math::defaultTolerance)
    {
        SPHStepState_.advectionTimeStepLimit_ = timeStep();
        SPHStepState_.advectionTimeStep_ = timeStep();
        SPHStepState_.timeInAdvectionStep_ = 0.0;
        SPHStepState_.advectionStepInterval_ = 1;
        return;
    }

    SPHStepState_.advectionTimeStepLimit_ = execution::calculateSPHAdvectionTimeStepFromLimits(SPHProperties_.smoothingLength_,
                                                                                               SPHProperties_.viscousTimeScale_,
                                                                                               SPHStepState_.observedMaximumVelocity_,
                                                                                               SPHStepState_.observedMaximumAcceleration_,
                                                                                               SPHProperties_.maximumVelocity_);
    if (!math::isFinite(SPHStepState_.advectionTimeStepLimit_) || SPHStepState_.advectionTimeStepLimit_ <= 0.0)
    {
        throw std::runtime_error("Failed to calculate a finite positive SPH advection time-step limit.");
    }

    const Real maximumSPHTimeStep = SPHStepState_.advectionTimeStepLimit_ < SPHStepState_.acousticTimeStepLimit_ ? SPHStepState_.advectionTimeStepLimit_ : SPHStepState_.acousticTimeStepLimit_;
    setSPHAcousticTimeStepFromLimit(maximumSPHTimeStep);
    const Real maximumAdvectionStepInterval = std::floor(SPHStepState_.advectionTimeStepLimit_ / SPHStepState_.acousticTimeStep_ + math::defaultTolerance);
    SPHStepState_.advectionStepInterval_ = static_cast<int>(std::clamp(maximumAdvectionStepInterval, Real{1.0}, Real{std::numeric_limits<int>::max()}));
    SPHStepState_.advectionTimeStep_ = SPHStepState_.advectionStepInterval_ * SPHStepState_.acousticTimeStep_;
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
    SPHStepState_.acousticStepInterval_ = static_cast<int>(std::clamp(maximumStepInterval, Real{1.0}, Real{std::numeric_limits<int>::max()}));
    SPHStepState_.acousticTimeStep_ = SPHStepState_.acousticStepInterval_ * DEMTimeStep;
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

void SPHDEM::applySPHJetVelocity(Real representedTime, gpu::mode, cudaStream_t stream)
{
    if (SPHStepState_.jetsCompleted_)
    {
        return;
    }

    for (const SPHJetConstraint& constraint : SPHJets_)
    {
        const bool active = representedTime < constraint.endTime_;
        const SPHJet& jet = constraint.value_;
        cuda::launchApplySPHJetVelocity(SPHParticles_,
                                        constraint.particleIndexBegin_,
                                        constraint.particleCount_,
                                        active,
                                        jet.outletCenter(),
                                        jet.direction(),
                                        jet.radius(),
                                        constraint.pipeLength_,
                                        jet.velocity(),
                                        stream);
    }
}

void SPHDEM::initializeSPHDevice(cudaStream_t stream)
{
    SPHParticles_.copyHostToDeviceAsync(stream);
    virtualParticles_.copyHostToDeviceAsync(stream);

    const Real supportRadius = SPHProperties_.smoothingLength_ > math::defaultTolerance ? 2.0 * SPHProperties_.smoothingLength_ : math::norm(maximumBoundary() - minimumBoundary());
    const Real particleMass = SPHProperties_.particleMass();
    const Vec3 boundaryPadding{supportRadius, supportRadius, supportRadius};
    SPHSpatialGrid_.setUniform(SPHParticles_, supportRadius, minimumBoundary() - boundaryPadding, maximumBoundary() + boundaryPadding, stream);
    virtualParticleSpatialGrid_.setUniform(virtualParticles_, supportRadius, minimumBoundary() - boundaryPadding, maximumBoundary() + boundaryPadding, stream);

    SPHStepState_.neighborSearchRadius_ = supportRadius;
    SPHStepState_.timeInAdvectionStep_ = 0.0;
    cuda::launchInitializeSPHParticles(SPHParticles_, SPHProperties_.referenceDensity_, SPHProperties_.soundSpeed_, stream);
    if (LSParticles().usesDevice())
    {
        cuda::launchUpdateVirtualParticlePositionAndNormal(virtualParticles_, LSParticles(), stream);
        cuda::launchClearVirtualParticleKinematics(virtualParticles_, stream);
        cuda::launchAccumulateVirtualParticleKinematics(virtualParticles_, LSParticles(), gravity(), stream);
    }
    cuda::launchBuildSPHSpatialGrid(SPHSpatialGrid_, SPHParticles_, stream);
    cuda::launchBuildVirtualParticleSpatialGrid(virtualParticleSpatialGrid_, virtualParticles_, stream);
    cuda::launchUpdateSPHPriorForceAndBoundaryForce(SPHParticles_,
                                                    virtualParticles_,
                                                    SPHSpatialGrid_.device(),
                                                    virtualParticleSpatialGrid_.device(),
                                                    SPHProperties_.smoothingLength_,
                                                    SPHStepState_.neighborSearchRadius_,
                                                    particleMass,
                                                    SPHProperties_.dynamicViscosity_,
                                                    stream);
    cuda::launchUpdateSPHPressureForceAndBoundaryForce(SPHParticles_,
                                                       virtualParticles_,
                                                       SPHSpatialGrid_.device(),
                                                       virtualParticleSpatialGrid_.device(),
                                                       SPHProperties_.smoothingLength_,
                                                       SPHStepState_.neighborSearchRadius_,
                                                       particleMass,
                                                       SPHProperties_.referenceDensity_,
                                                       SPHProperties_.soundSpeed_,
                                                       gravity(),
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

void SPHDEM::calculateSPHForceAndTorque(Real DEMTimeStep, gpu::mode, cudaStream_t stream)
{
    if (DEMTimeStep > 0.0)
    {
        ++SPHStepState_.pendingDEMSteps_;
        if (SPHStepState_.pendingDEMSteps_ >= SPHStepState_.acousticStepInterval_)
        {
            flushSPHToCurrentTime(stream, true);
        }
    }
    cuda::launchAddVirtualParticleForceAndTorque(mutableLSParticles(), virtualParticles_, stream);
}

void SPHDEM::calculateSPHForceAndTorque(Real DEMTimeStep, hybrid::mode, cudaStream_t stream)
{
    if (DEMTimeStep > 0.0)
    {
        ++SPHStepState_.pendingDEMSteps_;
        if (SPHStepState_.pendingDEMSteps_ >= SPHStepState_.acousticStepInterval_)
        {
            flushSPHToCurrentTime(stream, true);
        }
    }
    virtualParticleCoupling_.applyForceAndTorque(mutableLSParticles());
}

void SPHDEM::prepareSPHAdvectionStep(gpu::mode, cudaStream_t stream)
{
    const Real particleMass = SPHProperties_.particleMass();
    const cuda::SPHKinematicsStatistics boundaryStatistics = cuda::calculateSPHKinematicsStatistics(virtualParticles_, stream);
    if (boundaryStatistics.invalidValueCount_ > 0)
    {
        throw std::runtime_error("The SPH boundary contains a non-finite velocity or acceleration.");
    }

    const Real fluidVelocityScale = SPHStepState_.observedMaximumVelocity_ > SPHProperties_.maximumVelocity_ ? SPHStepState_.observedMaximumVelocity_ : SPHProperties_.maximumVelocity_;
    const Real timeHorizon = SPHStepState_.advectionTimeStep_ > 0.0 ? SPHStepState_.advectionTimeStep_ : SPHStepState_.acousticTimeStep_;
    const Real fluidSearchBuffer = 2.0 * fluidVelocityScale * timeHorizon + SPHStepState_.observedMaximumAcceleration_ * timeHorizon * timeHorizon;
    const Real boundarySearchBuffer = (fluidVelocityScale + boundaryStatistics.maximumVelocity_) * timeHorizon +
                                      0.5 * (SPHStepState_.observedMaximumAcceleration_ + boundaryStatistics.maximumAcceleration_) * timeHorizon * timeHorizon;
    const Real searchBuffer = fluidSearchBuffer > boundarySearchBuffer ? fluidSearchBuffer : boundarySearchBuffer;
    SPHStepState_.neighborSearchRadius_ = 2.0 * SPHProperties_.smoothingLength_ + searchBuffer;

    invalidateSPHNeighborhood();
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHFreeSurface(SPHParticles_,
                                     virtualParticles_,
                                     SPHSpatialGrid_.device(),
                                     virtualParticleSpatialGrid_.device(),
                                     SPHProperties_.smoothingLength_,
                                     SPHStepState_.neighborSearchRadius_,
                                     particleMass,
                                     stream);
    cuda::launchReinitializeSPHDensity(SPHParticles_,
                                       virtualParticles_,
                                       SPHSpatialGrid_.device(),
                                       virtualParticleSpatialGrid_.device(),
                                       SPHProperties_.smoothingLength_,
                                       SPHStepState_.neighborSearchRadius_,
                                       particleMass,
                                       SPHProperties_.referenceDensity_,
                                       SPHProperties_.latticeKernelSum3D_,
                                       SPHProperties_.soundSpeed_,
                                       stream);
    cuda::launchUpdateSPHPriorForceAndBoundaryForce(SPHParticles_,
                                                    virtualParticles_,
                                                    SPHSpatialGrid_.device(),
                                                    virtualParticleSpatialGrid_.device(),
                                                    SPHProperties_.smoothingLength_,
                                                    SPHStepState_.neighborSearchRadius_,
                                                    particleMass,
                                                    SPHProperties_.dynamicViscosity_,
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
    cuda::launchBuildSPHSpatialGrid(SPHSpatialGrid_, SPHParticles_, stream);
    cuda::launchBuildVirtualParticleSpatialGrid(virtualParticleSpatialGrid_, virtualParticles_, stream);
    SPHNeighborhood_.captureDevice(SPHParticles_, virtualParticles_, stream);
}

void SPHDEM::advanceSPH(Real timeStep, gpu::mode, cudaStream_t stream)
{
    const Real jetConstraintTime = SPHStepState_.representedTime_;
    applySPHJetVelocity(jetConstraintTime, gpu::mode{}, stream);
    if (SPHStepState_.timeInAdvectionStep_ <= math::defaultTolerance)
    {
        prepareSPHAdvectionStep(gpu::mode{}, stream);
    }

    const Real halfTimeStep = 0.5 * timeStep;
    const Real particleMass = SPHProperties_.particleMass();
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHDensity(SPHParticles_,
                                 virtualParticles_,
                                 SPHSpatialGrid_.device(),
                                 virtualParticleSpatialGrid_.device(),
                                 SPHProperties_.smoothingLength_,
                                 SPHStepState_.neighborSearchRadius_,
                                 particleMass,
                                 SPHProperties_.referenceDensity_,
                                 SPHProperties_.soundSpeed_,
                                 halfTimeStep,
                                 gravity(),
                                 stream);
    cuda::launchIntegrateSPHPosition(SPHParticles_, halfTimeStep, stream);
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHPressureForceAndBoundaryForce(SPHParticles_,
                                                       virtualParticles_,
                                                       SPHSpatialGrid_.device(),
                                                       virtualParticleSpatialGrid_.device(),
                                                       SPHProperties_.smoothingLength_,
                                                       SPHStepState_.neighborSearchRadius_,
                                                       particleMass,
                                                       SPHProperties_.referenceDensity_,
                                                       SPHProperties_.soundSpeed_,
                                                       gravity(),
                                                       stream);
    cuda::launchIntegrateSPHVelocity(SPHParticles_, gravity(), timeStep, stream);
    // Keep the boundary condition fixed through this complete acoustic substep.
    applySPHJetVelocity(jetConstraintTime, gpu::mode{}, stream);
    cuda::launchIntegrateSPHPosition(SPHParticles_, halfTimeStep, stream);
    ensureSPHNeighborhood(gpu::mode{}, stream);
    cuda::launchUpdateSPHDensity(SPHParticles_,
                                 virtualParticles_,
                                 SPHSpatialGrid_.device(),
                                 virtualParticleSpatialGrid_.device(),
                                 SPHProperties_.smoothingLength_,
                                 SPHStepState_.neighborSearchRadius_,
                                 particleMass,
                                 SPHProperties_.referenceDensity_,
                                 SPHProperties_.soundSpeed_,
                                 halfTimeStep,
                                 gravity(),
                                 stream);

    SPHStepState_.representedTime_ += timeStep;
    if (!SPHStepState_.jetsCompleted_ && SPHStepState_.representedTime_ >= SPHStepState_.latestSPHJetEndTime_)
    {
        applySPHJetVelocity(SPHStepState_.representedTime_, gpu::mode{}, stream);
    }
    updateSPHJetCompletion();
    SPHStepState_.timeInAdvectionStep_ += timeStep;
    const bool advectionStepCompleted = SPHStepState_.timeInAdvectionStep_ + math::defaultTolerance >= SPHStepState_.advectionTimeStep_;
    if (advectionStepCompleted)
    {
        SPHStepState_.timeInAdvectionStep_ = 0.0;
    }
    updateSPHAcousticTimeStep(gpu::mode{}, stream);
    if (advectionStepCompleted)
    {
        updateSPHAdvectionTimeStep();
    }
}

void SPHDEM::flushSPH(Real representedTime, gpu::mode, cudaStream_t stream)
{
    cuda::launchCacheVirtualParticleCouplingForce(virtualParticles_, stream);
    if (!SPHInitializationState_.boundaryStatic_)
    {
        cuda::launchAverageVirtualParticleKinematics(virtualParticles_, LSParticles(), 1.0 / Real(SPHStepState_.pendingDEMSteps_), stream);
    }
    advanceSPH(representedTime, gpu::mode{}, stream);
    cuda::launchRefreshVirtualParticleCouplingLoads(virtualParticles_, LSParticles(), stream);
    if (!SPHInitializationState_.boundaryStatic_)
    {
        cuda::launchClearVirtualParticleKinematics(virtualParticles_, stream);
    }
}

void SPHDEM::flushSPH(Real representedTime, hybrid::mode, cudaStream_t stream)
{
    if (!SPHInitializationState_.boundaryStatic_)
    {
        virtualParticleCoupling_.prepareDeviceKinematics(virtualParticles_, LSParticles(), SPHStepState_.pendingDEMSteps_, stream);
    }
    advanceSPH(representedTime, gpu::mode{}, stream);
    virtualParticleCoupling_.copyForceAndTorqueFromDevice(virtualParticles_, LSParticles(), stream);
    virtualParticleCoupling_.clearKinematics();
}

void SPHDEM::updateSPHAcousticTimeStep(gpu::mode, cudaStream_t stream)
{
    if (SPHParticles_.hostSize() == 0 || SPHProperties_.smoothingLength_ <= math::defaultTolerance)
    {
        SPHStepState_.acousticTimeStepLimit_ = timeStep();
        SPHStepState_.acousticTimeStep_ = timeStep();
        SPHStepState_.observedMaximumVelocity_ = 0.0;
        SPHStepState_.observedMaximumAcceleration_ = 0.0;
        SPHStepState_.acousticStepInterval_ = 1;
        return;
    }

    const cuda::SPHStateStatistics statistics = cuda::calculateSPHStateStatistics(SPHParticles_, gravity(), stream);
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

    SPHStepState_.observedMaximumVelocity_ = statistics.maximumVelocity_;
    SPHStepState_.observedMaximumAcceleration_ = statistics.maximumAcceleration_;
    SPHStepState_.acousticTimeStepLimit_ = execution::calculateSPHAcousticTimeStep(SPHProperties_.smoothingLength_, SPHProperties_.soundSpeed_, SPHStepState_.observedMaximumVelocity_);
    if (!math::isFinite(SPHStepState_.acousticTimeStepLimit_) || SPHStepState_.acousticTimeStepLimit_ <= 0.0)
    {
        throw std::runtime_error("Failed to calculate a finite positive SPH acoustic time-step limit.");
    }

    Real maximumSPHTimeStep = SPHStepState_.acousticTimeStepLimit_;
    if (SPHStepState_.timeInAdvectionStep_ > math::defaultTolerance)
    {
        const Real remainingAdvectionTime = SPHStepState_.advectionTimeStep_ - SPHStepState_.timeInAdvectionStep_;
        maximumSPHTimeStep = remainingAdvectionTime < maximumSPHTimeStep ? remainingAdvectionTime : maximumSPHTimeStep;
    }
    setSPHAcousticTimeStepFromLimit(maximumSPHTimeStep);
}

void SPHDEM::copySPHToHost(cudaStream_t stream)
{
    SPHParticles_.copyDeviceToHostAsync(stream);
    virtualParticles_.copyDeviceToHostAsync(stream);
}

void SPHDEM::copySystemDeviceToHost(cudaStream_t stream)
{
    LSDEM::copySystemDeviceToHost(stream);
    copySPHToHost(stream);
}

#else

void SPHDEM::initializeSPHDevice(cudaStream_t) { throw std::runtime_error("The SPH GPU backend is unavailable in this FunDEM build."); }
void SPHDEM::calculateSPHForceAndTorque(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::calculateSPHForceAndTorque(Real, hybrid::mode, cudaStream_t) {}
void SPHDEM::prepareSPHAdvectionStep(gpu::mode, cudaStream_t) {}
void SPHDEM::applySPHJetVelocity(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::advanceSPH(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::flushSPH(Real, gpu::mode, cudaStream_t) {}
void SPHDEM::flushSPH(Real, hybrid::mode, cudaStream_t) {}
void SPHDEM::updateSPHAcousticTimeStep(gpu::mode, cudaStream_t) {}
void SPHDEM::copySPHToHost(cudaStream_t) {}
void SPHDEM::copySystemDeviceToHost(cudaStream_t) {}

#endif

} // namespace fundem
