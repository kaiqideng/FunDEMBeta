#include "SphereDEM.h"

#include "data/energyOutput.h"
#include "execution/cpu/particleFunctions.h"

#include <stdexcept>
#include <utility>

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "data/HostAoSDeviceSoA.h"
#include "execution/cuda/bondForceKernel.cuh"
#include "execution/cuda/contactDetectionKernel.cuh"
#include "execution/cuda/contactForceKernel.cuh"
#include "execution/cuda/integrationKernel.cuh"
#include "execution/cuda/spatialGridKernel.cuh"
#endif

namespace fundem
{

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
namespace gpu
{

template <class ParticleStorage> void clearForceAndTorque(ParticleStorage& particles, cudaStream_t stream)
{
    const int particleCount = static_cast<int>(particles.deviceSize());
    if (particleCount <= 0)
    {
        return;
    }
    host_device_detail::checkCuda(cudaMemsetAsync(particles.template device<rigidBody::forceField>(), 0, particleCount * sizeof(math::Vec3), stream), "cudaMemsetAsync particle force");
    host_device_detail::checkCuda(cudaMemsetAsync(particles.template device<rigidBody::torqueField>(), 0, particleCount * sizeof(math::Vec3), stream), "cudaMemsetAsync particle torque");
}

} // namespace gpu
#endif

int SphereDEM::addSphere(const particle& value)
{
    if (initialized())
    {
        throw std::logic_error("Spheres cannot be added after solver initialization.");
    }
    if (!value.isValid())
    {
        throw std::invalid_argument("Cannot add an invalid sphere.");
    }
    if (!value.materialBelongsTo(materials()))
    {
        throw std::invalid_argument("The sphere must reference this solver's material container.");
    }
    if (value.getMaterial().isLevelSet())
    {
        throw std::invalid_argument("A sphere must use an ordinary material; LSMaterial is reserved for LSParticles.");
    }

    const int sphereIndex = static_cast<int>(spheres_.hostSize());
    spheres_.host().push_back(value);
    invalidateContinuation();
    return sphereIndex;
}

particle& SphereDEM::mutableSphere(int particleIndex)
{
    if (initialized())
        throw std::logic_error("Sphere motion cannot be changed while the solver is initialized.");
    if (particleIndex < 0 || particleIndex >= static_cast<int>(spheres_.hostSize()))
        throw std::out_of_range("The sphere index is outside this solver.");
    return spheres_.host()[particleIndex];
}

void SphereDEM::setSpherePosition(int particleIndex, const Vec3& value)
{
    particle& sphere = mutableSphere(particleIndex);
    if (!math::isFinite(value))
        throw std::invalid_argument("The sphere position must be finite.");
    if (value == sphere.position())
        return;
    sphere.setPosition(value);
    invalidateContinuation();
}

void SphereDEM::setSphereVelocity(int particleIndex, const Vec3& value)
{
    particle& sphere = mutableSphere(particleIndex);
    if (!math::isFinite(value))
        throw std::invalid_argument("The sphere velocity must be finite.");
    if (value == sphere.velocity())
        return;
    sphere.setVelocity(value);
    invalidateContinuation();
}

void SphereDEM::setSphereAngularVelocity(int particleIndex, const Vec3& value)
{
    particle& sphere = mutableSphere(particleIndex);
    if (!math::isFinite(value))
        throw std::invalid_argument("The sphere angular velocity must be finite.");
    if (value == sphere.angularVelocity())
        return;
    sphere.setAngularVelocity(value);
    invalidateContinuation();
}

int SphereDEM::addBond(const bond& value)
{
    validateBondForAddition(value);
    switch (value.type())
    {
    case bondType::sphereSphere:
    {
        if (!value.referencesParticleContainer(spheres_))
        {
            throw std::invalid_argument("The bond must reference this solver's sphere container.");
        }
        const int bondIndex = static_cast<int>(sphereInteractions_.bonds().hostSize());
        sphereInteractions_.bonds().host().push_back(value);
        invalidateContinuation();
        return bondIndex;
    }
    case bondType::sphereLSParticle:
    {
        if (!value.referencesParticleContainers(spheres_, LSParticles()))
        {
            throw std::invalid_argument("The bond must reference this solver's sphere and LSParticle containers.");
        }
        const int bondIndex = static_cast<int>(sphereLSInteractions_.bonds().hostSize());
        sphereLSInteractions_.bonds().host().push_back(value);
        invalidateContinuation();
        return bondIndex;
    }
    case bondType::LSParticleLSParticle:
        if (!value.referencesParticleContainer(LSParticles()))
        {
            throw std::invalid_argument("The bond must reference this solver's LSParticle container.");
        }
        return appendLSParticleBond(value);
    case bondType::undefined:
        throw std::invalid_argument("Cannot add a bond without a particle-pair type.");
    }
    throw std::invalid_argument("Unsupported bond particle-pair type.");
}

void SphereDEM::validateExecutionMode(executionMode value) const
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

void SphereDEM::executionModeChanged(executionMode value)
{
    spheres_.setUseDevice(value != executionMode::CPU);
    mutableLSParticles().setUseDevice(value == executionMode::GPU);
}

void SphereDEM::initialize(cpu::mode) { assembleForceAndTorque(0.0, cpu::mode{}); }

void SphereDEM::initialize(gpu::mode, cudaStream_t stream)
{
    initializeLSParticleDevice(stream);
    initializeSphereDevice(stream);
    assembleForceAndTorque(0.0, gpu::mode{}, stream);
}

void SphereDEM::initialize(hybrid::mode, cudaStream_t stream)
{
    initializeLSParticleDevice(stream);
    initializeSphereDevice(stream);
    mutableLSParticleSpatialGrid().set(LSParticles(), minimumBoundary(), maximumBoundary(), stream);
    assembleForceAndTorque(0.0, hybrid::mode{}, stream);
}

void SphereDEM::advanceCPU(Real timeStep)
{
    const Real halfTimeStep = 0.5 * timeStep;
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cpu::integrateVelocity(spheres_, gravity(), halfTimeStep);
    cpu::integratePosition(mutableLSParticles(), timeStep);
    cpu::integratePosition(spheres_, timeStep);
    assembleForceAndTorque(timeStep, cpu::mode{});
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cpu::integrateVelocity(spheres_, gravity(), halfTimeStep);
}

void SphereDEM::advanceGPU(Real timeStep)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const cudaStream_t stream = gpuStream();
    const Real halfTimeStep = 0.5 * timeStep;
    cuda::launchVelocityAndAngularVelocityIntegration(mutableLSParticles(), gravity(), halfTimeStep, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(spheres_, gravity(), halfTimeStep, stream);
    cuda::launchPositionAndOrientationIntegration(mutableLSParticles(), timeStep, stream);
    cuda::launchPositionAndOrientationIntegration(spheres_, timeStep, stream);
    assembleForceAndTorque(timeStep, gpu::mode{}, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(mutableLSParticles(), gravity(), halfTimeStep, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(spheres_, gravity(), halfTimeStep, stream);
#else
    (void)timeStep;
#endif
}

void SphereDEM::advanceHybrid(Real timeStep)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const cudaStream_t stream = gpuStream();
    const Real halfTimeStep = 0.5 * timeStep;
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cuda::launchVelocityAndAngularVelocityIntegration(spheres_, gravity(), halfTimeStep, stream);
    cpu::integratePosition(mutableLSParticles(), timeStep);
    cuda::launchPositionAndOrientationIntegration(spheres_, timeStep, stream);
    assembleForceAndTorque(timeStep, hybrid::mode{}, stream);
    cpu::integrateVelocity(mutableLSParticles(), gravity(), halfTimeStep);
    cuda::launchVelocityAndAngularVelocityIntegration(spheres_, gravity(), halfTimeStep, stream);
#else
    (void)timeStep;
#endif
}

void SphereDEM::setSphereVTUFields(std::vector<sphereVTUField> fields) { sphereVTUFields_ = std::move(fields); }

energyRecord SphereDEM::systemEnergy() const
{
    energyRecord result = LSDEM::systemEnergy();
    addParticleEnergy(result, spheres_, gravity());
    addContactEnergy(result, sphereInteractions_.contacts());
    addContactEnergy(result, sphereLSInteractions_.contacts());
    addBondEnergy(result, sphereInteractions_.bonds());
    addBondEnergy(result, sphereLSInteractions_.bonds());
    return result;
}

void SphereDEM::writeSystemVTU(int frameIndex)
{
    LSDEM::writeSystemVTU(frameIndex);
    writeSphereVTU(frameIndex);
}

void SphereDEM::writeSphereVTU(int frameIndex)
{
    if (spheres_.hostSize() > 0)
    {
        makeSphereVTUWriter(makeVTUFileName("particle", "sphere", frameIndex), spheres_, gravity(), sphereVTUFields_).write(vtuOutputFormat());
    }
    if (sphereInteractions_.contacts().hostSize() > 0)
    {
        makeContactVTUWriter(makeVTUFileName("interaction", "sphereContact", frameIndex), sphereInteractions_.contacts(), contactVTUFields()).write(vtuOutputFormat());
    }
    if (sphereLSInteractions_.contacts().hostSize() > 0)
    {
        makeContactVTUWriter(makeVTUFileName("interaction", "sphereLSParticleContact", frameIndex), sphereLSInteractions_.contacts(), contactVTUFields()).write(vtuOutputFormat());
    }
    if (sphereInteractions_.bonds().hostSize() > 0)
    {
        makeBondVTUWriter(makeVTUFileName("interaction", "sphereBond", frameIndex), sphereInteractions_.bonds(), spheres_, bondVTUFields()).write(vtuOutputFormat());
    }
    if (sphereLSInteractions_.bonds().hostSize() > 0)
    {
        makeBondVTUWriter(makeVTUFileName("interaction", "sphereLSParticleBond", frameIndex), sphereLSInteractions_.bonds(), spheres_, LSParticles(), bondVTUFields()).write(vtuOutputFormat());
    }
}

double SphereDEM::sphereDeviceMemoryGB() const noexcept
{
    return spheres_.deviceMemoryGB() + sphereInteractions_.deviceMemoryGB() + sphereLSInteractions_.deviceMemoryGB() + sphereSpatialGrid_.deviceMemoryGB();
}

double SphereDEM::systemDeviceMemoryGB() const noexcept { return LSDEM::systemDeviceMemoryGB() + sphereDeviceMemoryGB(); }

void SphereDEM::assembleForceAndTorque(Real historyTimeStep, cpu::mode)
{
    cpu::clearForceAndTorque(mutableLSParticles());
    cpu::clearForceAndTorque(spheres_);
    calculateLSParticleForceAndTorque(historyTimeStep, cpu::mode{});
    calculateSphereForceAndTorque(historyTimeStep, cpu::mode{});
    addLSParticleExternalForceAndTorque(mutableLSParticles().host());
    addSphereExternalForceAndTorque(spheres_.host());
}

void SphereDEM::assembleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    gpu::clearForceAndTorque(mutableLSParticles(), stream);
    gpu::clearForceAndTorque(spheres_, stream);
    cuda::launchBuildSpatialGrid(mutableLSParticleSpatialGrid(), LSParticles(), stream);
    calculateLSParticleForceAndTorque(historyTimeStep, gpu::mode{}, stream);
    calculateSphereForceAndTorque(historyTimeStep, gpu::mode{}, stream);
    addLSParticleExternalForceAndTorque(deviceFields(mutableLSParticles()), stream);
    addSphereExternalForceAndTorque(deviceFields(spheres_), stream);
#else
    (void)historyTimeStep;
    (void)stream;
#endif
}

void SphereDEM::assembleForceAndTorque(Real historyTimeStep, hybrid::mode, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    cpu::clearForceAndTorque(mutableLSParticles());
    gpu::clearForceAndTorque(spheres_, stream);
    calculateLSParticleForceAndTorque(historyTimeStep, cpu::mode{});
    calculateSphereForceAndTorque(historyTimeStep, hybrid::mode{}, stream);
    copyLSParticleForceAndTorqueToHost(stream);
    addLSParticleExternalForceAndTorque(mutableLSParticles().host());
    addSphereExternalForceAndTorque(deviceFields(spheres_), stream);
#else
    (void)historyTimeStep;
    (void)stream;
#endif
}

void SphereDEM::calculateSphereForceAndTorque(Real historyTimeStep, cpu::mode)
{
    sphereSearch_.findContacts(sphereInteractions_.contacts(), spheres_);
    cpu::addContactForceAndTorque(sphereInteractions_.contacts(), spheres_, spheres_, historyTimeStep);
    cpu::addBondForceAndTorque(sphereInteractions_.bonds(), spheres_);
    sphereLSSearch_.findContacts(sphereLSInteractions_.contacts(), spheres_, mutableLSParticles());
    cpu::addContactForceAndTorque(sphereLSInteractions_.contacts(), spheres_, mutableLSParticles(), historyTimeStep);
    cpu::addBondForceAndTorque(sphereLSInteractions_.bonds(), spheres_, mutableLSParticles());
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

void SphereDEM::initializeSphereDevice(cudaStream_t stream)
{
    spheres_.copyHostToDeviceAsync(stream);
    sphereInteractions_.initializeDevice(spheres_, stream);
    sphereSpatialGrid_.set(spheres_, minimumBoundary(), maximumBoundary(), stream);
    sphereLSInteractions_.initializeDevice(spheres_, stream);
}

void SphereDEM::calculateSphereForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream)
{
    cuda::launchBuildSpatialGrid(sphereSpatialGrid_, spheres_, stream);
    cuda::launchSphereSphereContactDetection(sphereInteractions_, spheres_, sphereSpatialGrid_.device(), stream);
    cuda::launchContactForce(sphereInteractions_.contacts(), spheres_, materials(), historyTimeStep, stream);
    cuda::launchBondForce(sphereInteractions_.bonds(), spheres_, stream);
    cuda::launchSphereLevelSetContactDetection(sphereLSInteractions_, spheres_, LSParticles(), geometries(), mutableLSParticleSpatialGrid().device(), stream);
    cuda::launchContactForce(sphereLSInteractions_.contacts(), spheres_, mutableLSParticles(), materials(), historyTimeStep, stream);
    cuda::launchBondForce(sphereLSInteractions_.bonds(), spheres_, mutableLSParticles(), stream);
}

void SphereDEM::calculateSphereForceAndTorque(Real historyTimeStep, hybrid::mode, cudaStream_t stream)
{
    mutableLSParticles().copyHostToDeviceAsync(stream);
    cuda::launchBuildSpatialGrid(mutableLSParticleSpatialGrid(), LSParticles(), stream);
    calculateSphereForceAndTorque(historyTimeStep, gpu::mode{}, stream);
}

void SphereDEM::copySphereDeviceToHost(cudaStream_t stream)
{
    spheres_.copyDeviceToHostAsync(stream);
    sphereInteractions_.copyDeviceToHostAsync(stream);
    sphereLSInteractions_.copyDeviceToHostAsync(stream);
}

void SphereDEM::copySystemDeviceToHost(cudaStream_t stream)
{
    LSDEM::copySystemDeviceToHost(stream);
    copySphereDeviceToHost(stream);
}

#else

void SphereDEM::initializeSphereDevice(cudaStream_t) { throw std::runtime_error("The SphereDEM GPU backend is unavailable in this FunDEM build."); }
void SphereDEM::calculateSphereForceAndTorque(Real, gpu::mode, cudaStream_t) {}
void SphereDEM::calculateSphereForceAndTorque(Real, hybrid::mode, cudaStream_t) {}
void SphereDEM::copySphereDeviceToHost(cudaStream_t) {}
void SphereDEM::copySystemDeviceToHost(cudaStream_t) {}

#endif

} // namespace fundem
