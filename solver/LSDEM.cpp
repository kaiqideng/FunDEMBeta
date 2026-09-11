#include "LSDEM.h"

#include "data/energyOutput.h"
#include "data/myLSObject.h"
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

int LSDEM::addMaterial(const LSMaterial& value) { return appendMaterial(value); }

int LSDEM::appendMaterial(const material& value)
{
    if (initialized())
    {
        throw std::logic_error("Materials cannot be added after solver initialization.");
    }
    if (!value.isValid())
    {
        throw std::invalid_argument("Cannot add an invalid material.");
    }

    const int materialIndex = static_cast<int>(materials_.hostSize());
    materials_.host().push_back(value);
    invalidateContinuation();
    return materialIndex;
}

int LSDEM::addGeometry(const levelset::LSInfo& value, bool isFixed)
{
    if (initialized())
    {
        throw std::logic_error("Geometries cannot be added after solver initialization.");
    }
    if (!value.valid())
    {
        throw std::invalid_argument("Cannot add an invalid level-set geometry.");
    }

    const int geometryIndex =
        geometries_.add(value.gridNodeOrigin(), value.gridNodeSpacing(), value.gridNodeSize3D(), value.gridNodeSFD(), value.surfaceNodePosition(), value.surfaceNodeConnectivity(), isFixed);
    invalidateContinuation();
    return geometryIndex;
}

int LSDEM::addLSParticle(const LSParticle& value)
{
    if (initialized())
    {
        throw std::logic_error("LSParticles cannot be added after solver initialization.");
    }
    if (!value.isValid())
    {
        throw std::invalid_argument("Cannot add an invalid LSParticle.");
    }
    if (!value.materialBelongsTo(materials_) || !value.geometryBelongsTo(geometries_))
    {
        throw std::invalid_argument("The LSParticle must reference this solver's material and geometry containers.");
    }

    const int particleIndex = static_cast<int>(LSParticles_.hostSize());
    LSParticles_.host().push_back(value);
    invalidateContinuation();
    return particleIndex;
}

LSParticle& LSDEM::mutableLSParticle(int particleIndex)
{
    if (initialized())
    {
        throw std::logic_error("LSParticle motion cannot be changed while the solver is initialized.");
    }
    if (particleIndex < 0 || particleIndex >= static_cast<int>(LSParticles_.hostSize()))
    {
        throw std::out_of_range("The LSParticle index is outside this solver.");
    }
    return LSParticles_.host()[particleIndex];
}

void LSDEM::setLSParticlePosition(int particleIndex, const Vec3& value)
{
    LSParticle& particle = mutableLSParticle(particleIndex);
    if (!math::isFinite(value))
    {
        throw std::invalid_argument("The LSParticle position must be finite.");
    }
    if (value == particle.position())
    {
        return;
    }
    particle.setPosition(value);
    invalidateContinuation();
}

void LSDEM::setLSParticleVelocity(int particleIndex, const Vec3& value)
{
    LSParticle& particle = mutableLSParticle(particleIndex);
    if (!math::isFinite(value))
    {
        throw std::invalid_argument("The LSParticle velocity must be finite.");
    }
    if (value == particle.velocity())
    {
        return;
    }
    particle.setVelocity(value);
    invalidateContinuation();
}

void LSDEM::setLSParticleAngularVelocity(int particleIndex, const Vec3& value)
{
    LSParticle& particle = mutableLSParticle(particleIndex);
    if (!math::isFinite(value))
    {
        throw std::invalid_argument("The LSParticle angular velocity must be finite.");
    }
    if (value == particle.angularVelocity())
    {
        return;
    }
    particle.setAngularVelocity(value);
    invalidateContinuation();
}

int LSDEM::addBond(const bond& value)
{
    validateBondForAddition(value);
    if (value.type() != bondType::LSParticleLSParticle)
    {
        throw std::invalid_argument("LSDEM only supports bonds between LSParticles.");
    }
    if (!value.referencesParticleContainer(LSParticles_))
    {
        throw std::invalid_argument("The bond must reference this solver's LSParticle container.");
    }
    return appendLSParticleBond(value);
}

void LSDEM::validateBondForAddition(const bond& value) const
{
    if (initialized())
    {
        throw std::logic_error("Bonds cannot be added after solver initialization.");
    }
    if (!value.isValid())
    {
        throw std::invalid_argument("Cannot add an invalid bond.");
    }
}

int LSDEM::appendLSParticleBond(const bond& value)
{
    bondContainer& bonds = LSParticleInteractions_.bonds();
    const int bondIndex = static_cast<int>(bonds.hostSize());
    bonds.host().push_back(value);
    invalidateContinuation();
    return bondIndex;
}

void LSDEM::setLSParticleVTUFields(std::vector<LSParticleVTUField> fields) { LSParticleVTUFields_ = std::move(fields); }

void LSDEM::setContactVTUFields(std::vector<contactVTUField> fields) { contactVTUFields_ = std::move(fields); }

void LSDEM::setBondVTUFields(std::vector<bondVTUField> fields) { bondVTUFields_ = std::move(fields); }

energyRecord LSDEM::systemEnergy() const
{
    energyRecord result;
    addParticleEnergy(result, LSParticles_, gravity());
    addContactEnergy(result, LSParticleInteractions_.contacts());
    addBondEnergy(result, LSParticleInteractions_.bonds());
    return result;
}

void LSDEM::writeSystemVTU(int frameIndex)
{
    bool haveLSParticles = false;
    bool haveFixedLSParticles = false;
    for (const LSParticle& value : LSParticles_.host())
    {
        const bool fixed = value.inverseMass() == 0.0;
        haveLSParticles = haveLSParticles || !fixed;
        haveFixedLSParticles = haveFixedLSParticles || fixed;
    }
    if (haveLSParticles)
    {
        makeLSParticleVTUWriter(makeVTUFileName("particle", "LSParticle", frameIndex), LSParticles_, gravity(), LSParticleVTUFields_).write(vtuOutputFormat());
    }
    if (haveFixedLSParticles)
    {
        makeFixedLSParticleVTUWriter(makeVTUFileName("particle", "fixedLSParticle", frameIndex), LSParticles_, gravity(), LSParticleVTUFields_).write(vtuOutputFormat());
    }
    if (LSParticleInteractions_.contacts().hostSize() > 0)
    {
        makeContactVTUWriter(makeVTUFileName("interaction", "LSParticleContact", frameIndex), LSParticleInteractions_.contacts(), contactVTUFields_).write(vtuOutputFormat());
    }
    if (LSParticleInteractions_.bonds().hostSize() > 0)
    {
        makeBondVTUWriter(makeVTUFileName("interaction", "LSParticleBond", frameIndex), LSParticleInteractions_.bonds(), LSParticles_, bondVTUFields_).write(vtuOutputFormat());
    }
}

void LSDEM::validateExecutionMode(executionMode value) const
{
    switch (value)
    {
    case executionMode::CPU:
    case executionMode::GPU:
        return;
    case executionMode::Hybrid:
        throw std::logic_error("LSDEM contains only LSParticles and does not define a Hybrid mode.");
    }
    throw std::invalid_argument("Invalid execution mode.");
}

void LSDEM::executionModeChanged(executionMode value) { LSParticles_.setUseDevice(value == executionMode::GPU); }

void LSDEM::initialize(cpu::mode) { assembleLSParticleForceAndTorque(0.0, cpu::mode{}); }

void LSDEM::initialize(gpu::mode, cudaStream_t stream)
{
    initializeLSParticleDevice(stream);
    assembleLSParticleForceAndTorque(0.0, gpu::mode{}, stream);
}

void LSDEM::initialize(hybrid::mode, cudaStream_t) { throw std::logic_error("LSDEM does not define a Hybrid execution template."); }

void LSDEM::advanceCPU(Real timeStep)
{
    const Real halfTimeStep = 0.5 * timeStep;
    cpu::integrateVelocity(LSParticles_, gravity(), halfTimeStep);
    cpu::integratePosition(LSParticles_, timeStep);
    assembleLSParticleForceAndTorque(timeStep, cpu::mode{});
    cpu::integrateVelocity(LSParticles_, gravity(), halfTimeStep);
}

void LSDEM::advanceGPU(Real timeStep)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const cudaStream_t stream = gpuStream();
    const Real halfTimeStep = 0.5 * timeStep;
    cuda::launchVelocityAndAngularVelocityIntegration(LSParticles_, gravity(), halfTimeStep, stream);
    cuda::launchPositionAndOrientationIntegration(LSParticles_, timeStep, stream);
    assembleLSParticleForceAndTorque(timeStep, gpu::mode{}, stream);
    cuda::launchVelocityAndAngularVelocityIntegration(LSParticles_, gravity(), halfTimeStep, stream);
#else
    (void)timeStep;
#endif
}

void LSDEM::advanceHybrid(Real) { throw std::logic_error("LSDEM does not define a Hybrid execution template."); }

void LSDEM::assembleLSParticleForceAndTorque(Real historyTimeStep, cpu::mode)
{
    cpu::clearForceAndTorque(LSParticles_);
    calculateLSParticleForceAndTorque(historyTimeStep, cpu::mode{});
    addLSParticleExternalForceAndTorque(LSParticles_.host());
}

void LSDEM::assembleLSParticleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const int particleCount = static_cast<int>(LSParticles_.deviceSize());
    if (particleCount > 0)
    {
        host_device_detail::checkCuda(cudaMemsetAsync(LSParticles_.device<rigidBody::forceField>(), 0, particleCount * sizeof(math::Vec3), stream), "cudaMemsetAsync LSParticle force");
        host_device_detail::checkCuda(cudaMemsetAsync(LSParticles_.device<rigidBody::torqueField>(), 0, particleCount * sizeof(math::Vec3), stream), "cudaMemsetAsync LSParticle torque");
    }
    cuda::launchBuildSpatialGrid(LSParticleSpatialGrid_, LSParticles_, stream);
    calculateLSParticleForceAndTorque(historyTimeStep, gpu::mode{}, stream);
    addLSParticleExternalForceAndTorque(deviceFields(LSParticles_), stream);
#else
    (void)historyTimeStep;
    (void)stream;
#endif
}

void LSDEM::calculateLSParticleForceAndTorque(Real historyTimeStep, cpu::mode)
{
    LSParticleSearch_.findContacts(LSParticleInteractions_.contacts(), LSParticles_);
    cpu::addContactForceAndTorque(LSParticleInteractions_.contacts(), LSParticles_, LSParticles_, historyTimeStep);
    cpu::addBondForceAndTorque(LSParticleInteractions_.bonds(), LSParticles_);
}

double LSDEM::systemDeviceMemoryGB() const noexcept
{
    return materials_.deviceMemoryGB() + geometries_.deviceMemoryGB() + LSParticles_.deviceMemoryGB() + LSParticleInteractions_.deviceMemoryGB() + LSParticleSpatialGrid_.deviceMemoryGB();
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

void LSDEM::initializeLSParticleDevice(cudaStream_t stream)
{
    materials_.copyHostToDeviceAsync(stream);
    geometries_.copyHostToDeviceAsync(stream);
    LSParticles_.copyHostToDeviceAsync(stream);
    if (LSParticles_.usesDevice())
    {
        LSParticleInteractions_.initializeDevice(LSParticles_, geometries_, stream);
        LSParticleSpatialGrid_.set(LSParticles_, minimumBoundary(), maximumBoundary(), stream);
    }
}

void LSDEM::calculateLSParticleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream)
{
    cuda::launchLevelSetLevelSetContactDetection(LSParticleInteractions_, LSParticles_, geometries_, LSParticleSpatialGrid_.device(), stream);
    cuda::launchContactForce(LSParticleInteractions_.contacts(), LSParticles_, materials_, historyTimeStep, stream);
    cuda::launchBondForce(LSParticleInteractions_.bonds(), LSParticles_, stream);
}

void LSDEM::copyLSParticleForceAndTorqueToHost(cudaStream_t stream)
{
    LSParticles_.copyDeviceFieldToHostAsync<rigidBody::forceField>(stream);
    LSParticles_.copyDeviceFieldToHostAsync<rigidBody::torqueField>(stream);
    host_device_detail::synchronize(stream);
}

void LSDEM::copySystemDeviceToHost(cudaStream_t stream)
{
    if (LSParticles_.usesDevice())
    {
        LSParticles_.copyDeviceToHostAsync(stream);
        LSParticleInteractions_.copyDeviceToHostAsync(stream);
    }
}

#else

void LSDEM::initializeLSParticleDevice(cudaStream_t) { throw std::runtime_error("The GPU backend is unavailable in this FunDEM build."); }
void LSDEM::calculateLSParticleForceAndTorque(Real, gpu::mode, cudaStream_t) {}
void LSDEM::copyLSParticleForceAndTorqueToHost(cudaStream_t) {}
void LSDEM::copySystemDeviceToHost(cudaStream_t) {}

#endif

} // namespace fundem
