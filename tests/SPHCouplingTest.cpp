#include "data/myLSObject.h"
#include "solver/SPHDEM.h"
#include "interaction/virtualParticleCoupling.h"

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "execution/cuda/sphInteractionKernel.cuh"
#endif

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using namespace fundem;
using math::Real;
using math::Vec3;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void requireNear(const Vec3& value, const Vec3& expected, const char* message, Real tolerance = 1.0e-10)
{
    if (math::norm(value - expected) > tolerance * (1.0 + math::norm(expected)))
    {
        std::cerr << message << ": got " << value.x << ',' << value.y << ',' << value.z << "; expected " << expected.x << ',' << expected.y << ',' << expected.z << '\n';
        throw std::runtime_error(message);
    }
}

LSParticle syntheticOwner()
{
    LSParticle owner;
    static_cast<rigidBody&>(owner).setMassAndInertia(2.0, math::Mat3::diagonal({3.0, 3.0, 3.0}));
    return owner;
}

void testHeldWorldTorqueAndImpulse()
{
    LSParticleContainer owners;
    owners.host().push_back(syntheticOwner());
    virtualParticleContainer samples;
    samples.host().emplace_back(Vec3::unitX(), Vec3::unitX(), 1.0, 0);
    virtualParticleCoupling coupling;
    coupling.initialize(samples, owners, Vec3::zero());
    samples.host()[0].setForce({0.0, 1.0, 0.0});
    coupling.collectForceAndTorque(samples, owners);
    const auto oldState = coupling.captureState();

    LSParticle& owner = owners.host()[0];
    const Real halfAngle = math::pi / 8.0;
    owner.setOrientation({std::cos(halfAngle), 0.0, 0.0, std::sin(halfAngle)});
    coupling.applyForceAndTorque(owners);
    requireNear(owner.torque(), {0.0, 0.0, 1.0}, "A rotated owner must retain the last world torque until fluid refresh");

    constexpr int stepCount = 5;
    constexpr Real step = 0.01;
    constexpr Real interval = stepCount * step;
    const Vec3 externalForce{0.3, 0.0, 0.0};
    const Vec3 oldForce = owner.force();
    const Vec3 oldTorque = owner.torque();
    // Actual DEM half-kicks with held loads: all but the final half kick use F0.
    owner.setForce(oldForce + externalForce);
    owner.updateVelocityAndAngularVelocity(Vec3::zero(), interval - 0.5 * step);
    samples.host()[0].setPosition(owner.position() + math::rotateUnit(owner.orientation(), samples.host()[0].localPosition()));
    samples.host()[0].setForce({0.0, 4.0, 0.0});
    coupling.collectForceAndTorque(samples, owners);
    const Vec3 newTorque = math::cross(samples.host()[0].position() - owner.position(), samples.host()[0].force());
    owner.setForce(externalForce);
    owner.setTorque(Vec3::zero());
    coupling.applyForceAndTorque(owners);
    owner.updateVelocityAndAngularVelocity(Vec3::zero(), 0.5 * step);
    coupling.applyImpulseCorrection(owners, interval - 0.5 * step);
    requireNear(owner.velocity() / owner.inverseMass(),
                interval * (samples.host()[0].force() + externalForce),
                "Transient wall impulse must match the fluid impulse without duplicating external impulse");
    requireNear(3.0 * owner.angularVelocity(), interval * newTorque, "Centroidal angular impulse must match the refreshed torque for isotropic inertia");
    require(math::norm(newTorque - oldTorque) > 0.1, "Impulse test must change the held torque");

    coupling.restoreState(oldState);
    owner.setForce(Vec3::zero());
    owner.setTorque(Vec3::zero());
    coupling.applyForceAndTorque(owners);
    requireNear(owner.force(), oldForce, "Restoring coupling state must restore the held force");
    requireNear(owner.torque(), oldTorque, "Restoring coupling state must restore the held torque");
}

class forcedBoundarySolver : public SPHDEM
{
public:
    explicit forcedBoundarySolver(executionMode mode = executionMode::CPU) : SPHDEM(mode) { setSPHBoundaryMotionExpected(true); }
    Vec3 externalForce_{Vec3::zero()};
    void prescribeVelocity(const Vec3& velocity)
    {
        mutableLSParticles().host()[0].setVelocity(velocity);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        if (LSParticles().usesDevice())
            mutableLSParticles().copyHostFieldToDeviceAsync<pointMass::velocityField>(gpuStream());
#endif
    }

protected:
    void addLSParticleExternalForceAndTorque(LSParticleContainer::host_container_type& owners) override
    {
        for (LSParticle& owner : owners)
            owner.addForce(externalForce_);
    }
};

void configureBoundary(SPHDEM& simulation, bool nearbyFluid, bool rotating = false, bool infinite = false)
{
    constexpr Real spacing = 0.02;
    simulation.setBoundary({-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    simulation.setGravity(Vec3::zero());
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(spacing, 1.3 * spacing, 1000.0, 0.0);
    simulation.setSPHMaximumVelocity(1.0);
    simulation.setSPHSoundSpeed(10.0);
    simulation.addSPHParticle(SPHParticle{nearbyFluid ? Vec3{0.0, 0.0, 0.055} : Vec3{0.2, 0.2, 0.2}, nearbyFluid ? Vec3{0.0, 0.0, -0.1} : Vec3::zero()});
    const int material = simulation.addMaterial(LSMaterial{0.0, 0.0, 0.0, 1.0, 1000.0});
    levelset::BoxParticle box{{0.08, 0.08, 0.08}};
    box.buildLSGrid(spacing, 3);
    const int geometry = simulation.addGeometry(box);
    LSParticle owner;
    owner.setMaterial(simulation.materials(), material);
    owner.setGeometry(simulation.geometries(), geometry);
    if (rotating)
        owner.setAngularVelocity({0.0, 0.0, 2.0});
    if (infinite)
        owner.setInfiniteMass();
    simulation.addLSParticle(owner);
    simulation.initialize();
    require(infinite == (simulation.LSParticles().host()[0].inverseMass() == 0.0), "Coupling test boundary mass classification is incorrect");
    require(simulation.SPHStepInterval() > 2, "Coupling test must exercise deferred acoustic updates");
}

Vec3 momentum(const SPHDEM& simulation)
{
    Vec3 result = Vec3::zero();
    for (const SPHParticle& particle : simulation.SPHParticles().host())
        result += particle.mass() * particle.velocity();
    for (const LSParticle& owner : simulation.LSParticles().host())
        result += owner.velocity() / owner.inverseMass();
    return result;
}

void testCompleteEndpointAcceleration()
{
    forcedBoundarySolver simulation;
    simulation.externalForce_ = {0.04, 0.0, 0.0};
    configureBoundary(simulation, false);
    const Vec3 expected = simulation.LSParticles().host()[0].inverseMass() * simulation.externalForce_;
    simulation.step();
    const Vec3 velocityBefore = simulation.LSParticles().host()[0].velocity();
    simulation.observeCurrentState(
        [&](const solver&)
        {
            require(simulation.virtualParticles().hostSize() > 0, "Acceleration test requires virtual boundary samples");
            for (const virtualParticle& sample : simulation.virtualParticles().host())
                requireNear(sample.acceleration(), expected, "Boundary acceleration must include the complete previous endpoint external force");
        });
    requireNear(simulation.LSParticles().host()[0].velocity(), velocityBefore, "Observation must not apply a rigid-body impulse");
}

void testPrescribedAccelerationHistory(executionMode mode)
{
    forcedBoundarySolver simulation{mode};
    configureBoundary(simulation, false, false, true);
    simulation.step();
    simulation.observeCurrentState(
        [&](const solver&)
        {
            for (const virtualParticle& sample : simulation.virtualParticles().host())
                requireNear(sample.acceleration(), Vec3::zero(), "Constant prescribed velocity must have zero sampled acceleration");
        });
    simulation.prescribeVelocity({2.0 * simulation.timeStep(), 0.0, 0.0});
    simulation.step();
    simulation.observeCurrentState(
        [&](const solver&)
        {
            for (const virtualParticle& sample : simulation.virtualParticles().host())
            {
                if (mode == executionMode::GPU)
                    // Device velocity/acceleration are cleared accumulator storage
                    // after a flush; validate their separate history reference here.
                    requireNear(sample.previousKinematicsVelocity(), {2.0 * simulation.timeStep(), 0.0, 0.0}, "GPU observation must preserve the current prescribed velocity reference");
                else
                    requireNear(sample.acceleration(), {1.0, 0.0, 0.0}, "Prescribed acceleration must average the complete sampled velocity differences");
            }
        });
    simulation.step();
    simulation.observeCurrentState(
        [&](const solver&)
        {
            for (const virtualParticle& sample : simulation.virtualParticles().host())
            {
                if (mode == executionMode::GPU)
                    requireNear(sample.previousKinematicsVelocity(), {2.0 * simulation.timeStep(), 0.0, 0.0}, "GPU restore must preserve the prior prescribed velocity reference");
                else
                    requireNear(sample.acceleration(), {2.0 / 3.0, 0.0, 0.0}, "Observation restore must preserve the prior prescribed velocity reference");
            }
        });
}

void testPrescribedRotationDifference()
{
    LSParticleContainer owners;
    owners.host().push_back(syntheticOwner());
    owners.host()[0].setInfiniteMass();
    owners.host()[0].setAngularVelocity({0.0, 0.0, 1.0});
    virtualParticleContainer samples;
    samples.host().emplace_back(Vec3::unitX(), Vec3::unitX(), 1.0, 0);
    virtualParticleCoupling coupling;
    coupling.initialize(samples, owners, Vec3::zero());
    constexpr Real step = 0.01;
    const Vec3 oldVelocity{0.0, 1.0, 0.0};
    owners.host()[0].setOrientation({std::cos(0.5 * step), 0.0, 0.0, std::sin(0.5 * step)});
    const Vec3 currentVelocity = math::cross(owners.host()[0].angularVelocity(), math::rotateUnit(owners.host()[0].orientation(), Vec3::unitX()));
    coupling.accumulateKinematics(samples, owners, Vec3::zero(), step);
    coupling.averageKinematics(samples, owners, 1);
    requireNear(samples.host()[0].acceleration(), (currentVelocity - oldVelocity) / step, "Prescribed rotation must not double count centripetal acceleration");
}

void testCoupledMomentumAndObservations()
{
    forcedBoundarySolver uninterrupted;
    forcedBoundarySolver observed;
    uninterrupted.externalForce_ = observed.externalForce_ = {0.001, 0.0, 0.0};
    configureBoundary(uninterrupted, true);
    configureBoundary(observed, true);
    const Vec3 initialMomentum = momentum(uninterrupted);
    const Vec3 initialWallForce = uninterrupted.LSParticles().host()[0].force();
    bool loadChanged = false;
    for (int acoustic = 0; acoustic < 3; ++acoustic)
    {
        const int steps = uninterrupted.SPHStepInterval();
        for (int step = 0; step < steps; ++step)
        {
            uninterrupted.step();
            observed.step();
            if (step % 3 == 0)
            {
                const LSParticle before = observed.LSParticles().host()[0];
                observed.observeCurrentState([](const solver&) {});
                const LSParticle& after = observed.LSParticles().host()[0];
                requireNear(after.position(), before.position(), "An observation-only fluid projection must not move the rigid body", 0.0);
                requireNear(after.velocity(), before.velocity(), "An observation-only fluid projection must not change rigid velocity", 0.0);
                requireNear(after.angularVelocity(), before.angularVelocity(), "An observation-only fluid projection must not change rigid angular velocity", 0.0);
            }
        }
        uninterrupted.observeCurrentState([](const solver&) {});
        observed.observeCurrentState([](const solver&) {});
        requireNear(momentum(uninterrupted), initialMomentum + uninterrupted.time() * uninterrupted.externalForce_, "A real acoustic update must match coupled fluid/solid linear momentum");
        requireNear(momentum(observed), momentum(uninterrupted), "Observation cadence must preserve coupled momentum");
        requireNear(observed.SPHParticles().host()[0].position(), uninterrupted.SPHParticles().host()[0].position(), "Observation cadence must preserve coupled fluid trajectory");
        requireNear(observed.LSParticles().host()[0].velocity(), uninterrupted.LSParticles().host()[0].velocity(), "Observation cadence must preserve solid endpoint corrections");
        loadChanged = loadChanged || math::norm(uninterrupted.LSParticles().host()[0].force() - initialWallForce) > 1.0e-6;
    }
    require(loadChanged, "Coupled momentum test must include transient wall loads");
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
void testDevicePrescribedAcceleration()
{
    LSParticleContainer owners;
    owners.setUseDevice(true);
    owners.host().push_back(syntheticOwner());
    owners.host()[0].setInfiniteMass();
    owners.host()[0].setVelocity({1.0, 0.0, 0.0});
    virtualParticleContainer samples;
    samples.setUseDevice(true);
    samples.host().emplace_back(Vec3::unitX(), Vec3::unitX(), 1.0, 0);
    owners.copyHostToDeviceAsync();
    samples.copyHostToDeviceAsync();
    cuda::launchResetVirtualParticleKinematicsReference(samples, owners);
    const auto sample = [&]
    {
        cuda::launchClearVirtualParticleKinematics(samples);
        cuda::launchAccumulateVirtualParticleKinematics(samples, owners, Vec3::zero(), Real{0.01});
        cuda::launchAverageVirtualParticleKinematics(samples, owners, 1.0);
        samples.copyDeviceToHostAsync();
        host_device_detail::synchronize(nullptr);
    };
    sample();
    requireNear(samples.host()[0].acceleration(), Vec3::zero(), "GPU constant prescribed motion must have zero difference acceleration");
    owners.host()[0].setVelocity({1.02, 0.0, 0.0});
    owners.host()[0].setAngularVelocity({0.0, 0.0, 0.03});
    owners.copyHostFieldToDeviceAsync<pointMass::velocityField>();
    owners.copyHostFieldToDeviceAsync<rigidBody::angularVelocityField>();
    sample();
    requireNear(samples.host()[0].acceleration(), {2.0, 3.0, 0.0}, "GPU prescribed linear and angular acceleration must use velocity differences");
    const Vec3 previous = samples.host()[0].velocity();
    owners.host()[0].setOrientation({std::cos(0.005), 0.0, 0.0, std::sin(0.005)});
    owners.copyHostFieldToDeviceAsync<rigidBody::orientationField>();
    const Vec3 current = owners.host()[0].velocity() + math::cross(owners.host()[0].angularVelocity(), math::rotateUnit(owners.host()[0].orientation(), Vec3::unitX()));
    sample();
    requireNear(samples.host()[0].acceleration(), (current - previous) / 0.01, "GPU prescribed rotation must not double count centripetal acceleration");
}

void testDeviceHeldTorqueAndImpulse()
{
    LSParticleContainer owners;
    owners.setUseDevice(true);
    owners.host().push_back(syntheticOwner());
    virtualParticleContainer samples;
    samples.setUseDevice(true);
    samples.host().emplace_back(Vec3::unitX(), Vec3::unitX(), 1.0, 0);
    samples.host()[0].setForce({0.0, 1.0, 0.0});
    owners.copyHostToDeviceAsync();
    samples.copyHostToDeviceAsync();
    cuda::launchRefreshVirtualParticleCouplingLoads(samples, owners);
    const Real halfAngle = math::pi / 8.0;
    owners.host()[0].setOrientation({std::cos(halfAngle), 0.0, 0.0, std::sin(halfAngle)});
    owners.copyHostFieldToDeviceAsync<rigidBody::orientationField>();
    cuda::launchAddVirtualParticleForceAndTorque(owners, samples);
    owners.copyDeviceToHostAsync();
    host_device_detail::synchronize(nullptr);
    requireNear(owners.host()[0].torque(), {0.0, 0.0, 1.0}, "GPU must hold world torque while owner orientation changes");
    cuda::launchCacheVirtualParticleCouplingForce(samples);
    samples.host()[0].setForce({0.0, 4.0, 0.0});
    samples.copyHostFieldToDeviceAsync<pointMass::forceField>();
    cuda::launchRefreshVirtualParticleCouplingLoads(samples, owners);
    cuda::launchApplyVirtualParticleImpulseCorrection(owners, samples, 0.045);
    owners.copyDeviceToHostAsync();
    host_device_detail::synchronize(nullptr);
    const Vec3 newTorque = math::cross(math::rotateUnit(owners.host()[0].orientation(), Vec3::unitX()), Vec3{0.0, 4.0, 0.0});
    requireNear(2.0 * owners.host()[0].velocity(), 0.045 * Vec3{0.0, 3.0, 0.0}, "GPU force impulse correction must match the host convention");
    requireNear(3.0 * owners.host()[0].angularVelocity(), 0.045 * (newTorque - Vec3{0.0, 0.0, 1.0}), "GPU torque impulse correction must match the host convention");
}

void testRotatingBoundaryBackendParity(executionMode mode)
{
    SPHDEM reference;
    SPHDEM device{mode};
    configureBoundary(reference, true, true);
    configureBoundary(device, true, true);
    const Vec3 initialMomentum = momentum(reference);
    for (int acoustic = 0; acoustic < 2; ++acoustic)
    {
        const int steps = reference.SPHStepInterval();
        for (int step = 0; step < steps; ++step)
        {
            reference.step();
            device.step();
        }
        device.observeCurrentState([](const solver&) {});
        requireNear(momentum(device), initialMomentum, "Device acoustic boundary must conserve linear momentum", 1.0e-8);
        requireNear(device.LSParticles().host()[0].angularVelocity(),
                    reference.LSParticles().host()[0].angularVelocity(),
                    "Rotating finite-mass boundary angular velocity must agree across backends",
                    1.0e-8);
        requireNear(device.LSParticles().host()[0].velocity(), reference.LSParticles().host()[0].velocity(), "Finite-mass boundary velocity must agree across backends", 1.0e-8);
    }
}
#endif
} // namespace

int main()
{
    try
    {
        testHeldWorldTorqueAndImpulse();
        testCompleteEndpointAcceleration();
        testPrescribedAccelerationHistory(executionMode::CPU);
        testPrescribedRotationDifference();
        testCoupledMomentumAndObservations();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testDevicePrescribedAcceleration();
        testDeviceHeldTorqueAndImpulse();
        testRotatingBoundaryBackendParity(executionMode::GPU);
        testRotatingBoundaryBackendParity(executionMode::Hybrid);
        testPrescribedAccelerationHistory(executionMode::GPU);
        testPrescribedAccelerationHistory(executionMode::Hybrid);
#endif
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
