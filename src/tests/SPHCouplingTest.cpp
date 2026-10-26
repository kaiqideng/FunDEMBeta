#include "data/myLSObject.h"
#include "solver/SPHDEM.h"
#include "interaction/virtualParticleCoupling.h"

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "execution/cuda/integrationKernel.cuh"
#include "execution/cuda/sphInteractionKernel.cuh"
#endif

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

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
    if (!math::isFinite(value) || !math::isFinite(expected) || math::norm(value - expected) > tolerance * (1.0 + math::norm(expected)))
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

void testHeldWorldTorqueAndKDKBudget()
{
    LSParticleContainer owners;
    owners.host().push_back(syntheticOwner());
    virtualParticleContainer samples;
    samples.host().emplace_back(Vec3::unitX(), Vec3::unitX(), 1.0, 0);
    virtualParticleCoupling coupling;
    coupling.initialize(samples, owners, Vec3::zero());
    samples.host()[0].setForce({0.0, 1.0, 0.0});
    coupling.collectForceAndTorque(samples, owners);
    const auto oldState = coupling;

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
    // Refresh the acoustic reaction before the FIRST solid half kick. Every
    // KDK substep must use this same world wrench on both sides of its drift.
    samples.host()[0].setPosition(owner.position() + math::rotateUnit(owner.orientation(), samples.host()[0].localPosition()));
    samples.host()[0].setForce({0.0, 4.0, 0.0});
    coupling.collectForceAndTorque(samples, owners);
    const Vec3 newTorque = math::cross(samples.host()[0].position() - owner.position(), samples.host()[0].force());
    SPHParticleContainer fluid;
    fluid.host().emplace_back(Vec3::zero());
    fluid.host()[0].setMass(3.0);
    fluid.host()[0].setForce(-samples.host()[0].force());
    fluid.host()[0].updatePosition(0.5 * interval);
    fluid.host()[0].updateVelocity(Vec3::zero(), interval);
    fluid.host()[0].updatePosition(0.5 * interval);
    for (int index = 0; index < stepCount; ++index)
    {
        owner.setForce(externalForce);
        owner.setTorque(Vec3::zero());
        coupling.applyForceAndTorque(owners);
        owner.updateVelocityAndAngularVelocity(Vec3::zero(), 0.5 * step);
        owner.updatePositionAndOrientation(step);
        owner.setForce(externalForce);
        owner.setTorque(Vec3::zero());
        coupling.applyForceAndTorque(owners);
        owner.updateVelocityAndAngularVelocity(Vec3::zero(), 0.5 * step);
    }
    requireNear(owner.velocity() / owner.inverseMass(),
                interval * (samples.host()[0].force() + externalForce),
                "Held KDK loads must match the full fluid reaction without duplicating external impulse");
    requireNear(fluid.host()[0].mass() * fluid.host()[0].velocity() + owner.velocity() / owner.inverseMass(),
                interval * externalForce, "Fluid DKD and solid KDK kicks must exchange equal and opposite momentum");
    requireNear(3.0 * owner.angularVelocity(), interval * newTorque, "Every solid half kick must use the same held world torque");
    require(math::norm(newTorque - oldTorque) > 0.1, "Held-wrench test must change the prior reaction");

    coupling = oldState;
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
protected:
    void addLSParticleExternalForceAndTorque(LSParticleContainer::host_container_type& owners) override
    {
        for (LSParticle& owner : owners)
            owner.addForce(externalForce_);
    }
};

void configureBoundary(SPHDEM& simulation, bool nearbyFluid, bool rotating = false, bool infinite = false, Real dynamicViscosity = 0.0)
{
    constexpr Real spacing = 0.02;
    simulation.setBoundary({-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    simulation.setGravity(Vec3::zero());
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(spacing, 1.3 * spacing, 1000.0, dynamicViscosity);
    simulation.setSPHMaximumVelocity(1.0);
    simulation.setSPHSoundSpeed(10.0);
    simulation.addSPHParticle(SPHParticle{nearbyFluid ? Vec3{0.0, 0.0, 0.055} : Vec3{0.2, 0.2, 0.2}, nearbyFluid ? Vec3{0.0, 0.0, -0.1} : Vec3::zero()});
    const int materialIndex = simulation.addMaterial(material{1.0e9, 0.25, 1.0, 1000.0});
    levelset::BoxParticle box{{0.08, 0.08, 0.08}};
    box.buildLSGrid(spacing, 3);
    const int geometry = simulation.addGeometry(box);
    LSParticle owner;
    owner.setMaterial(simulation.materials(), materialIndex);
    owner.setGeometry(simulation.geometries(), geometry);
    if (rotating)
        owner.setAngularVelocity({0.0, 0.0, 2.0});
    if (infinite)
        owner.setInfiniteMass();
    simulation.addLSParticle(owner);
    simulation.initialize();
    require(infinite == (simulation.LSParticles().host()[0].inverseMass() == 0.0), "Coupling test boundary mass classification is incorrect");
    require(simulation.SPHStepInterval() > 2, "Coupling test must exercise multirate acoustic intervals");
}

Vec3 momentum(const SPHDEM& simulation)
{
    Vec3 result = Vec3::zero();
    for (const SPHParticle& particle : simulation.SPHParticles().host())
        result += particle.mass() * particle.velocity();
    for (const LSParticle& owner : simulation.LSParticles().host())
        if (owner.inverseMass() > 0.0)
            result += owner.velocity() / owner.inverseMass();
    return result;
}

struct synchronizedState
{
    SPHParticleContainer::host_container_type fluid;
    LSParticleContainer::host_container_type owners;
    Vec3 momentum;
    Real time;
    int stepCount;
};

synchronizedState currentState(SPHDEM& simulation)
{
    synchronizedState result;
    simulation.observeCurrentState([&](const solver& state)
    {
        result.fluid = simulation.SPHParticles().host();
        result.owners = simulation.LSParticles().host();
        result.momentum = momentum(simulation);
        result.time = state.time();
        result.stepCount = state.stepCount();
    });
    return result;
}

void requireSameState(const synchronizedState& first, const synchronizedState& second, const char* message, Real tolerance = 1.0e-10)
{
    require(first.fluid.size() == second.fluid.size() && first.owners.size() == second.owners.size(), message);
    require(first.stepCount == second.stepCount && std::abs(first.time - second.time) <= tolerance, message);
    for (std::size_t index = 0; index < first.fluid.size(); ++index)
    {
        const SPHParticle& a = first.fluid[index];
        const SPHParticle& b = second.fluid[index];
        requireNear(a.position(), b.position(), message, tolerance);
        requireNear(a.velocity(), b.velocity(), message, tolerance);
        requireNear({a.density(), a.pressure(), a.densityRate()}, {b.density(), b.pressure(), b.densityRate()}, message, tolerance);
    }
    for (std::size_t index = 0; index < first.owners.size(); ++index)
    {
        const LSParticle& a = first.owners[index];
        const LSParticle& b = second.owners[index];
        requireNear(a.position(), b.position(), message, tolerance);
        requireNear(a.velocity(), b.velocity(), message, tolerance);
        requireNear(a.angularVelocity(), b.angularVelocity(), message, tolerance);
        requireNear({a.orientation().x, a.orientation().y, a.orientation().z}, {b.orientation().x, b.orientation().y, b.orientation().z}, message, tolerance);
        require(std::abs(a.orientation().w - b.orientation().w) <= tolerance, message);
    }
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

void testPrescribedAccelerationHistory()
{
    LSParticleContainer owners;
    owners.host().push_back(syntheticOwner());
    owners.host()[0].setInfiniteMass();
    virtualParticleContainer samples;
    samples.host().emplace_back(Vec3::unitX(), Vec3::unitX(), 1.0, 0);
    virtualParticleCoupling coupling;
    coupling.initialize(samples, owners, Vec3::zero());
    constexpr Real step = 0.01;
    coupling.accumulateKinematics(samples, owners, Vec3::zero(), step);
    coupling.prepareKinematics(samples, owners, Vec3::zero(), 1);
    requireNear(samples.host()[0].acceleration(), Vec3::zero(), "Constant prescribed velocity must have zero sampled acceleration");
    owners.host()[0].setVelocity({2.0 * step, 0.0, 0.0});
    coupling.accumulateKinematics(samples, owners, Vec3::zero(), step);
    coupling.prepareKinematics(samples, owners, Vec3::zero(), 2);
    requireNear(samples.host()[0].acceleration(), {1.0, 0.0, 0.0}, "Prepared boundary acceleration must average the completed velocity differences");
    const auto saved = coupling;
    coupling.clearKinematics();
    coupling = saved;
    coupling.accumulateKinematics(samples, owners, Vec3::zero(), step);
    coupling.prepareKinematics(samples, owners, Vec3::zero(), 3);
    requireNear(samples.host()[0].acceleration(), {2.0 / 3.0, 0.0, 0.0}, "Restoring a boundary snapshot must preserve its finite-difference history");
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
    coupling.prepareKinematics(samples, owners, Vec3::zero(), 1);
    requireNear(samples.host()[0].acceleration(), (currentVelocity - oldVelocity) / step, "Prescribed rotation must not double count centripetal acceleration");
}

void testCoupledMomentumAndObservations(executionMode mode, Real dynamicViscosity = 0.0)
{
    forcedBoundarySolver uninterrupted{mode};
    forcedBoundarySolver observed{mode};
    // The CPU hook is deliberately a known non-fluid impulse budget. Device
    // cases use no user force so they exercise only the shared FSI exchange.
    if (mode != executionMode::GPU)
        uninterrupted.externalForce_ = observed.externalForce_ = {0.001, 0.0, 0.0};
    configureBoundary(uninterrupted, true, false, false, dynamicViscosity);
    configureBoundary(observed, true, false, false, dynamicViscosity);
    const auto initial = currentState(uninterrupted);
    const Vec3 initialWallForce = initial.owners[0].force();
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
                const int intervalBefore = observed.SPHStepInterval();
                const auto before = currentState(observed);
                observed.observeCurrentState([&](const solver& outer)
                {
                    observed.observeCurrentState([&](const solver& nested)
                    {
                        require(&outer == &nested && nested.stepCount() == before.stepCount, "Nested coupling observations must share one physical time");
                    });
                });
                bool propagated = false;
                try
                {
                    observed.observeCurrentState([](const solver&) { throw std::runtime_error("coupling observer failure"); });
                }
                catch (const std::runtime_error& error)
                {
                    propagated = std::string(error.what()) == "coupling observer failure";
                }
                require(propagated, "A throwing coupling observer must propagate its exception");
                requireSameState(before, currentState(observed), "Nested and throwing coupling observations must not alter current physical state");
                require(observed.SPHStepInterval() == intervalBefore, "Coupling observations must not re-quantize the acoustic interval");
            }
        }
        const auto reference = currentState(uninterrupted);
        const auto actual = currentState(observed);
        requireNear(reference.momentum, initial.momentum + uninterrupted.time() * uninterrupted.externalForce_, "A real acoustic interval must match the combined fluid/solid momentum budget", 1.0e-8);
        requireSameState(reference, actual, "Observation cadence must preserve coupled fluid and solid trajectories");
        loadChanged = loadChanged || math::norm(reference.owners[0].force() - initialWallForce) > 1.0e-6;
    }
    require(loadChanged, "Coupled momentum test must include transient wall loads");
}

void testCoupledSolveContinuation(executionMode mode)
{
    SPHDEM uninterrupted{mode};
    SPHDEM segmented{mode};
    SPHDEM singleStepSolves{mode};
    SPHDEM singleStepPauses{mode};
    configureBoundary(uninterrupted, true, true);
    configureBoundary(segmented, true, true);
    configureBoundary(singleStepSolves, true, true);
    configureBoundary(singleStepPauses, true, true);
    const int interval = uninterrupted.SPHStepInterval();
    const int firstPart = interval / 2 + 1;
    const int total = 2 * interval + 7;
    const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto directory = std::filesystem::temp_directory_path() / ("fundem-sph-coupling-" + std::to_string(stamp));
    segmented.setOutputDirectory(directory.string());
    segmented.setVTUOutputEnabled(false);
    segmented.setOutputStepInterval(3);
    uninterrupted.solve(total);
    segmented.solve(firstPart);
    const auto paused = currentState(segmented);
    require(paused.stepCount == firstPart && std::abs(paused.time - firstPart * segmented.timeStep()) < math::defaultTolerance,
            "A paused coupled solve must expose exactly its visible DEM time");
    segmented.writeOutput();
    requireSameState(paused, currentState(segmented), "Writing coupled output during a partial acoustic interval must not mutate continuation");
    segmented.solve(total - firstPart);
    requireSameState(currentState(uninterrupted), currentState(segmented), "Splitting a finite-mass coupled solve must preserve fluid, rigid-body and thermodynamic state");
    for (int step = 0; step < total; ++step)
    {
        singleStepSolves.solve(1);
        singleStepPauses.step();
        singleStepPauses.solve(0);
        requireSameState(currentState(singleStepSolves), currentState(singleStepPauses),
                         "Repeated one-step solves and manual one-step pauses must expose the same current-time state");
    }
    const auto expected = currentState(uninterrupted);
    requireSameState(expected, currentState(singleStepSolves), "Pausing and resuming at every DEM step must not restart the leading SPH acoustic interval");
    requireSameState(expected, currentState(singleStepPauses), "Manual step and solve(0) at every DEM step must preserve the uninterrupted coupled trajectory");
    uninterrupted.solve(9);
    segmented.solve(9);
    singleStepSolves.solve(9);
    singleStepPauses.solve(9);
    requireSameState(currentState(uninterrupted), currentState(segmented), "Coupled pause/output/resume must preserve later integration too");
    requireSameState(currentState(uninterrupted), currentState(singleStepSolves), "One-step solve pauses must retain the subsequent acoustic schedule");
    requireSameState(currentState(uninterrupted), currentState(singleStepPauses), "Manual one-step pauses must retain the subsequent acoustic schedule");
    std::filesystem::remove_all(directory);
}

void testAcceleratedCurrentTimeState(executionMode mode)
{
    SPHDEM simulation{mode};
    const Vec3 initialPosition{0.01, -0.02, 0.03};
    const Vec3 initialVelocity{0.4, 0.2, -0.1};
    const Vec3 acceleration{0.0, -2.0, 0.0};
    simulation.setBoundary({-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    simulation.setGravity(acceleration);
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(0.02, 0.026, 1000.0, 0.0);
    simulation.setSPHSoundSpeed(10.0);
    simulation.addSPHParticle(SPHParticle{initialPosition, initialVelocity});
    simulation.initialize();
    require(simulation.SPHStepInterval() > 2, "Accelerated current-time regression requires a multirate acoustic interval");
    const int totalSteps = 2 * simulation.SPHStepInterval() + 7;
    const auto verify = [&]
    {
        const auto& particle = simulation.SPHParticles().host()[0];
        const Real time = simulation.time();
        requireNear(particle.position(), initialPosition + time * initialVelocity + 0.5 * time * time * acceleration,
                    "Current-time SPH dense output must reproduce the exact constant-acceleration DKD trajectory", 1.0e-11);
        requireNear(particle.velocity(), initialVelocity + time * acceleration,
                    "Current-time SPH dense output must expose velocity at the DEM clock, not the leading acoustic endpoint", 1.0e-11);
    };
    for (int step = 1; step <= totalSteps; ++step)
    {
        simulation.step();
        simulation.observeCurrentState([&](const solver&) { verify(); });
        if (mode != executionMode::CPU)
        {
            simulation.copyDeviceToHost();
            verify();
        }
        if (step == totalSteps / 2)
        {
            simulation.solve(0);
            verify();
        }
    }
    simulation.solve(0);
    verify();
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
        cuda::launchPrepareVirtualParticleKinematics(samples, owners, Vec3::zero(), 1);
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

void testDeviceHeldTorqueAndKDKBudget()
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
    samples.host()[0].setForce({0.0, 4.0, 0.0});
    samples.copyHostFieldToDeviceAsync<pointMass::forceField>();
    cuda::launchRefreshVirtualParticleCouplingLoads(samples, owners);

    SPHParticleContainer fluid;
    fluid.setUseDevice(true);
    fluid.host().emplace_back(Vec3::zero());
    fluid.host()[0].setMass(3.0);
    fluid.host()[0].setForce({0.0, -4.0, 0.0});
    fluid.copyHostToDeviceAsync();
    constexpr int stepCount = 5;
    constexpr Real step = 0.01;
    constexpr Real interval = stepCount * step;
    cuda::launchIntegrateSPHPosition(fluid, 0.5 * interval);
    cuda::launchIntegrateSPHVelocity(fluid, Vec3::zero(), interval);
    cuda::launchIntegrateSPHPosition(fluid, 0.5 * interval);
    for (int index = 0; index < stepCount; ++index)
    {
        host_device_detail::checkCuda(cudaMemsetAsync(owners.device<pointMass::forceField>(), 0, sizeof(Vec3), nullptr), "Clear held-wrench test force");
        host_device_detail::checkCuda(cudaMemsetAsync(owners.device<rigidBody::torqueField>(), 0, sizeof(Vec3), nullptr), "Clear held-wrench test torque");
        cuda::launchAddVirtualParticleForceAndTorque(owners, samples);
        cuda::launchVelocityAndAngularVelocityIntegration(owners, Vec3::zero(), 0.5 * step);
        cuda::launchPositionAndOrientationIntegration(owners, step);
        host_device_detail::checkCuda(cudaMemsetAsync(owners.device<pointMass::forceField>(), 0, sizeof(Vec3), nullptr), "Clear held-wrench test force");
        host_device_detail::checkCuda(cudaMemsetAsync(owners.device<rigidBody::torqueField>(), 0, sizeof(Vec3), nullptr), "Clear held-wrench test torque");
        cuda::launchAddVirtualParticleForceAndTorque(owners, samples);
        cuda::launchVelocityAndAngularVelocityIntegration(owners, Vec3::zero(), 0.5 * step);
    }
    owners.copyDeviceToHostAsync();
    fluid.copyDeviceToHostAsync();
    host_device_detail::synchronize(nullptr);
    const Vec3 newTorque = math::cross(math::rotateUnit({std::cos(halfAngle), 0.0, 0.0, std::sin(halfAngle)}, Vec3::unitX()), Vec3{0.0, 4.0, 0.0});
    requireNear(2.0 * owners.host()[0].velocity(), interval * Vec3{0.0, 4.0, 0.0}, "GPU KDK must use the refreshed load from its first half kick");
    requireNear(3.0 * owners.host()[0].angularVelocity(), interval * newTorque, "GPU KDK must hold the evaluated world torque while the owner rotates");
    requireNear(fluid.host()[0].mass() * fluid.host()[0].velocity() + 2.0 * owners.host()[0].velocity(), Vec3::zero(),
                "GPU fluid DKD and solid KDK must exchange equal and opposite momentum");
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
            if (step == 1 || step + 1 == steps)
            {
                virtualParticleContainer::host_container_type referenceSamples;
                reference.observeCurrentState([&](const solver&) { referenceSamples = reference.virtualParticles().host(); });
                device.observeCurrentState([&](const solver&)
                {
                    const auto& deviceSamples = device.virtualParticles().host();
                    require(deviceSamples.size() == referenceSamples.size(), "Current-time boundary observations must retain identical sample ownership across backends");
                    for (std::size_t index = 0; index < referenceSamples.size(); ++index)
                    {
                        requireNear(deviceSamples[index].position(), referenceSamples[index].position(), "Current-time virtual boundary positions must match CPU at interior and completed acoustic steps", 1.0e-8);
                        requireNear(deviceSamples[index].velocity(), referenceSamples[index].velocity(), "Current-time virtual boundary velocity must match CPU at interior and completed acoustic steps", 1.0e-8);
                        requireNear(deviceSamples[index].acceleration(), referenceSamples[index].acceleration(), "Current-time virtual boundary acceleration must match CPU at interior and completed acoustic steps", 1.0e-8);
                    }
                });
            }
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
        testHeldWorldTorqueAndKDKBudget();
        testCompleteEndpointAcceleration();
        testPrescribedAccelerationHistory();
        testPrescribedRotationDifference();
        testCoupledMomentumAndObservations(executionMode::CPU);
        testCoupledMomentumAndObservations(executionMode::CPU, 0.2);
        testCoupledSolveContinuation(executionMode::CPU);
        testAcceleratedCurrentTimeState(executionMode::CPU);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testDevicePrescribedAcceleration();
        testDeviceHeldTorqueAndKDKBudget();
        testCoupledMomentumAndObservations(executionMode::GPU);
        testCoupledMomentumAndObservations(executionMode::Hybrid);
        testCoupledMomentumAndObservations(executionMode::GPU, 0.2);
        testCoupledMomentumAndObservations(executionMode::Hybrid, 0.2);
        testCoupledSolveContinuation(executionMode::GPU);
        testCoupledSolveContinuation(executionMode::Hybrid);
        testAcceleratedCurrentTimeState(executionMode::GPU);
        testAcceleratedCurrentTimeState(executionMode::Hybrid);
        testRotatingBoundaryBackendParity(executionMode::GPU);
        testRotatingBoundaryBackendParity(executionMode::Hybrid);
#endif
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
