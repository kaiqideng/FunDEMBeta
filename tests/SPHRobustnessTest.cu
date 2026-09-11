#include "execution/cuda/sphInteractionKernel.cuh"
#include "execution/cuda/spatialGridKernel.cuh"
#include "solver/SPHDEM.h"
#include "interaction/SPHNeighborhood.h"
#include "interaction/sphInteraction.h"

#include <cmath>
#include <iostream>
#include <limits>
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

void requireNear(const Vec3& value, const Vec3& expected, const char* message)
{
    if (math::norm(value - expected) > 1.0e-9 * (1.0 + math::norm(expected)))
        throw std::runtime_error(message);
}

void testDeviceDisplacement()
{
    SPHParticleContainer particles;
    virtualParticleContainer boundaries;
    particles.setUseDevice(true);
    boundaries.setUseDevice(true);
    particles.host().emplace_back(Vec3::zero());
    particles.host().emplace_back(Vec3{0.2, 0.0, 0.0});
    boundaries.host().emplace_back(Vec3::zero(), Vec3::unitX(), 1.0, 0);
    particles.copyHostToDeviceAsync();
    boundaries.copyHostToDeviceAsync();
    SPHNeighborhood reference;
    reference.captureDevice(particles, boundaries, nullptr);
    const auto expired = [&] { return reference.exceedsSkin(cuda::maximumSPHNeighborDisplacementSquared(particles, boundaries, reference.positions(), nullptr), 0.02); };
    require(!expired(), "An unchanged device neighborhood must remain valid");
    particles.host()[1].setPosition({0.191, 0.0, 0.0});
    particles.copyHostFieldToDeviceAsync<pointMass::positionField>();
    require(!expired(), "Device displacement below half the skin must not rebuild");
    particles.host()[1].setPosition({0.025, 0.0, 0.0});
    particles.copyHostFieldToDeviceAsync<pointMass::positionField>();
    require(expired(), "Device cross-cell motion must expire the search skin");
    reference.captureDevice(particles, boundaries, nullptr);
    boundaries.host()[0].setPosition({0.011, 0.0, 0.0});
    boundaries.copyHostFieldToDeviceAsync<pointMass::positionField>();
    require(expired(), "Device moving-wall displacement must expire the search skin");
    reference.invalidate();
    require(expired(), "A restored output snapshot must invalidate device reference positions");
}

void testInvalidComputedDensityParity()
{
    SPHParticleContainer cpuParticles;
    SPHParticleContainer deviceParticles;
    virtualParticleContainer boundaries;
    deviceParticles.setUseDevice(true);
    boundaries.setUseDevice(true);
    cpuParticles.host().emplace_back(Vec3::zero());
    cpuParticles.host().emplace_back(Vec3{0.01, 0.0, 0.0});
    for (auto& particle : cpuParticles.host())
    {
        particle.setMass(1.0);
        particle.setDensity(1000.0);
        particle.setPressure(0.0);
        particle.setDensityRate(0.0);
    }
    cpuParticles.host()[1].*SPHParticle::pressureField::member = std::numeric_limits<Real>::quiet_NaN();
    deviceParticles.host() = cpuParticles.host();
    deviceParticles.copyHostToDeviceAsync();
    boundaries.copyHostToDeviceAsync();
    spatialGridContainer fluidGrid;
    spatialGridContainer wallGrid;
    fluidGrid.setUniform(deviceParticles, 0.06, {-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    wallGrid.setUniform(boundaries, 0.06, {-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    cuda::launchBuildSPHSpatialGrid(fluidGrid, deviceParticles);
    cuda::launchBuildVirtualParticleSpatialGrid(wallGrid, boundaries);
    cpu::SPHInteraction interaction;
    interaction.buildNeighborhood(cpuParticles, boundaries, {-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3}, 0.06, 0.08);
    const cpu::SPHInteractionParameters parameters{0.03, 1.0, 1000.0, 1.0, 10.0, 0.001, Vec3::zero()};
    const auto update = [&]
    {
        interaction.updateDensity(cpuParticles, boundaries, parameters, 1.0e-5);
        cuda::launchUpdateSPHDensity(deviceParticles, boundaries, fluidGrid.device(), wallGrid.device(), 0.03, 0.08, 1.0, 1000.0, 10.0, 1.0e-5, Vec3::zero());
        deviceParticles.copyDeviceToHostAsync();
        host_device_detail::synchronize(nullptr);
        require(!math::isFinite(cpuParticles.host()[0].densityRate()) && !math::isFinite(deviceParticles.host()[0].densityRate()), "Computed invalid density rates must survive on CPU and device");
        require(cpuParticles.host()[0].density() == deviceParticles.host()[0].density(), "CPU and GPU emergency density handling must agree");
        const auto cpuStatistics = interaction.stateStatistics(cpuParticles, Vec3::zero());
        const auto gpuStatistics = cuda::calculateSPHStateStatistics(deviceParticles, Vec3::zero());
        require(cpuStatistics.invalidValueCount_ > 0 && cpuStatistics.invalidValueCount_ == gpuStatistics.invalidValueCount_,
                "Computed-rate diagnostics must match even when density/pressure became finite");
    };
    update();
    require(cpuParticles.host()[0].density() == 100.0, "Invalid density integration must retain the shared emergency floor");
    // A later split stage may have finite input again; it must not erase failure.
    for (auto& particle : cpuParticles.host())
        particle.setPressure(0.0);
    for (auto& particle : deviceParticles.host())
        particle.setPressure(0.0);
    deviceParticles.copyHostFieldToDeviceAsync<SPHParticle::pressureField>();
    update();
    for (auto& particle : cpuParticles.host())
        particle.setConstrained(true);
    for (auto& particle : deviceParticles.host())
        particle.setConstrained(true);
    deviceParticles.copyHostFieldToDeviceAsync<SPHParticle::constrainedField>();
    update();
    interaction.initializeParticles(cpuParticles, 1000.0, 10.0);
    cuda::launchInitializeSPHParticles(deviceParticles, 1000.0, 10.0);
    require(interaction.stateStatistics(cpuParticles, Vec3::zero()).invalidValueCount_ == 0 && cuda::calculateSPHStateStatistics(deviceParticles, Vec3::zero()).invalidValueCount_ == 0,
            "Explicit initialization must reset the diagnostic latch on both backends");
}

class displacedSolver : public SPHDEM
{
public:
    explicit displacedSolver(executionMode mode) : SPHDEM(mode) {}
    void crossCellsDuringAdvection()
    {
        auto& particles = mutableSPHParticles();
        particles.host()[1].setPosition({0.025, 0.0, 0.0});
        for (auto& particle : particles.host())
        {
            particle.setDensity(1000.01);
            particle.setPressure(1.0);
        }
        particles.host()[0].setPriorForce({0.001, 0.0, 0.0});
        particles.host()[1].setPriorForce({-0.001, 0.0, 0.0});
        if (particles.usesDevice())
        {
            particles.copyHostFieldToDeviceAsync<pointMass::positionField>(gpuStream());
            particles.copyHostFieldToDeviceAsync<SPHParticle::densityField>(gpuStream());
            particles.copyHostFieldToDeviceAsync<SPHParticle::pressureField>(gpuStream());
            particles.copyHostFieldToDeviceAsync<SPHParticle::priorForceField>(gpuStream());
        }
    }
};

void configure(displacedSolver& simulation)
{
    simulation.setBoundary({-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    simulation.setGravity(Vec3::zero());
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(0.02, 0.03, 1000.0, 0.0);
    simulation.setSPHSoundSpeed(10.0);
    simulation.addSPHParticle(SPHParticle{Vec3::zero()});
    simulation.addSPHParticle(SPHParticle{{0.2, 0.0, 0.0}});
    simulation.initialize();
    require(simulation.SPHAdvectionStepInterval() >= 3, "Cross-cell test must reuse the same advection interval");
}

void testActualCrossCellGuard(executionMode mode)
{
    displacedSolver cpu{executionMode::CPU};
    displacedSolver device{mode};
    configure(cpu);
    configure(device);
    const int firstSteps = cpu.SPHStepInterval();
    for (int step = 0; step < firstSteps; ++step)
    {
        cpu.step();
        device.step();
    }
    // The initial advection neighborhood excludes the distant pair. Move it
    // into support without invalidating or starting another advection interval.
    cpu.crossCellsDuringAdvection();
    device.crossCellsDuringAdvection();
    const int nextSteps = cpu.SPHStepInterval();
    for (int step = 0; step < nextSteps; ++step)
    {
        cpu.step();
        device.step();
    }
    device.observeCurrentState([](const solver&) {});
    const auto& cpuParticles = cpu.SPHParticles().host();
    const auto& deviceParticles = device.SPHParticles().host();
    require(math::norm(cpuParticles[0].force() - cpuParticles[0].priorForce()) > 1.0e-6, "CPU actual-displacement guard must recover the approaching pressure pair");
    require(math::norm(deviceParticles[0].force() - deviceParticles[0].priorForce()) > 1.0e-6, "GPU actual-displacement guard must recover the approaching pressure pair");
    requireNear(deviceParticles[0].force(), cpuParticles[0].force(), "Pressure force after mid-advection grid rebuild must match CPU");
    requireNear(deviceParticles[0].velocity(), cpuParticles[0].velocity(), "Cross-cell CPU/device velocity must agree");
    requireNear(deviceParticles[0].velocity() + deviceParticles[1].velocity(), Vec3::zero(), "Recovered device pair must preserve momentum symmetry");
    requireNear(deviceParticles[0].priorForce(), {0.001, 0.0, 0.0}, "An early search rebuild must not refresh frozen prior forces");
    requireNear(cpuParticles[0].priorForce(), {0.001, 0.0, 0.0}, "CPU search rebuild must not reset the advection split");
}
} // namespace

int main()
{
    try
    {
        testDeviceDisplacement();
        testInvalidComputedDensityParity();
        testActualCrossCellGuard(executionMode::GPU);
        testActualCrossCellGuard(executionMode::Hybrid);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
