#include "data/myLSObject.h"
#include "execution/sphJetFunctions.h"
#include "execution/sphSourceDiscretization.h"
#include "solver/SPHDEM.h"
#include "particle/SPHJet.h"

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "execution/cuda/sphJetKernel.cuh"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using namespace fundem;

class testSPHDEM : public SPHDEM
{
public:
    using SPHDEM::SPHDEM;
    SPHParticleContainer& particles() noexcept { return mutableSPHParticles(); }
    void repositionParticle(int index, const math::Vec3& position)
    {
        mutableSPHParticles().host()[index].setPosition(position);
        invalidateContinuation();
    }
};

bool nearlyEqual(math::Real first, math::Real second, math::Real tolerance = 1.0e-12) noexcept { return std::abs(first - second) <= tolerance * (1.0 + std::max(std::abs(first), std::abs(second))); }

void require(bool condition, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void testValueType()
{
    const SPHJet jet{{1.0, 2.0, 3.0}, {3.0, 4.0, 0.0}, 0.25, 2.0, 1.5};
    require(jet.isValid(), "A finite positive SPH jet must be valid.");
    require(math::nearlyEqual(jet.direction(), {0.6, 0.8, 0.0}), "SPHJet must store a unit direction.");
    require(math::nearlyEqual(jet.velocity(), {1.2, 1.6, 0.0}), "SPHJet velocity must equal speed times unit direction.");
    require(nearlyEqual(jet.pipeLength(), 3.0), "The nominal finite jet stroke must equal speed times duration, not the bounded inlet-buffer length.");

    require(!SPHJet{}.isValid(), "A default SPH jet must be invalid.");
    require(!SPHJet{{0.0, 0.0, 0.0}, math::Vec3::zero(), 0.25, 2.0, 1.5}.isValid(), "A zero-direction SPH jet must be rejected.");
    require(!SPHJet{{0.0, 0.0, 0.0}, math::Vec3::unitX(), 0.0, 2.0, 1.5}.isValid(), "A zero-radius SPH jet must be rejected.");
    require(!SPHJet{{0.0, 0.0, 0.0}, math::Vec3::unitX(), 0.25, 0.0, 1.5}.isValid(), "A zero-speed SPH jet must be rejected.");
    require(!SPHJet{{0.0, 0.0, 0.0}, math::Vec3::unitX(), 0.25, 2.0, 0.0}.isValid(), "A zero-duration SPH jet must be rejected.");
}

void testPipeGeometry()
{
    const math::Vec3 outlet{1.0, 2.0, 3.0};
    require(execution::isInsideSPHJetPipe({0.0, 2.0, 3.0}, outlet, math::Vec3::unitX(), 1.0, 2.0), "An upstream axial point must lie inside the jet pipe.");
    require(execution::isInsideSPHJetPipe({0.0, 3.0, 3.0}, outlet, math::Vec3::unitX(), 1.0, 2.0), "The cylindrical radial boundary must be included.");
    require(!execution::isInsideSPHJetPipe(outlet, outlet, math::Vec3::unitX(), 1.0, 2.0), "The outlet plane must not retain emitted particles.");
    require(!execution::isInsideSPHJetPipe({1.1, 2.0, 3.0}, outlet, math::Vec3::unitX(), 1.0, 2.0), "A downstream point must not lie inside the jet pipe.");
    require(!execution::isInsideSPHJetPipe({-1.1, 2.0, 3.0}, outlet, math::Vec3::unitX(), 1.0, 2.0), "A point beyond the upstream end must not lie inside the jet pipe.");
    require(!execution::isInsideSPHJetPipe({0.0, 3.01, 3.0}, outlet, math::Vec3::unitX(), 1.0, 2.0), "A point outside the jet radius must not lie inside the jet pipe.");
}

void testDiscretization()
{
    execution::SPHJetDiscretization discretization;
    require(execution::trySPHJetDiscretization(0.125, 0.1875, 0.5, discretization), "A representable SPH jet discretization must be accepted.");
    require(discretization.radialIndex_ == 1 && discretization.axialCount_ == 4, "The jet lattice bounds must match radius and length.");
    require(discretization.particleCount_ == 20, "The jet particle count must be correct.");
    require(discretization.contains(1, 0) && !discretization.contains(1, 1), "The transverse square lattice must be clipped to a circle.");

    require(!execution::trySPHJetDiscretization(0.125, 0.06, 0.5, discretization), "A jet narrower than one centered particle must be rejected.");
    require(!execution::trySPHJetDiscretization(0.125, 0.1875, 0.1, discretization), "A jet shorter than one axial layer must be rejected.");
}

void testOpeningGeometry()
{
    const execution::SPHJetOpening opening{{math::Vec3::zero(), math::Vec3::unitX(), 0.02, 0.08}};
    constexpr math::Real clearance = 0.01;
    require(execution::isInsideSPHJetOpening({-0.08, 0.0, 0.0}, opening.value_, clearance), "The finite upstream inlet buffer must be open.");
    require(execution::isInsideSPHJetOpening(math::Vec3::zero(), opening.value_, clearance), "The outlet plane must be open to virtual boundary samples.");
    require(execution::isInsideSPHJetOpening({0.005, 0.0, 0.0}, opening.value_, clearance), "The outlet cut must include half-spacing clearance.");
    require(!execution::isInsideSPHJetOpening({0.011, 0.0, 0.0}, opening.value_, clearance), "The opening must not cut a downstream target body.");
    require(!execution::isInsideSPHJetOpening({-0.091, 0.0, 0.0}, opening.value_, clearance), "The opening must not extend behind the short upstream buffer.");
    require(!execution::isInsideSPHJetOpening({0.0, 0.031, 0.0}, opening.value_, clearance), "Wall samples next to the inlet rim must remain active.");
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
void testDeviceOpeningMask()
{
    struct scopedStream {
        cudaStream_t value_{nullptr};
        scopedStream() { host_device_detail::checkCuda(cudaStreamCreateWithFlags(&value_, cudaStreamNonBlocking), "Create jet-opening test stream"); }
        ~scopedStream() { cudaStreamDestroy(value_); }
    } stream;
    constexpr math::Real nominalVolume = 0.001;
    constexpr math::Real clearance = 0.01;
    const math::Vec3 force{1.0, 2.0, 3.0};
    const math::Vec3 priorForce{4.0, 5.0, 6.0};
    const math::Vec3 torque{7.0, 8.0, 9.0};
    virtualParticleContainer samples;
    samples.setUseDevice(true);
    for (const math::Vec3 position : {math::Vec3{-0.04, 0.0, 0.0}, math::Vec3{-0.04, 0.025, 0.0}, math::Vec3{-0.04, 0.05, 0.0},
                                     math::Vec3{0.03, 0.0, 0.0}, math::Vec3{-0.04, 0.15, 0.0}})
    {
        samples.host().emplace_back(position, math::Vec3::unitX(), nominalVolume, 0);
        virtualParticle& sample = samples.host().back();
        sample.setPosition(position);
        sample.setForce(force);
        sample.setPriorForce(priorForce);
        sample.setCouplingTorque(torque);
    }
    execution::SPHJetOpeningContainer openings;
    openings.setUseDevice(true);
    for (math::Real y : {0.0, 0.03})
        openings.host().push_back(execution::SPHJetOpening{{{0.0, y, 0.0}, math::Vec3::unitX(), 0.02, 0.08}});
    samples.copyHostToDeviceAsync(stream.value_);
    openings.copyHostToDeviceAsync(stream.value_);
    require(cuda::launchUpdateSPHJetOpenings(samples, openings, nominalVolume, clearance, stream.value_), "The device must report its initial opening-mask change.");
    samples.copyDeviceToHost(stream.value_);
    const auto requireSample = [&](int index, bool masked, bool cleared)
    {
        const virtualParticle& sample = samples.host()[index];
        require(sample.volume() == (masked ? 0.0 : nominalVolume), "The device opening kernel produced an incorrect effective volume.");
        require(sample.force() == (cleared ? math::Vec3::zero() : force) && sample.priorForce() == (cleared ? math::Vec3::zero() : priorForce) &&
                    sample.couplingTorque() == (cleared ? math::Vec3::zero() : torque),
                "An opening-mask change must clear held load caches without changing untouched samples.");
    };
    requireSample(0, true, true);
    requireSample(1, true, true);
    requireSample(2, true, true);
    requireSample(3, false, false);
    requireSample(4, false, false);
    require(!cuda::launchUpdateSPHJetOpenings(samples, openings, nominalVolume, clearance, stream.value_), "An unchanged device opening mask must not trigger a refresh.");

    // Only positions change: the kernel itself must restore/cut volume and clear stale loads.
    samples.host()[0].setPosition({-0.04, 0.2, 0.0});
    samples.host()[3].setPosition({-0.04, 0.015, 0.0});
    samples.copyHostFieldToDeviceAsync<pointMass::positionField>(stream.value_);
    require(cuda::launchUpdateSPHJetOpenings(samples, openings, nominalVolume, clearance, stream.value_), "Moving samples across an inlet rim must update the actual device mask.");
    samples.copyDeviceToHost(stream.value_);
    requireSample(0, false, true);
    requireSample(1, true, true);
    requireSample(2, true, true);
    requireSample(3, true, true);
    requireSample(4, false, false);
    require(!cuda::launchUpdateSPHJetOpenings(samples, openings, nominalVolume, clearance, stream.value_), "A repeated moving-wall mask must remain unchanged.");
}
#endif

constexpr math::Real openingTestSpacing = 0.02;
constexpr math::Real openingTestDuration = 0.02;

void configureOpening(testSPHDEM& simulation, bool moving)
{
    simulation.setBoundary({-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    simulation.setGravity(math::Vec3::zero());
    simulation.setTimeStep(1.0e-4);
    simulation.setSPHProperties(openingTestSpacing, 1.3 * openingTestSpacing, 1000.0, 1.0e-3);
    simulation.setSPHMaximumVelocity(3.0);
    simulation.setVTUOutputEnabled(false);
    for (math::Real y : {-0.04, 0.04})
        simulation.addSPHJet(SPHJet{{0.0, y, 0.0}, math::Vec3::unitX(), 0.5 * openingTestSpacing, 1.0, openingTestDuration});

    const int materialIndex = simulation.addMaterial(material{1.0e9, 0.25, 1.0, 1000.0});
    const auto addWall = [&](const math::Vec3& size, const math::Vec3& position, const math::Vec3& velocity)
    {
        levelset::BoxParticle box{size};
        box.buildLSGrid(openingTestSpacing, 3);
        const int geometryIndex = simulation.addGeometry(box);
        LSParticle wall;
        wall.setMaterial(simulation.materials(), materialIndex);
        wall.setGeometry(simulation.geometries(), geometryIndex);
        wall.setInfiniteMass();
        wall.setPosition(position);
        wall.setVelocity(velocity);
        simulation.addLSParticle(wall);
    };
    addWall({0.08, 0.2, 0.08}, {-0.04, 0.0, 0.0}, moving ? math::Vec3{0.0, 2.0, 0.0} : math::Vec3::zero());
    addWall({0.04, 0.04, 0.04}, {0.18, 0.0, 0.0}, math::Vec3::zero());
    simulation.initialize();
}

std::vector<math::Real> openingVolumes(SPHDEM& simulation)
{
    std::vector<math::Real> volumes;
    simulation.observeCurrentState([&](const solver&)
    {
        int masked = 0;
        int neighbouringWall = 0;
        int downstreamWall = 0;
        const math::Real cellVolume = openingTestSpacing * openingTestSpacing * openingTestSpacing;
        for (const virtualParticle& sample : simulation.virtualParticles().host())
        {
            bool insideOpening = false;
            for (math::Real y : {-0.04, 0.04})
            {
                const execution::SPHJetOpening opening{{{0.0, y, 0.0}, math::Vec3::unitX(), 0.5 * openingTestSpacing, openingTestSpacing}};
                insideOpening = insideOpening || execution::isInsideSPHJetOpening(sample.position(), opening.value_, 0.5 * openingTestSpacing);
            }
            require(nearlyEqual(sample.volume(), insideOpening ? 0.0 : cellVolume), "Virtual inlet openings must match the finite union of source buffers.");
            if (sample.ownerLSParticleIndex() == 0)
            {
                masked += insideOpening ? 1 : 0;
                neighbouringWall += insideOpening ? 0 : 1;
            }
            else
            {
                require(!insideOpening && sample.volume() > 0.0, "A downstream target wall must not be masked by the upstream nozzle opening.");
                ++downstreamWall;
            }
            volumes.push_back(sample.volume());
        }
        require(masked > 0 && neighbouringWall > 0 && downstreamWall > 0, "The opening fixture must exercise cut, neighbouring, and downstream wall samples.");
    });
    return volumes;
}

void requireSameSynchronizedState(SPHDEM& first, SPHDEM& second);

void testVirtualOpeningAndContinuation(executionMode mode, bool moving)
{
    testSPHDEM reference{mode};
    testSPHDEM observed{mode};
    testSPHDEM segmented{mode};
    for (testSPHDEM* simulation : {&reference, &observed, &segmented})
        configureOpening(*simulation, moving);
    const auto initialVolumes = openingVolumes(reference);
    const auto originalGrid = reference.geometries().hostGridNodes();
    const auto originalSamples = reference.virtualParticles().host();
    constexpr int finalStep = 240;
    reference.solve(finalStep);
    for (int step = 1; step <= finalStep; ++step)
    {
        observed.step();
        if (step % 19 == 0)
            openingVolumes(observed);
    }
    observed.solve(0);
    segmented.solve(73);
    openingVolumes(segmented);
    segmented.solve(88);
    openingVolumes(segmented);
    segmented.solve(finalStep - 161);

    requireSameSynchronizedState(reference, observed);
    requireSameSynchronizedState(reference, segmented);
    const auto finalVolumes = openingVolumes(reference);
    require(finalVolumes == openingVolumes(observed) && finalVolumes == openingVolumes(segmented), "Output or pause/resume changed the virtual opening mask.");
    require(reference.SPHJetsCompleted(), "Opening persistence must be checked after finite source completion.");
    require(finalVolumes.size() == initialVolumes.size(), "Opening masks must preserve virtual-particle slots and owner mapping.");
    if (moving)
    {
        bool restored = false;
        bool newlyMasked = false;
        for (std::size_t index = 0; index < finalVolumes.size(); ++index)
        {
            restored = restored || (initialVolumes[index] == 0.0 && finalVolumes[index] > 0.0);
            newlyMasked = newlyMasked || (initialVolumes[index] > 0.0 && finalVolumes[index] == 0.0);
        }
        require(restored && newlyMasked, "Moving walls must both restore samples leaving an opening and mask samples entering it.");
    }
    else
        require(finalVolumes == initialVolumes, "A static inlet opening must remain open after emission ends.");
    for (std::size_t index = 0; index < originalSamples.size(); ++index)
    {
        const auto& sample = reference.virtualParticles().host()[index];
        require(sample.ownerLSParticleIndex() == originalSamples[index].ownerLSParticleIndex() && sample.localPosition() == originalSamples[index].localPosition(),
                "Opening updates must not reorder samples or alter owner-local geometry.");
    }
    const auto& currentGrid = reference.geometries().hostGridNodes();
    require(currentGrid.size() == originalGrid.size(), "Opening masks must not resize level-set geometry.");
    for (std::size_t index = 0; index < originalGrid.size(); ++index)
        require(currentGrid[index].signedDistance_ == originalGrid[index].signedDistance_, "Opening masks must not modify signed-distance grid values.");
}

void testGenerationConstraintAndCompletion(executionMode mode)
{
    constexpr math::Real spacing = 1.0e-3;
    constexpr math::Real duration = 1.0e-3;
    testSPHDEM simulation{mode};
    simulation.setBoundary({-0.01, -0.01, -0.01}, {0.01, 0.01, 0.01});
    simulation.setGravity({0.0, 2.0, 0.0});
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(spacing, 1.3 * spacing, 1000.0, 1.0e-3);
    simulation.setSPHMaximumVelocity(1.0);

    const SPHJet jet{{0.0, 0.0, 0.0}, math::Vec3::unitX(), 0.5 * spacing, 1.0, duration};
    require(simulation.addSPHJet(jet) == 0, "The first jet must own the first generated SPH particle.");
    require(simulation.SPHParticles().hostSize() == 1, "A one-layer, one-column jet must generate exactly one particle.");
    const SPHParticle& generated = simulation.SPHParticles().host()[0];
    require(math::nearlyEqual(generated.position(), {-0.5 * spacing, 0.0, 0.0}), "Jet particles must be generated upstream from the outlet.");
    require(math::nearlyEqual(generated.velocity(), jet.velocity()), "Generated jet particles must start at the prescribed velocity.");
    simulation.particles().host()[0].setDensity(1100.0);

    simulation.solve(0);
    require(!simulation.SPHJetsCompleted(), "A finite-duration jet must remain active at its start time.");
    require(simulation.SPHParticles().host()[0].isConstrained(), "A jet particle inside the virtual pipe must carry the constraint marker.");
    simulation.solve(1);
    require(math::nearlyEqual(simulation.SPHParticles().host()[0].velocity(), jet.velocity()), "The virtual pipe must restore prescribed velocity through a complete acoustic substep.");
    require(nearlyEqual(simulation.SPHParticles().host()[0].density(), 1000.0), "The inlet must impose reference density rather than retain a user-modified reservoir density.");
    simulation.solve(100);
    require(simulation.SPHJetsCompleted(), "The jet completion marker must be set after its duration has elapsed.");
    require(!simulation.SPHParticles().host()[0].isConstrained(), "The particle constraint must clear when the finite-duration jet ends.");
    simulation.solve(1);
    require(simulation.SPHJetsCompleted(), "A completed jet must remain permanently marked complete across later solve calls.");
}

void requireSameState(const SPHDEM& first, const SPHDEM& second)
{
    require(first.SPHJetsCompleted() == second.SPHJetsCompleted(), "Jet completion must survive output and solve continuation.");
    require(first.SPHParticles().hostSize() == second.SPHParticles().hostSize(), "Jet snapshots must preserve the particle count.");
    for (std::size_t index = 0; index < first.SPHParticles().hostSize(); ++index)
    {
        const SPHParticle& a = first.SPHParticles().host()[index];
        const SPHParticle& b = second.SPHParticles().host()[index];
        require(math::nearlyEqual(a.position(), b.position()) && math::nearlyEqual(a.velocity(), b.velocity()) &&
                    nearlyEqual(a.density(), b.density()) && a.isConstrained() == b.isConstrained(),
                "Jet output or split solve changed the particle continuation state.");
    }
}

void requireSameSynchronizedState(SPHDEM& first, SPHDEM& second)
{
    first.observeCurrentState([&](const solver&)
    {
        second.observeCurrentState([&](const solver&) { requireSameState(first, second); });
    });
}

void testMultipleJetSnapshotContinuation(executionMode mode)
{
    constexpr math::Real spacing = 1.0e-3;
    constexpr math::Real timeStep = 1.0e-5;
    constexpr math::Real lastEndTime = 156.5 * timeStep;
    testSPHDEM reference{mode};
    testSPHDEM observed{mode};
    testSPHDEM split{mode};
    for (testSPHDEM* simulation : {&reference, &observed, &split})
    {
        require(simulation->SPHJetsCompleted(), "A solver without jets must report completion.");
        simulation->setBoundary({-0.01, -0.01, -0.01}, {0.01, 0.01, 0.01});
        simulation->setGravity({0.0, 2.0, 0.0});
        simulation->setTimeStep(timeStep);
        simulation->setSPHProperties(spacing, 1.3 * spacing, 1000.0, 1.0e-3);
        simulation->setSPHMaximumVelocity(1.0);
        // Add the longer jet first: the completion bound must be the maximum,
        // independent of insertion order.
        require(simulation->addSPHJet(SPHJet{{0.0, 0.004, 0.0}, math::Vec3::unitX(), 0.5 * spacing, 1.0, lastEndTime}) == 0,
                "The longer jet must own the first particle.");
        require(simulation->addSPHJet(SPHJet{{0.0, -0.004, 0.0}, math::Vec3::unitX(), 0.5 * spacing, 1.0, 103.5 * timeStep}) == 1,
                "The shorter jet must own a separate particle.");
    }

    const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const std::filesystem::path outputDirectory = std::filesystem::temp_directory_path() / ("fundem-sph-jet-state-" + std::to_string(stamp));
    observed.setOutputDirectory((outputDirectory / "observed").string());
    split.setOutputDirectory((outputDirectory / "split").string());
    observed.solve(0);
    for (int step = 1; step <= 200; ++step)
    {
        reference.step();
        observed.step();
        if (step == 104 || step == 156 || step == 157)
        {
            require(observed.SPHJetsCompleted() == (step == 157), "Only the latest jet end may complete all jets at the visible time.");
            observed.observeCurrentState([&](const solver&)
            {
                require(observed.SPHJetsCompleted() == (step == 157), "A current-state snapshot must report the visible jet completion time.");
                require(!observed.SPHParticles().host()[0].isConstrained(), "A seed that has crossed the outlet must remain free before and after source completion.");
            });
            requireSameSynchronizedState(reference, observed);
            observed.writeOutput();
        }
        requireSameSynchronizedState(reference, observed);
    }
    reference.solve(0);
    observed.solve(0);

    split.solve(104);
    require(!split.SPHJetsCompleted(), "The shorter jet ending must not complete the longer jet.");
    split.writeOutput();
    split.solve(52);
    require(!split.SPHJetsCompleted(), "The final jet must remain active just before its end.");
    split.solve(1);
    require(split.SPHJetsCompleted() && !split.SPHParticles().host()[0].isConstrained(), "Solve completion must expose the cleared constraint after the final jet ends.");
    split.writeOutput();
    split.solve(43);
    requireSameSynchronizedState(reference, observed);
    requireSameSynchronizedState(reference, split);

    reference.solve(19);
    observed.solve(19);
    split.solve(19);
    requireSameSynchronizedState(reference, observed);
    requireSameSynchronizedState(reference, split);
    require(split.SPHJetsCompleted(), "All jets must remain complete after further solve continuation.");
    for (const SPHParticle& particle : split.SPHParticles().host())
        require(!particle.isConstrained(), "No completed jet may retain a particle constraint after continuation.");
    std::filesystem::remove_all(outputDirectory);
}

void configureGrowingJets(testSPHDEM& simulation)
{
    constexpr math::Real spacing = 1.0e-3;
    simulation.setBoundary({-0.05, -0.05, -0.05}, {0.05, 0.05, 0.05});
    simulation.setGravity(math::Vec3::zero());
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(spacing, 1.3 * spacing, 1000.0, 1.0e-3);
    simulation.setSPHMaximumVelocity(1.0);
    require(simulation.addSPHBlock({0.02, 0.02, 0.02}, {1, 1, 1}, {0.125, 0.0, 0.0}) == 0,
            "The ordinary prefix must retain its insertion index.");
    require(simulation.addSPHJet(SPHJet{{0.0, 0.006, 0.0}, math::Vec3::unitX(), 0.5 * spacing, 1.0, 0.0125}) == 1,
            "A long jet must append its short initial buffer after the block.");
    require(simulation.SPHJetParticleIndices(1).size() == 4, "The first jet must seed ceil(2h/spacing)+1 layers, not its complete duration.");
    require(simulation.addSPHJet(SPHJet{{0.0, -0.006, 0.0}, math::Vec3::unitX(), 0.5 * spacing, 0.5, 0.0165}) == 5,
            "The second jet must be placed immediately after the first short buffer.");
    require(simulation.SPHJetParticleIndices(5).size() == 4, "Each jet must own only its own short seed buffer.");
    require(simulation.addSPHBlock({-0.02, 0.02, 0.02}, {1, 1, 1}, {-0.125, 0.0, 0.0}) == 9,
            "Ordinary particles appended after source setup must keep stable IDs.");
    require(simulation.SPHParticles().hostSize() == 10, "Future source layers must not exist in the initial particle container.");
}

void requireJetOwnership(const SPHDEM& simulation, int firstIndex, int expectedCount, math::Real speed, math::Real outletY)
{
    constexpr math::Real spacing = 1.0e-3;
    const auto& indices = simulation.SPHJetParticleIndices(firstIndex);
    require(static_cast<int>(indices.size()) == expectedCount, "A finite jet emitted an incorrect total number of particles.");
    for (int layer = 0; layer < expectedCount; ++layer)
    {
        const int index = indices[layer];
        require(index >= 0 && index < static_cast<int>(simulation.SPHParticles().hostSize()), "A source contains an invalid stable particle index.");
        const SPHParticle& particle = simulation.SPHParticles().host()[index];
        require(math::nearlyEqual(particle.position(), {speed * simulation.time() - (layer + 0.5) * spacing, outletY, 0.0}),
                "Source birth phase or device capacity growth changed the regular emitted-particle trajectory.");
        require(math::nearlyEqual(particle.velocity(), {speed, 0.0, 0.0}) && nearlyEqual(particle.density(), 1000.0),
                "Source growth changed an existing particle or did not initialize an emitted particle.");
        require(!particle.isConstrained(), "Particles must be free after their finite source has completed.");
    }
}

void testGrowingJetOwnershipAndContinuation(executionMode mode)
{
    constexpr int finalStep = 1701;
    testSPHDEM uninterrupted{mode};
    testSPHDEM observed{mode};
    testSPHDEM segmented{mode};
    for (testSPHDEM* simulation : {&uninterrupted, &observed, &segmented})
        configureGrowingJets(*simulation);

    const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const std::filesystem::path outputDirectory = std::filesystem::temp_directory_path() / ("fundem-growing-jets-" + std::to_string(stamp));
    observed.setOutputDirectory((outputDirectory / "observed").string());
    observed.setVTUOutputEnabled(false);
    segmented.setOutputDirectory((outputDirectory / "segmented").string());
    segmented.setVTUOutputEnabled(false);

    uninterrupted.solve(finalStep);
    for (int step = 1; step <= finalStep; ++step)
    {
        observed.step();
        if (step == 99)
        {
            observed.observeCurrentState([&](const solver&)
            {
                require(observed.SPHParticles().hostSize() == 10 && observed.SPHJetParticleIndices(1).size() == 4,
                        "Dense output exposed a source layer before its birth time.");
            });
        }
        if (step % 17 == 0)
        {
            const std::size_t emitted = observed.SPHJetParticleIndices(1).size();
            observed.observeCurrentState([](const solver&) {});
            require(observed.SPHJetParticleIndices(1).size() == emitted, "An output observation emitted particles or advanced a source counter.");
        }
        if (step == 99 || step == 101 || step == 201 || step == 997)
            observed.writeOutput();
    }
    observed.solve(0);

    segmented.solve(99);
    require(segmented.SPHParticles().hostSize() == 10 && segmented.SPHJetParticleIndices(1).size() == 4,
            "Pausing before the first scheduled birth must not generate a future layer.");
    segmented.solve(2);
    require(segmented.SPHJetParticleIndices(1).size() == 5 && segmented.SPHJetParticleIndices(5).size() == 4,
            "The first scheduled layer must be appended once, to its own source.");
    segmented.writeOutput();
    segmented.solve(98);
    require(segmented.SPHJetParticleIndices(1).size() == 5 && segmented.SPHJetParticleIndices(5).size() == 4,
            "The next two source layers must not be emitted early.");
    segmented.solve(2);
    require(segmented.SPHJetParticleIndices(1).size() == 6 && segmented.SPHJetParticleIndices(5).size() == 5,
            "Simultaneous source births must preserve both source ownership lists.");
    segmented.writeOutput();
    segmented.solve(finalStep - 201);

    requireSameSynchronizedState(uninterrupted, observed);
    requireSameSynchronizedState(uninterrupted, segmented);
    for (const testSPHDEM* simulation : {&uninterrupted, &observed, &segmented})
    {
        require(simulation->SPHParticles().hostSize() == 22, "A short-buffer source must grow to its finite stroke count, without adding a trailing reservoir.");
        require(simulation->SPHJetsCompleted(), "All finite sources must complete at the visible time.");
        requireJetOwnership(*simulation, 1, 12, 1.0, 0.006);
        requireJetOwnership(*simulation, 5, 8, 0.5, -0.006);
        std::vector<int> ownership(simulation->SPHParticles().hostSize(), 0);
        for (int firstIndex : {1, 5})
            for (int index : simulation->SPHJetParticleIndices(firstIndex))
                ++ownership[index];
        for (int index = 0; index < static_cast<int>(ownership.size()); ++index)
            require(ownership[index] == (index == 0 || index == 9 ? 0 : 1), "A source claimed another source's or an ordinary block's particle.");
        require(math::nearlyEqual(simulation->SPHParticles().host()[0].position(), {0.0205 + 0.125 * simulation->time(), 0.0205, 0.0205}),
                "Dynamic device growth overwrote the ordinary block prefix.");
        require(math::nearlyEqual(simulation->SPHParticles().host()[9].position(), {-0.0195 - 0.125 * simulation->time(), 0.0205, 0.0205}),
                "Source range bookkeeping changed an ordinary block inserted between source setup and emission.");
    }
    segmented.solve(31);
    require(segmented.SPHParticles().hostSize() == 22, "A completed source resumed generating particles.");
    std::filesystem::remove_all(outputDirectory);
}

void testReleasedParticleCannotReenterConstraint(executionMode mode)
{
    testSPHDEM simulation{mode};
    configureGrowingJets(simulation);
    simulation.solve(60);
    require(simulation.SPHParticles().host()[1].position().x > 0.0 && !simulation.SPHParticles().host()[1].isConstrained(),
            "The leading seed must cross the outlet and be released while the source is still active.");
    simulation.repositionParticle(1, {-0.00025, 0.006, 0.0});
    simulation.solve(1);
    require(!simulation.SPHParticles().host()[1].isConstrained(), "A released particle was constrained again after returning into the inlet pipe.");
}

void testNonGridBirthsAfterTimeStepChange(executionMode mode)
{
    constexpr math::Real spacing = 1.0e-3;
    constexpr math::Real speed = 0.7;
    testSPHDEM simulation{mode};
    simulation.setBoundary({-0.03, -0.03, -0.03}, {0.03, 0.03, 0.03});
    simulation.setGravity(math::Vec3::zero());
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(spacing, 1.3 * spacing, 1000.0, 1.0e-3);
    simulation.setSPHMaximumVelocity(1.0);
    require(simulation.addSPHJet(SPHJet{math::Vec3::zero(), math::Vec3::unitX(), 0.5 * spacing, speed, 0.012}) == 0,
            "The non-grid source must retain its initial stable handle.");
    require(simulation.SPHJetParticleIndices(0).size() == 4, "A non-grid source must initially contain only its four buffer layers.");
    simulation.solve(99);
    require(simulation.SPHParticles().hostSize() == 4, "A source emitted a layer before its non-grid birth time.");
    simulation.setTimeStep(1.5e-5);
    // First birth = spacing/speed = 0.001428571..., between new grid nodes
    // 0.001425 and 0.001440. It may be appended only when the latter starts
    // the next acoustic interval, not while observing the former.
    simulation.solve(29);
    simulation.observeCurrentState([&](const solver&)
    {
        require(simulation.SPHParticles().hostSize() == 4 && simulation.SPHJetParticleIndices(0).size() == 4,
                "Changing DEM time step or observing a non-grid birth exposed a future layer.");
    });
    simulation.solve(2);
    require(simulation.SPHParticles().hostSize() == 5 && simulation.SPHJetParticleIndices(0).size() == 5,
            "The first non-grid layer must be emitted exactly once at the next common DEM node.");
    const auto& firstNew = simulation.SPHParticles().host()[4];
    require(math::nearlyEqual(firstNew.position(), {speed * simulation.time() - 4.5 * spacing, 0.0, 0.0}),
            "A delayed layer did not compensate its insertion position for the non-grid birth phase.");
    simulation.solve(0);
    simulation.solve(0);
    require(simulation.SPHJetParticleIndices(0).size() == 5, "Repeated pause completion duplicated a non-grid source layer.");
    simulation.solve(740 - 31);
    require(simulation.SPHParticles().hostSize() == 8 && simulation.SPHJetsCompleted(),
            "A changed DEM time step changed the finite source supply or completion time.");
    requireJetOwnership(simulation, 0, 8, speed, 0.0);
    std::vector<int> ownership(8, 0);
    for (int index : simulation.SPHJetParticleIndices(0))
    {
        require(index >= 0 && index < 8, "A non-grid source returned an invalid appended index.");
        ++ownership[index];
    }
    for (int count : ownership)
        require(count == 1, "Non-grid source timing created a duplicate or skipped a stable particle index.");
}

} // namespace

int main()
{
    try
    {
        testValueType();
        testPipeGeometry();
        testDiscretization();
        testOpeningGeometry();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testDeviceOpeningMask();
#endif
        const auto verifyMode = [](executionMode mode)
        {
            testGenerationConstraintAndCompletion(mode);
            testMultipleJetSnapshotContinuation(mode);
            testGrowingJetOwnershipAndContinuation(mode);
            testReleasedParticleCannotReenterConstraint(mode);
            testNonGridBirthsAfterTimeStepChange(mode);
            testVirtualOpeningAndContinuation(mode, false);
            testVirtualOpeningAndContinuation(mode, true);
        };
        verifyMode(executionMode::CPU);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        verifyMode(executionMode::GPU);
        verifyMode(executionMode::Hybrid);
#endif
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
