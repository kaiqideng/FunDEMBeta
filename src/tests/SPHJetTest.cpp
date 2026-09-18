#include "execution/sphJetFunctions.h"
#include "execution/sphSourceDiscretization.h"
#include "solver/SPHDEM.h"
#include "particle/SPHJet.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

using namespace fundem;

class testSPHDEM : public SPHDEM
{
public:
    using SPHDEM::SPHDEM;
    SPHParticleContainer& particles() noexcept { return mutableSPHParticles(); }
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
    require(nearlyEqual(jet.pipeLength(), 3.0), "SPHJet pipe length must equal speed times duration.");

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
    require(nearlyEqual(simulation.SPHParticles().host()[0].density(), 1100.0), "A constrained jet particle must skip density reinitialization and integration.");
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

void testMultipleJetSnapshotContinuation()
{
    constexpr math::Real spacing = 1.0e-3;
    constexpr math::Real timeStep = 1.0e-5;
    constexpr math::Real lastEndTime = 156.5 * timeStep;
    testSPHDEM reference;
    testSPHDEM observed;
    testSPHDEM split;
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
        // Keep the tail at the upstream pipe boundary so its constraint remains
        // active in the deferred state immediately before the final jet ends.
        simulation->particles().host()[0].setPosition({-lastEndTime, 0.004, 0.0});
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
            require(observed.SPHParticles().host()[0].isConstrained(), "The final jet must still be constrained in the deferred state at the snapshot boundary.");
            observed.observeCurrentState([&](const solver&)
            {
                require(observed.SPHJetsCompleted() == (step == 157), "A current-state snapshot must report the visible jet completion time.");
                require(observed.SPHParticles().host()[0].isConstrained() == (step != 157), "A snapshot crossing the final jet end must clear its constraint.");
            });
            requireSameState(reference, observed);
            observed.writeOutput();
        }
        requireSameState(reference, observed);
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
    requireSameState(reference, observed);
    requireSameState(reference, split);

    reference.solve(19);
    observed.solve(19);
    split.solve(19);
    requireSameState(reference, observed);
    requireSameState(reference, split);
    require(split.SPHJetsCompleted(), "All jets must remain complete after further solve continuation.");
    for (const SPHParticle& particle : split.SPHParticles().host())
        require(!particle.isConstrained(), "No completed jet may retain a particle constraint after continuation.");
    std::filesystem::remove_all(outputDirectory);
}

} // namespace

int main()
{
    try
    {
        testValueType();
        testPipeGeometry();
        testDiscretization();
        testGenerationConstraintAndCompletion(executionMode::CPU);
        testMultipleJetSnapshotContinuation();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testGenerationConstraintAndCompletion(executionMode::GPU);
#endif
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
