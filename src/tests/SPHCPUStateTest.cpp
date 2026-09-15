#include "solver/SPHDEM.h"

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

void configure(SPHDEM& solver)
{
    solver.setBoundary({-0.2, -0.2, -0.2}, {0.2, 0.2, 0.2});
    solver.setGravity(math::Vec3::zero());
    solver.setTimeStep(1.0e-4);
    solver.setSPHProperties(0.02, 0.03, 1000.0, 1.0e-3);
    solver.setSPHSoundSpeed(10.0);
    solver.addSPHParticle(SPHParticle{{-0.015, 0.0, 0.0}, {0.1, 0.0, 0.0}});
    solver.addSPHParticle(SPHParticle{{0.015, 0.0, 0.0}, {-0.1, 0.0, 0.0}});
}

void configureSingleParticle(SPHDEM& solver, const SPHParticle& particle, math::Real timeStep, const math::Vec3& gravity, math::Real maximumVelocity, math::Real soundSpeed)
{
    solver.setBoundary({-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3});
    solver.setGravity(gravity);
    solver.setTimeStep(timeStep);
    solver.setSPHProperties(0.02, 0.03, 1000.0, 1.0e-3);
    solver.setSPHMaximumVelocity(maximumVelocity);
    solver.setSPHSoundSpeed(soundSpeed);
    solver.addSPHParticle(particle);
}

bool nearlyEqual(math::Real first, math::Real second, math::Real tolerance = 1.0e-11) noexcept { return std::abs(first - second) <= tolerance * (1.0 + std::max(std::abs(first), std::abs(second))); }

void requireEqual(const SPHDEM& first, const SPHDEM& second, const char* message)
{
    const auto& firstParticles = first.SPHParticles().host();
    const auto& secondParticles = second.SPHParticles().host();
    if (firstParticles.size() != secondParticles.size())
    {
        throw std::runtime_error(message);
    }
    for (int index = 0; index < static_cast<int>(firstParticles.size()); ++index)
    {
        const SPHParticle& a = firstParticles[index];
        const SPHParticle& b = secondParticles[index];
        if (!nearlyEqual(a.position().x, b.position().x) || !nearlyEqual(a.position().y, b.position().y) || !nearlyEqual(a.position().z, b.position().z) ||
            !nearlyEqual(a.velocity().x, b.velocity().x) || !nearlyEqual(a.velocity().y, b.velocity().y) || !nearlyEqual(a.velocity().z, b.velocity().z) || !nearlyEqual(a.density(), b.density()))
        {
            throw std::runtime_error(message);
        }
    }
}

class observedOutputSPH final : public SPHDEM
{
public:
    math::Real outputPosition_{0.0};

protected:
    void writeSystemVTU(int frameIndex) override
    {
        outputPosition_ = SPHParticles().host()[0].position().x;
        SPHDEM::writeSystemVTU(frameIndex);
    }
};

void verifyCurrentStateObservation(const std::filesystem::path& outputDirectory)
{
    observedOutputSPH observed;
    SPHDEM reference;
    configureSingleParticle(observed, SPHParticle{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, 1.0e-5, math::Vec3::zero(), 1.0, 10.0);
    configureSingleParticle(reference, SPHParticle{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, 1.0e-5, math::Vec3::zero(), 1.0, 10.0);
    observed.setOutputDirectory((outputDirectory / "observation").string());
    observed.solve(0);
    observed.initialize();
    reference.initialize();
    const int initialInterval = observed.SPHStepInterval();
    if (initialInterval <= 1)
        throw std::runtime_error("The observation regression requires a deferred acoustic step.");

    for (int step = 1; step <= 2 * observed.SPHAdvectionStepInterval() + initialInterval; ++step)
    {
        observed.step();
        reference.step();
        const int intervalBefore = observed.SPHStepInterval();
        observed.observeCurrentState([&](const solver& state)
        {
            const auto& particles = static_cast<const SPHDEM&>(state).SPHParticles().host();
            if (state.stepCount() != step || !nearlyEqual(particles[0].position().x, state.time()))
                throw std::runtime_error("A current-time SPH observation exposed deferred particle positions.");
            observed.observeCurrentState([&](const solver& nested)
            {
                if (&state != &nested || !nearlyEqual(observed.SPHParticles().host()[0].position().x, nested.time()))
                    throw std::runtime_error("A nested SPH observation lost the outer snapshot.");
            });
            if (step == 1)
            {
                int rejectedLifecycleChanges = 0;
                for (const auto& mutation : {std::function<void()>{[&] { observed.initialize(); }},
                                             std::function<void()>{[&] { observed.step(); }},
                                             std::function<void()>{[&] { observed.solve(0); }}})
                {
                    try
                    {
                        mutation();
                    }
                    catch (const std::logic_error&)
                    {
                        ++rejectedLifecycleChanges;
                    }
                }
                if (rejectedLifecycleChanges != 3)
                    throw std::runtime_error("An SPH observer was allowed to change the integration lifecycle.");
                observed.writeOutput();
                if (!nearlyEqual(observed.outputPosition_, state.time()) || !nearlyEqual(particles[0].position().x, state.time()))
                    throw std::runtime_error("VTU and nested host observation did not share the current SPH state.");
            }
        });
        if (observed.SPHStepInterval() != intervalBefore)
            throw std::runtime_error("Observing SPH changed the deferred acoustic schedule.");
        requireEqual(reference, observed, "Observing SPH changed its deferred continuation state.");

        bool caught = false;
        try
        {
            observed.observeCurrentState([](const solver&) { throw std::runtime_error("observer failure"); });
        }
        catch (const std::runtime_error& error)
        {
            caught = std::string(error.what()) == "observer failure";
        }
        if (!caught)
            throw std::runtime_error("An SPH observation did not propagate its callback exception.");
        requireEqual(reference, observed, "A throwing SPH observation did not restore continuation state.");
    }

    reference.solve(7);
    observed.solve(7);
    observed.observeCurrentState([&](const solver&) { requireEqual(reference, observed, "Observing completed solve() changed its current host state."); });
    requireEqual(reference, observed, "Observing completed solve() discarded its current host state.");
    reference.solve(11);
    observed.solve(11);
    requireEqual(reference, observed, "Observation after solve() changed subsequent continuation.");

    SPHDEM interactingReference;
    SPHDEM interactingObserved;
    configure(interactingReference);
    configure(interactingObserved);
    for (int step = 0; step < 125; ++step)
    {
        interactingReference.step();
        interactingObserved.step();
        interactingObserved.observeCurrentState([](const solver&) {});
        requireEqual(interactingReference, interactingObserved, "SPH observation changed interacting-particle continuation across an advection boundary.");
    }
    interactingReference.solve(0);
    interactingObserved.solve(0);
    requireEqual(interactingReference, interactingObserved, "Observed interacting SPH particles ended at a different current state.");
}

} // namespace

int main()
{
    namespace fs = std::filesystem;
    try
    {
        SPHDEM uninterrupted;
        SPHDEM withOutput;
        SPHDEM split;
        SPHDEM manualOutput;
        configure(uninterrupted);
        configure(withOutput);
        configure(split);
        configure(manualOutput);

        const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        const fs::path outputDirectory = fs::temp_directory_path() / ("fundem-sph-state-" + std::to_string(stamp));
        withOutput.setOutputDirectory(outputDirectory.string());
        withOutput.setOutputStepInterval(3);
        manualOutput.setOutputDirectory((outputDirectory / "manual").string());

        verifyCurrentStateObservation(outputDirectory);

        uninterrupted.solve(25);
        withOutput.solve(25);
        split.solve(7);
        split.solve(18);
        manualOutput.solve(7);
        manualOutput.writeOutput();
        manualOutput.solve(18);

        requireEqual(uninterrupted, withOutput, "SPH output cadence changed the numerical state.");
        requireEqual(uninterrupted, split, "Splitting solve() changed the SPH numerical state.");
        requireEqual(uninterrupted, manualOutput, "Writing output between CPU SPH solve() calls changed the continuation state.");

        SPHDEM changedConfiguration;
        configureSingleParticle(changedConfiguration, SPHParticle{{0.0, 0.0, 0.0}, {0.2, 0.0, 0.0}}, 1.0e-4, math::Vec3::zero(), 0.5, 5.0);
        changedConfiguration.solve(7);
        const SPHParticle continuationParticle = changedConfiguration.SPHParticles().host()[0];

        SPHDEM restartedFromCurrentState;
        configureSingleParticle(restartedFromCurrentState, continuationParticle, 2.0e-4, {0.0, -1.0, 0.0}, 0.7, 8.0);
        restartedFromCurrentState.setBoundary({-0.4, -0.4, -0.4}, {0.4, 0.4, 0.4});
        changedConfiguration.setBoundary({-0.4, -0.4, -0.4}, {0.4, 0.4, 0.4});
        changedConfiguration.setGravity({0.0, -1.0, 0.0});
        changedConfiguration.setTimeStep(2.0e-4);
        changedConfiguration.setSPHMaximumVelocity(0.7);
        changedConfiguration.setSPHSoundSpeed(8.0);
        changedConfiguration.solve(11);
        restartedFromCurrentState.solve(11);
        requireEqual(changedConfiguration, restartedFromCurrentState, "Changing solver or SPH configuration after a partial SPH step reused stale deferred state.");

        SPHDEM stationary;
        configureSingleParticle(stationary, SPHParticle{{0.0, 0.0, 0.0}}, 1.0e-4, math::Vec3::zero(), 1.0, 10.0);
        stationary.solve(25);
        const SPHParticle& stationaryParticle = stationary.SPHParticles().host()[0];
        if (stationaryParticle.position() != math::Vec3::zero() || stationaryParticle.velocity() != math::Vec3::zero() || !nearlyEqual(stationaryParticle.density(), 1000.0))
        {
            throw std::runtime_error("An isolated stationary SPH particle without gravity must remain stationary at its reference density.");
        }

        bool rejectedUnsafeTimeStep = false;
        try
        {
            SPHDEM unsafeTimeStep;
            configureSingleParticle(unsafeTimeStep, SPHParticle{{0.0, 0.0, 0.0}}, 1.0e-2, math::Vec3::zero(), 1.0, 10.0);
            unsafeTimeStep.solve(0);
        }
        catch (const std::runtime_error&)
        {
            rejectedUnsafeTimeStep = true;
        }
        if (!rejectedUnsafeTimeStep)
        {
            throw std::runtime_error("A DEM time step above the SPH stability limit must be rejected.");
        }
        fs::remove_all(outputDirectory);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
