#include "execution/sphFunctions.h"
#include "execution/sphSourceDiscretization.h"
#include "solver/SPHDEM.h"
#include "solver/SPHNeighborhood.h"
#include "solver/sphInteraction.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
using namespace fundem;
void require(bool value, const char* message)
{
    if (!value)
        throw std::runtime_error(message);
}

void testComputedDensityFailure()
{
    SPHParticleContainer particles;
    virtualParticleContainer boundaries;
    particles.host().push_back(SPHParticle{{0.0, 0.0, 0.0}});
    particles.host().push_back(SPHParticle{{0.01, 0.0, 0.0}});
    for (auto& particle : particles.host())
    {
        particle.setMass(1.0);
        particle.setDensity(1000.0);
        particle.setPressure(0.0);
    }
    cpu::SPHInteraction interaction;
    const cpu::SPHInteractionParameters parameters{0.03, 1.0, 1000.0, 1.0, 10.0, 0.001, math::Vec3::zero()};
    interaction.buildNeighborhood(particles, boundaries, {-0.2, -0.2, -0.2}, {0.2, 0.2, 0.2}, 0.06, 0.08);
    particles.host()[1].*SPHParticle::pressureField::member = std::numeric_limits<math::Real>::quiet_NaN();
    interaction.updateDensity(particles, boundaries, parameters, 1.0e-5);
    require(!math::isFinite(particles.host()[0].densityRate()), "A computed invalid rate was replaced by old CPU state.");
    require(interaction.stateStatistics(particles, math::Vec3::zero()).invalidValueCount_ > 0, "Invalid density rate escaped existing diagnostics.");
    require(particles.host()[0].density() == 100.0, "CPU emergency density handling differs from the shared device formula.");
    for (auto& particle : particles.host())
    {
        particle.setPressure(0.0);
        particle.setDensity(1000.0);
    }
    interaction.updateDensity(particles, boundaries, parameters, 1.0e-5);
    require(!math::isFinite(particles.host()[0].densityRate()), "A later finite split stage erased the earlier failure.");
    particles.host()[0].setConstrained(true);
    interaction.updateDensity(particles, boundaries, parameters, 1.0e-5);
    require(!math::isFinite(particles.host()[0].densityRate()), "An inlet constraint erased the earlier failure.");
    interaction.initializeParticles(particles, 1000.0, 10.0);
    require(particles.host()[0].densityRate() == 0.0, "Initialization did not clear the diagnostic latch.");
}

void testSearchSkin()
{
    SPHParticleContainer particles;
    virtualParticleContainer boundaries;
    particles.host().push_back(SPHParticle{{0.0, 0.0, 0.0}});
    particles.host().push_back(SPHParticle{{0.2, 0.0, 0.0}});
    for (auto& particle : particles.host())
    {
        particle.setMass(1.0);
        particle.setDensity(1010.0);
        particle.setPressure(1000.0);
    }
    SPHNeighborhood neighborhood;
    neighborhood.captureHost(particles, boundaries);
    require(!neighborhood.needsRebuild(particles, boundaries, 0.02), "Unmoved particles should reuse neighbors.");
    particles.host()[1].setPosition({0.191, 0.0, 0.0});
    require(!neighborhood.needsRebuild(particles, boundaries, 0.02), "Displacement inside half the skin should reuse neighbors.");
    particles.host()[1].setPosition({0.025, 0.0, 0.0});
    require(neighborhood.needsRebuild(particles, boundaries, 0.02), "An approaching pair crossing cells did not expire the search skin.");
    cpu::SPHInteraction interaction;
    interaction.buildNeighborhood(particles, boundaries, {-0.3, -0.3, -0.3}, {0.3, 0.3, 0.3}, 0.06, 0.08);
    neighborhood.captureHost(particles, boundaries);
    const cpu::SPHInteractionParameters parameters{0.03, 1.0, 1000.0, 1.0, 10.0, 0.001, math::Vec3::zero()};
    interaction.updatePressureForceAndBoundaryForce(particles, boundaries, parameters);
    require(math::norm(particles.host()[0].force()) > 0.0, "The approaching neighbor was not recovered by rebuilding.");
    require(math::norm(particles.host()[0].force() + particles.host()[1].force()) < math::defaultTolerance, "Rebuilt fluid pair lost force symmetry.");
    boundaries.host().emplace_back();
    neighborhood.captureHost(particles, boundaries);
    boundaries.host()[0].setPosition({0.011, 0.0, 0.0});
    require(neighborhood.needsRebuild(particles, boundaries, 0.02), "Moving boundary displacement must also expire the skin.");
    neighborhood.invalidate();
    require(neighborhood.needsRebuild(particles, boundaries, 0.02), "Restored observation must invalidate search references.");
}

void testExtremeConfiguration()
{
    execution::SPHJetDiscretization lattice;
    require(!execution::trySPHJetDiscretization(0.001, 0.0005, 3.0e6, lattice), "Oversized axial lattice was converted to int.");
    require(!execution::trySPHJetDiscretization(0.001, 1.0e9, 0.001, lattice), "Oversized radial lattice entered enumeration.");
    require(!execution::trySPHJetDiscretization(1.0, 20000.5, 2.0, lattice), "An oversized circular lattice product was accepted.");
    require(execution::trySPHJetDiscretization(0.001, 0.0005, 0.003, lattice) && lattice.particleCount_ == 3, "Ordinary jet discretization changed.");
    require(execution::trySPHJetDiscretization(0.5, 3.2, 1.0, lattice), "Small circular lattice was rejected.");
    int expectedCount = 0;
    for (int y = -lattice.radialIndex_; y <= lattice.radialIndex_; ++y)
        for (int x = -lattice.radialIndex_; x <= lattice.radialIndex_; ++x)
            expectedCount += lattice.contains(x, y) ? lattice.axialCount_ : 0;
    require(lattice.particleCount_ == expectedCount, "Fast circular row counting disagrees with particle generation.");
    SPHDEM solver;
    solver.setBoundary({-0.1, -0.1, -0.1}, {0.1, 0.1, 0.1});
    solver.setGravity(math::Vec3::zero());
    solver.setTimeStep(1.0e-13);
    solver.setSPHProperties(0.02, 0.03, 1000.0, 0.001);
    solver.addSPHParticle(SPHParticle{{0.0, 0.0, 0.0}});
    solver.initialize();
    require(solver.SPHStepInterval() == std::numeric_limits<int>::max(), "Large acoustic ratio was not safely capped.");
    require(solver.SPHTimeStep() <= solver.SPHTimeStepLimit(), "Interval capping exceeded the acoustic limit.");
    SPHDEM unsafeSolver;
    unsafeSolver.setBoundary({-0.1, -0.1, -0.1}, {0.1, 0.1, 0.1});
    unsafeSolver.setGravity(math::Vec3::zero());
    unsafeSolver.setTimeStep(1.0e-12);
    unsafeSolver.setSPHProperties(0.02, 0.03, 1000.0, 0.001);
    unsafeSolver.setSPHSoundSpeed(1.0e12);
    unsafeSolver.addSPHParticle(SPHParticle{{0.0, 0.0, 0.0}});
    bool rejected = false;
    try
    {
        unsafeSolver.initialize();
    }
    catch (const std::runtime_error&)
    {
        rejected = true;
    }
    require(rejected, "An absolute time tolerance admitted a DEM step far above the acoustic limit.");
}
} // namespace

int main()
{
    try
    {
        testComputedDensityFailure();
        testSearchSkin();
        testExtremeConfiguration();
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
