/** @file tutorial3.cpp @brief A small SPH block released inside a fixed tank. */
#include "solver/solvers.h"
#include <cmath> // Standard C++ square root for the velocity estimate below.

int main()
{
    using namespace fundem;
    using Real = math::Real;
    using Vec3 = math::Vec3;

    // 1. Choose fluid resolution, tank dimensions and run duration in SI units.
    constexpr Real spacing = 0.02, smoothingLength = 1.3 * spacing;
    constexpr Real gravity = 9.81, timeStep = 2.0e-5;
    constexpr int steps = 10000;
    constexpr int3 fluidCount{12, 8, 10};
    constexpr Vec3 tankSize{0.50, 0.20, 0.30};

    // 2. Set spacing, smoothing length, density and dynamic viscosity.
    // The velocity estimate controls the automatic sound-speed choice.
    SPHDEM simulation(executionMode::CPU);
    simulation.setSPHProperties(spacing, smoothingLength, 1000.0, 1.0e-3);
    simulation.setSPHMaximumVelocity(std::sqrt(2.0 * gravity * fluidCount.z * spacing));
    simulation.addSPHBlock(Vec3::zero(), fluidCount); // 960 resting particles on a cell-centered lattice.

    // 3. Build a closed box and reverse its SDF: solid boundary outside, fluid
    // cavity inside. Translate the centered geometry by half the tank size.
    const int materialIndex = simulation.addMaterial(material{1.0e9, 0.25, 1.0, 1000.0});
    levelset::BoxWall tank{tankSize};
    tank.buildLSGrid(spacing, 3, true);
    tank.reverseSDFSign();
    const int geometryIndex = simulation.addGeometry(tank);
    LSParticle wall;
    wall.setPosition(0.5 * tankSize);
    wall.setMaterial(simulation.materials(), materialIndex);
    wall.setGeometry(simulation.geometries(), geometryIndex);
    wall.setInfiniteMass();
    simulation.addLSParticle(wall);

    // 4. Configure the search region, gravity and DEM time step. The solver
    // chooses SPH acoustic/advection intervals as integer DEM-step multiples.
    simulation.setBoundary({-0.08, -0.08, -0.08}, tankSize + Vec3{0.08});
    simulation.setGravity({0.0, 0.0, -gravity});
    simulation.setTimeStep(timeStep);
    simulation.setOutputDirectory("tutorial3_files");
    simulation.setOutputStepInterval(500); // 0.01 s between frames.
    simulation.setVTUOutputEnabled(true);
    simulation.setVTUOutputFormat(vtuFormat::binary);
    simulation.solve(steps); // 0.2 s, 21 frames. energy.dat covers solid DEM only.
    return 0;
}
