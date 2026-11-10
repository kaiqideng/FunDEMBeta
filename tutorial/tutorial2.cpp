/** @file tutorial2.cpp @brief A five-sphere bonded beam bending under gravity. */
#include "solver/solvers.h"

int main()
{
    using namespace fundem;
    using Real = math::Real;

    // 1. Choose dimensions/material properties in SI units and run settings.
    constexpr int particleCount = 5;
    constexpr Real radius = 0.01, pitch = 0.025, density = 1000.0;
    constexpr Real young = 5.0e5, poisson = 0.3, timeStep = 2.0e-5;
    constexpr int steps = 5000;
    SphereDEM simulation(executionMode::CPU);
    const int materialIndex = simulation.addMaterial(material{young, poisson, 1.0, density});

    // 2. Disable ordinary sphere contact stiffness to isolate the bonds.
    materialPair contacts{materialIndex, materialIndex};
    contacts.setLinearStiffness(0.0, 0.0);
    simulation.setMaterialPair(contacts);

    // 3. Add five spheres; infinite mass fixes the first one as the beam root.
    for (int index = 0; index < particleCount; ++index)
    {
        particle value;
        value.setPosition({pitch * index, 0.0, 0.0});
        value.setRadius(radius);
        value.setMaterial(simulation.materials(), materialIndex);
        if (index == 0)
            value.setInfiniteMass();
        simulation.addSphere(value);
    }

    // 4. Connect adjacent particle indices. Circular-section axial/shear
    // stiffnesses use N/m; bending/torsional stiffnesses use N m.
    constexpr Real area = math::pi * radius * radius;
    constexpr Real secondMoment = 0.25 * area * radius * radius;
    constexpr Real shearModulus = young / (2.0 * (1.0 + poisson));
    for (int index = 0; index + 1 < particleCount; ++index)
    {
        bond connection{pitch};
        connection.setCrossSectionArea(0.0); // Disable fracture, not elastic stiffness.
        if (!connection.setConnection(simulation.spheres(), index, index + 1) ||
            !connection.setStiffness(young * area / pitch,
                12.0 * young * secondMoment / (pitch * pitch * pitch),
                young * secondMoment / pitch, 2.0 * shearModulus * secondMoment / pitch))
            return 1; // Both setters report failure through their return values.
        simulation.addBond(connection);
    }

    // 5. Set search bounds, gravity and output; solve initializes automatically.
    simulation.setBoundary({-0.03, -0.03, -0.15}, {0.15, 0.03, 0.03});
    simulation.setGravity({0.0, 0.0, -9.81});
    simulation.setTimeStep(timeStep);
    simulation.setOutputDirectory("tutorial2_files");
    simulation.setOutputStepInterval(250); // 0.005 s between frames.
    simulation.setVTUOutputEnabled(true);
    simulation.setVTUOutputFormat(vtuFormat::binary);
    simulation.solve(steps); // 0.1 s, 21 frames; an illustrative beam transient.
    return 0;
}
