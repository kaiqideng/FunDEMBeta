/** @file tutorial0.cpp @brief One sphere rebounding from a fixed plane. */
#include "solver/solvers.h"

int main()
{
    using namespace fundem;
    using math::Vec3;

    // 1. Choose the backend. GPU/Hybrid need a CUDA-enabled core.
    SphereDEM simulation(executionMode::CPU);

    // 2. material(E in Pa, Poisson ratio, restitution, density in kg/m^3).
    const int materialIndex = simulation.addMaterial(material{1.0e7, 0.25, 0.8, 2500.0});

    // 3. Build the wall's SDF grid, then register the reusable geometry.
    levelset::PlaneWall plane{Vec3::unitZ(), 0.4};
    plane.buildLSGrid(0.02, 2);
    const int planeIndex = simulation.addGeometry(plane);

    // 4. Set radius/material before adding the sphere to determine mass/inertia.
    particle sphere;
    sphere.setRadius(0.02);
    sphere.setMaterial(simulation.materials(), materialIndex);
    sphere.setPosition({0.0, 0.0, 0.06});
    simulation.addSphere(sphere);

    LSParticle wall;
    wall.setMaterial(simulation.materials(), materialIndex);
    wall.setGeometry(simulation.geometries(), planeIndex);
    simulation.addLSParticle(wall); // PlaneWall has zero volume, so the wall is fixed.

    // 5. Configure the material pair; stiffness/damping remain automatic.
    materialPair contact{materialIndex, materialIndex};
    contact.setFrictionCoefficients(0.1);
    simulation.setMaterialPair(contact);

    // 6. Set search bounds (not walls), gravity, time step and output.
    simulation.setBoundary({-0.25, -0.25, -0.25}, {0.25, 0.25, 0.25});
    simulation.setGravity({0.0, 0.0, -9.81});
    simulation.setTimeStep(1.0e-4);
    simulation.setOutputDirectory("tutorial0_files"); // Dedicated folder beside the executable.
    simulation.setOutputStepInterval(100); // Save every 0.01 s.
    simulation.setVTUOutputEnabled(true);
    simulation.setVTUOutputFormat(vtuFormat::binary);
    simulation.solve(3000); // Initializes automatically; 0.3 s, 31 frames including t=0.
    return 0;
}
