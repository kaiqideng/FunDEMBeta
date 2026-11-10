// One FunDEM include supplies solvers, materials, particles and output types.
#include "solver/solvers.h"

// Standard C++ utilities for creating result folders and printing their paths.
#include <filesystem>
#include <iostream>
#include <string>

int main()
{
    using namespace fundem;

    // 1. Choose model and run parameters in SI units: m, kg, s, Pa.
    // Restitution, Poisson ratio and friction are dimensionless.
    const math::Real radius = 0.02, density = 2500.0;
    const math::Real young = 1e9, poisson = 0.25, restitution = 0.8;
    const math::Real friction = 0.0, speed = 0.5, halfSpacing = 0.021;
    const math::Real timeStep = 1e-6;
    const int steps = 6000, outputEvery = 100;

    // 2. SphereDEM owns the model and advances sphere/sphere and sphere/LS contacts.
    // Choose CPU, GPU or Hybrid here; GPU/Hybrid require a CUDA-enabled core.
    SphereDEM simulation(executionMode::CPU);

    // 3. material(E, Poisson ratio, restitution, density) defines bulk properties.
    // addMaterial() stores it and returns the index used by particles and pairs.
    const int materialIndex = simulation.addMaterial(material{young, poisson, restitution, density});

    // 4. materialPair selects contact settings for two material indices.
    // The same index twice configures contacts within this material.
    materialPair contact{materialIndex, materialIndex};
    contact.setFrictionCoefficients(friction); // Sliding; optional rolling/torsional values default to zero.
    // Leave stiffness/damping setters out for Hertz-Mindlin stiffness and
    // restitution-derived damping. Explicit zero damping would disable dissipation.
    simulation.setMaterialPair(contact); // Stores/replaces the complete pair settings.

    // 5. particle describes one sphere. Set radius and material first: these
    // determine its mass and inertia. Position is its center; velocity is in m/s.
    for (const double side : {-1.0, 1.0})
    {
        particle sphere;
        sphere.setRadius(radius);
        sphere.setMaterial(simulation.materials(), materialIndex);
        sphere.setPosition({side * halfSpacing, 0, 0});
        sphere.setVelocity({-side * speed, 0, 0});
        simulation.addSphere(sphere); // Copies the sphere; returns an index for later edits.
    }

    // 6. Define contact-search bounds, gravity (m/s^2), and the time step (s).
    // setBoundary() creates no physical walls; add wall geometry separately.
    simulation.setBoundary({-0.08, -0.04, -0.04}, {0.08, 0.04, 0.04});
    simulation.setGravity({0, 0, 0});
    simulation.setTimeStep(timeStep);

    // 7. Reserve a fresh example_result/run_N beside the executable.
    // outputDirectory() returns the resolved path; create_directory() reserves
    // it only if this call creates it, so earlier/concurrent results are safe.
    std::filesystem::path output;
    int run = 0;
    do
    {
        simulation.setOutputDirectory("example_result/run_" + std::to_string(++run));
        output = simulation.outputDirectory();
        std::filesystem::create_directories(output.parent_path());
    } while (!std::filesystem::create_directory(output));

    // Save particle/contact VTU files and energy.dat every outputEvery steps.
    // Binary VTU is compact and can be opened in ParaView.
    simulation.setOutputStepInterval(outputEvery);
    simulation.setVTUOutputEnabled(true);
    simulation.setVTUOutputFormat(vtuFormat::binary);
    std::cout << "Output: " << output << '\n';

    // 8. solve() initializes automatically, then advances this many steps.
    // Here: 0.006 s, 61 scheduled frames including t=0. Calling solve() again
    // continues from the current time; no separate initialize() is needed.
    simulation.solve(steps);
    return 0;
}
