/** @file dem15-hooke-oblique.cpp
 * @brief Direct native-core LAMMPS comparison: dem15-hooke-oblique.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem15-hooke-oblique-3d-si.yaml
 * Historical PASS is provenance; current acceptance uses only the declared checks.
 */
#include "comparisonCommon.h"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;
        const auto options = comparison::parseOptions(argc, argv, "dem15-hooke-oblique");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem15-hooke-oblique", "LAMMPS", "PASS",
            1e-08, 7000, {0.0, 0.0, 0.0},
            {
                {7000, "vx", "left", "relative", -0.69999999999999996, 0.01},
                {7000, "vy", "left", "relative", 1.575, 0.01},
                {7000, "wz", "left", "relative", -425.0, 0.01},
                {7000, "vx", "right", "relative", 0.69999999999999996, 0.01},
                {7000, "vy", "right", "relative", -1.575, 0.01},
                {7000, "wz", "right", "relative", -425.0, 0.01},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{10000.0, 0.25, 0.69999999999999996, 2700.0});
        // Initial body states, in the original stable insertion order.
        particle leftBody;
        leftBody.setRadius(0.0025000000000000001);
        leftBody.setMaterial(simulation.materials(), particleMaterial);
        leftBody.setPosition({-0.0025500000000000002, -0.0001, 0.0});
        leftBody.setVelocity({1.0, 2.0, 0.0});
        leftBody.setAngularVelocity({0.0, 0.0, 0.0});
        leftBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int leftIndex = simulation.addSphere(leftBody);
        bodies.push_back({"left", true, leftIndex});

        particle rightBody;
        rightBody.setRadius(0.0025000000000000001);
        rightBody.setMaterial(simulation.materials(), particleMaterial);
        rightBody.setPosition({0.0025500000000000002, 0.0001, 0.0});
        rightBody.setVelocity({-1.0, -2.0, 0.0});
        rightBody.setAngularVelocity({0.0, 0.0, 0.0});
        rightBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int rightIndex = simulation.addSphere(rightBody);
        bodies.push_back({"right", true, rightIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{particleMaterial, particleMaterial};
        contact0.setLinearStiffness(300000000.0, 80000000.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(true);
        contact0.setFrictionCoefficients(0.25, 0.0, 0.0);
        contact0.setDampingCoefficients(36.732751400992463, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(true);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.00762, -0.0052399999999999999, -0.0050000000000000001}, {0.00762, 0.0052399999999999999, 0.0050000000000000001});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem15-hooke-oblique: " << error.what() << '\n';
        return 1;
    }
}
