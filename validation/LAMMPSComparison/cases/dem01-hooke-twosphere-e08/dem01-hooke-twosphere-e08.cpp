/** @file dem01-hooke-twosphere-e08.cpp
 * @brief Direct native-core LAMMPS comparison: dem01-hooke-twosphere-e08.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem01-hooke-twosphere-e08-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem01-hooke-twosphere-e08");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem01-hooke-twosphere-e08", "LAMMPS", "PASS",
            2.9999999999999997e-08, 2795, {0.0, 0.0, 0.0},
            {
                {2795, "momentum_x", "left", "absolute", 0.0, 1.4137166941154069e-05},
                {2795, "restitution", "left", "relative", 0.80000000000000004, 0.02},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{10000.0, 0.25, 0.80000000000000004, 2700.0});
        // Initial body states, in the original stable insertion order.
        particle leftBody;
        leftBody.setRadius(0.0025000000000000001);
        leftBody.setMaterial(simulation.materials(), particleMaterial);
        leftBody.setPosition({-0.0025500000000000002, 0.0, 0.0});
        leftBody.setVelocity({2.0, 0.0, 0.0});
        leftBody.setAngularVelocity({0.0, 0.0, 0.0});
        leftBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int leftIndex = simulation.addSphere(leftBody);
        bodies.push_back({"left", true, leftIndex});

        particle rightBody;
        rightBody.setRadius(0.0025000000000000001);
        rightBody.setMaterial(simulation.materials(), particleMaterial);
        rightBody.setPosition({0.0025500000000000002, 0.0, 0.0});
        rightBody.setVelocity({-2.0, 0.0, 0.0});
        rightBody.setAngularVelocity({0.0, 0.0, 0.0});
        rightBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int rightIndex = simulation.addSphere(rightBody);
        bodies.push_back({"right", true, rightIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{particleMaterial, particleMaterial};
        contact0.setLinearStiffness(1000000.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(true);
        contact0.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.0077177000000000009, -0.0050000000000000001, -0.0050000000000000001}, {0.0077177000000000009, 0.0050000000000000001, 0.0050000000000000001});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem01-hooke-twosphere-e08: " << error.what() << '\n';
        return 1;
    }
}
