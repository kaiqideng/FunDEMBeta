/** @file dem02-hertz-twosphere.cpp
 * @brief Direct native-core LAMMPS comparison: dem02-hertz-twosphere.
 * Reference: https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem02-hertz-twosphere-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem02-hertz-twosphere");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem02-hertz-twosphere", "LAMMPS", "PASS",
            2e-08, 2508, {0.0, 0.0, 0.0},
            {
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int sphereMaterial = simulation.addMaterial(material{49999950000.0, 0.0, 1.0, 2800.0});
        // Initial body states, in the original stable insertion order.
        particle leftBody;
        leftBody.setRadius(0.01);
        leftBody.setMaterial(simulation.materials(), sphereMaterial);
        leftBody.setPosition({-0.0101, 0.0, 0.0});
        leftBody.setVelocity({10.0, 0.0, 0.0});
        leftBody.setAngularVelocity({0.0, 0.0, 0.0});
        leftBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int leftIndex = simulation.addSphere(leftBody);
        bodies.push_back({"left", true, leftIndex});

        particle rightBody;
        rightBody.setRadius(0.01);
        rightBody.setMaterial(simulation.materials(), sphereMaterial);
        rightBody.setPosition({0.0101, 0.0, 0.0});
        rightBody.setVelocity({-10.0, 0.0, 0.0});
        rightBody.setAngularVelocity({0.0, 0.0, 0.0});
        rightBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int rightIndex = simulation.addSphere(rightBody);
        bodies.push_back({"right", true, rightIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{sphereMaterial, sphereMaterial};
        contact0.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(false);
        contact0.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.0306016, -0.02, -0.02}, {0.0306016, 0.02, 0.02});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem02-hertz-twosphere: " << error.what() << '\n';
        return 1;
    }
}
