/** @file dem02-hertz-peak.cpp
 * @brief Direct native-core LAMMPS comparison: dem02-hertz-peak.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem02-hertz-peak-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem02-hertz-peak");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem02-hertz-peak", "LAMMPS", "PASS",
            2e-08, 2508, {0.0, 0.0, 0.0},
            {
                {1508, "overlap", "left", "relative", 0.00027411268366588467, 0.01},
                {1508, "normal_force", "left", "relative", 10696.88970293861, 0.01},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{46874953125.0, 0.25, 1.0, 2800.0});
        // Initial body states, in the original stable insertion order.
        particle leftBody;
        leftBody.setRadius(0.01);
        leftBody.setMaterial(simulation.materials(), particleMaterial);
        leftBody.setPosition({-0.0101, 0.0, 0.0});
        leftBody.setVelocity({10.0, 0.0, 0.0});
        leftBody.setAngularVelocity({0.0, 0.0, 0.0});
        leftBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int leftIndex = simulation.addSphere(leftBody);
        bodies.push_back({"left", true, leftIndex});

        particle rightBody;
        rightBody.setRadius(0.01);
        rightBody.setMaterial(simulation.materials(), particleMaterial);
        rightBody.setPosition({0.0101, 0.0, 0.0});
        rightBody.setVelocity({-10.0, 0.0, 0.0});
        rightBody.setAngularVelocity({0.0, 0.0, 0.0});
        rightBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int rightIndex = simulation.addSphere(rightBody);
        bodies.push_back({"right", true, rightIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{particleMaterial, particleMaterial};
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
        std::cerr << "dem02-hertz-peak: " << error.what() << '\n';
        return 1;
    }
}
