/** @file dem06-spin-identical.cpp
 * @brief Direct native-core LAMMPS comparison: dem06-spin-identical.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem06-spin-identical-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem06-spin-identical");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem06-spin-identical", "LAMMPS", "PASS",
            9.9999999999999995e-08, 2400, {0.0, 0.0, 0.0},
            {
                {2400, "wy", "lower", "relative", 22.859999999999999, 1e-08},
                {2400, "vx", "lower", "absolute", 0.0, 1e-08},
                {2400, "vy", "lower", "absolute", 0.0, 1e-08},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{70000000000.0, 0.33000000000000002, 0.5, 2700.0});
        // Initial body states, in the original stable insertion order.
        particle lowerBody;
        lowerBody.setRadius(0.0025000000000000001);
        lowerBody.setMaterial(simulation.materials(), particleMaterial);
        lowerBody.setPosition({0.0, 0.0, -0.0025249999999999999});
        lowerBody.setVelocity({0.0, 0.0, 0.20000000000000001});
        lowerBody.setAngularVelocity({0.0, 22.859999999999999, 0.0});
        lowerBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int lowerIndex = simulation.addSphere(lowerBody);
        bodies.push_back({"lower", true, lowerIndex});

        particle upperBody;
        upperBody.setRadius(0.0025000000000000001);
        upperBody.setMaterial(simulation.materials(), particleMaterial);
        upperBody.setPosition({0.0, 0.0, 0.0025249999999999999});
        upperBody.setVelocity({0.0, 0.0, -0.20000000000000001});
        upperBody.setAngularVelocity({0.0, -22.859999999999999, 0.0});
        upperBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int upperIndex = simulation.addSphere(upperBody);
        bodies.push_back({"upper", true, upperIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{particleMaterial, particleMaterial};
        contact0.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(false);
        contact0.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.0050000000000000001, -0.0050000000000000001, -0.0075729999999999999}, {0.0050000000000000001, 0.0050000000000000001, 0.0075729999999999999});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem06-spin-identical: " << error.what() << '\n';
        return 1;
    }
}
