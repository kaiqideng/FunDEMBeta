/** @file dem06-hookehistory-spin.cpp
 * @brief Direct native-core LAMMPS comparison: dem06-hookehistory-spin.
 * Reference: https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem06-hookehistory-spin-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem06-hookehistory-spin");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem06-hookehistory-spin", "LAMMPS", "PASS",
            9.9999999999999995e-07, 1800, {0.0, 0.0, 0.0},
            {
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int smallMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 0.5, 2700.0});
        const int largeMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 0.5, 2700000.0});
        // Initial body states, in the original stable insertion order.
        particle smallBody;
        smallBody.setRadius(0.10000000000000001);
        smallBody.setMaterial(simulation.materials(), smallMaterial);
        smallBody.setPosition({0.0, 0.0, 0.60001000000000004});
        smallBody.setVelocity({0.0, 0.0, -0.20000000000000001});
        smallBody.setAngularVelocity({0.0, 22.859999999999999, 0.0});
        smallBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int smallIndex = simulation.addSphere(smallBody);
        bodies.push_back({"small", true, smallIndex});

        particle largeBody;
        largeBody.setRadius(0.5);
        largeBody.setMaterial(simulation.materials(), largeMaterial);
        largeBody.setPosition({0.0, 0.0, 0.0});
        largeBody.setVelocity({0.0, 0.0, 0.0});
        largeBody.setAngularVelocity({0.0, 0.0, 0.0});
        largeBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int largeIndex = simulation.addSphere(largeBody);
        bodies.push_back({"large", true, largeIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{smallMaterial, largeMaterial};
        contact0.setLinearStiffness(120000000.0, 34000000.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(true);
        contact0.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact0.setDampingCoefficients(15874.466906486981, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(true);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-1.0, -1.0, -1.0}, {1.0, 1.0, 1.2003699999999999});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem06-hookehistory-spin: " << error.what() << '\n';
        return 1;
    }
}
