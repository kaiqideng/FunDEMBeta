/** @file dem05-slip-mu02.cpp
 * @brief Direct native-core LAMMPS comparison: dem05-slip-mu02.
 * Reference: https://github.com/lammps/lammps/blob/f9720ce2d176d615495d1d6c31349b1f71be1c73/unittest/granular/tests/dem05-slip-mu02-3d-si.yaml
 * Historical FAIL is provenance; current acceptance uses only the declared checks.
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
        const auto options = comparison::parseOptions(argc, argv, "dem05-slip-mu02");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem05-slip-mu02", "LAMMPS", "FAIL",
            1.0000000000000001e-05, 150000, {0.0, 0.0, -1.0},
            {
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int sphereMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 2700.0});
        const int wallMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 2700.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::PlaneWall wallShape{{0.0, 0.0, 1.0}, 3.0};
        const std::vector<math::Vec3> wallSurfaceNodes{
            {-1.5, 1.5, 0.0},
            {-1.5, -1.5, 0.0},
            {1.5, -1.5, 0.0},
            {1.5, 1.5, 0.0},
        };
        const std::vector<int3> wallSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto wallPrepared = comparison::prepareWall(
            wallShape, {-1.6000000000000001, -1.6000000000000001, -1.6000000000000001}, 0.050000000000000003, {65, 65, 65},
            2.1213203435596428, false, wallSurfaceNodes, wallSurfaceFaces);
        const int wallGeometry = simulation.addGeometry(*wallPrepared);

        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(0.0025000000000000001);
        sphereBody.setMaterial(simulation.materials(), sphereMaterial);
        sphereBody.setPosition({0.0, 0.0, 0.0024393000000000001});
        sphereBody.setVelocity({1.0, 0.0, 0.0});
        sphereBody.setAngularVelocity({0.0, 0.0, 0.0});
        sphereBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int sphereIndex = simulation.addSphere(sphereBody);
        bodies.push_back({"sphere", true, sphereIndex});

        LSParticle wallBody;
        wallBody.setMaterial(simulation.materials(), wallMaterial);
        wallBody.setGeometry(simulation.geometries(), wallGeometry);
        wallBody.setInfiniteMass();
        wallBody.setPosition({0.65000000000000002, 0.0, 0.0});
        wallBody.setVelocity({0.0, 0.0, 0.0});
        wallBody.setAngularVelocity({0.0, 0.0, 0.0});
        wallBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int wallIndex = simulation.addLSParticle(wallBody);
        bodies.push_back({"wall", false, wallIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{sphereMaterial, wallMaterial};
        contact0.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(false);
        contact0.setFrictionCoefficients(0.20000000000000001, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-3.6238203435596432, -4.2426406871192857, -5.3676406871192857}, {4.892640687119286, 4.2426406871192857, 5.3676406871192857});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem05-slip-mu02: " << error.what() << '\n';
        return 1;
    }
}
