/** @file dem13-twist-sds.cpp
 * @brief Direct native-core LAMMPS comparison: dem13-twist-sds.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem13-twist-sds-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem13-twist-sds");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem13-twist-sds", "LAMMPS", "PASS",
            9.9999999999999995e-07, 50000, {0.0, 0.0, -1.0},
            {
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int sphereMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 2700.0});
        const int wallMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 2700.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::PlaneWall wallShape{{0.0, 0.0, 1.0}, 0.10000000000000001};
        const std::vector<math::Vec3> wallSurfaceNodes{
            {-0.050000000000000003, 0.050000000000000003, 0.0},
            {-0.050000000000000003, -0.050000000000000003, 0.0},
            {0.050000000000000003, -0.050000000000000003, 0.0},
            {0.050000000000000003, 0.050000000000000003, 0.0},
        };
        const std::vector<int3> wallSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto wallPrepared = comparison::prepareWall(
            wallShape, {-0.059999999999999998, -0.059999999999999998, -0.059999999999999998}, 0.0050000000000000001, {25, 25, 25},
            0.070710678118654766, false, wallSurfaceNodes, wallSurfaceFaces);
        const int wallGeometry = simulation.addGeometry(*wallPrepared);

        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(0.0025000000000000001);
        sphereBody.setMaterial(simulation.materials(), sphereMaterial);
        sphereBody.setPosition({0.0, 0.0, 0.0024999998});
        sphereBody.setVelocity({0.0, 0.0, 0.0});
        sphereBody.setAngularVelocity({0.0, 0.0, 10.0});
        sphereBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int sphereIndex = simulation.addSphere(sphereBody);
        bodies.push_back({"sphere", true, sphereIndex});

        LSParticle wallBody;
        wallBody.setMaterial(simulation.materials(), wallMaterial);
        wallBody.setGeometry(simulation.geometries(), wallGeometry);
        wallBody.setInfiniteMass();
        wallBody.setPosition({0.0, 0.0, 0.0});
        wallBody.setVelocity({0.0, 0.0, 0.0});
        wallBody.setAngularVelocity({0.0, 0.0, 0.0});
        wallBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int wallIndex = simulation.addLSParticle(wallBody);
        bodies.push_back({"wall", false, wallIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{sphereMaterial, wallMaterial};
        contact0.setLinearStiffness(1000000.0, 280000.0, 0.0, 40000.0);
        contact0.setLinearStiffnessEnabled(true);
        contact0.setFrictionCoefficients(0.5, 0.0, 0.02);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.040000000000000001);
        contact0.setDampingEnabled(true);
        simulation.setMaterialPair(contact0);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.14142135623730953, -0.14142135623730953, -0.14267135623730953}, {0.14142135623730953, 0.14142135623730953, 0.14267135623730953});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem13-twist-sds: " << error.what() << '\n';
        return 1;
    }
}
