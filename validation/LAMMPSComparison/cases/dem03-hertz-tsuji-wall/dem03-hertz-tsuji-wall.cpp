/** @file dem03-hertz-tsuji-wall.cpp
 * @brief Direct native-core LAMMPS comparison: dem03-hertz-tsuji-wall.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem03-hertz-tsuji-wall-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem03-hertz-tsuji-wall");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem03-hertz-tsuji-wall", "LAMMPS", "PASS",
            1e-08, 7000, {0.0, 0.0, 0.0},
            {
                {7000, "vx", "sphere", "relative", -1.0, 0.0050000000000000001},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{70000000000.0, 0.33000000000000002, 0.49991645217278552, 2700.0});
        const int wallMaterial = simulation.addMaterial(material{70000000000.0, 0.33000000000000002, 0.49991645217278552, 2700.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::PlaneWall wallShape{{-1.0, 0.0, 0.0}, 0.10000000000000001};
        const std::vector<math::Vec3> wallSurfaceNodes{
            {0.0, 0.050000000000000003, -0.050000000000000003},
            {0.0, -0.050000000000000003, -0.050000000000000003},
            {0.0, -0.050000000000000003, 0.050000000000000003},
            {0.0, 0.050000000000000003, 0.050000000000000003},
        };
        const std::vector<int3> wallSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto wallPrepared = comparison::prepareWall(
            wallShape, {-0.055, -0.055, -0.055}, 0.0025000000000000001, {45, 45, 45},
            0.070710678118654766, false, wallSurfaceNodes, wallSurfaceFaces);
        const int wallGeometry = simulation.addGeometry(*wallPrepared);

        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(0.0025000000000000001);
        sphereBody.setMaterial(simulation.materials(), particleMaterial);
        sphereBody.setPosition({-0.0025999999999999999, 0.0, 0.0});
        sphereBody.setVelocity({2.0, 0.0, 0.0});
        sphereBody.setAngularVelocity({0.0, 0.0, 0.0});
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
        materialPair contact0{particleMaterial, particleMaterial};
        contact0.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(false);
        contact0.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        materialPair contact1{particleMaterial, wallMaterial};
        contact1.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact1.setLinearStiffnessEnabled(false);
        contact1.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact1.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact1.setDampingEnabled(false);
        simulation.setMaterialPair(contact1);

        materialPair contact2{wallMaterial, wallMaterial};
        contact2.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact2.setLinearStiffnessEnabled(false);
        contact2.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact2.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact2.setDampingEnabled(false);
        simulation.setMaterialPair(contact2);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.14142135623730953, -0.14142135623730953, -0.14142135623730953}, {0.14142135623730953, 0.14142135623730953, 0.14142135623730953});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem03-hertz-tsuji-wall: " << error.what() << '\n';
        return 1;
    }
}
