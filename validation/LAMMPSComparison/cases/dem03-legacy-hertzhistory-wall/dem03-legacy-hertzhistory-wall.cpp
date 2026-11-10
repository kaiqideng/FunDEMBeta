/** @file dem03-legacy-hertzhistory-wall.cpp
 * @brief Direct native-core LAMMPS comparison: dem03-legacy-hertzhistory-wall.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem03-legacy-hertzhistory-wall-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem03-legacy-hertzhistory-wall");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem03-legacy-hertzhistory-wall", "LAMMPS", "PASS",
            1e-08, 3038, {0.0, 0.0, 0.0},
            {
                {2138, "translational_plus_contact_energy", "sphere", "relative", 0.0034842226023719296, 0.01},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{99999984375.0, 0.25, 1.0, 7000.0});
        const int wallMaterial = simulation.addMaterial(material{99999984375.0, 0.25, 1.0, 7000.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::PlaneWall floorShape{{0.0, 0.0, 1.0}, 0.10000000000000001};
        const std::vector<math::Vec3> floorSurfaceNodes{
            {-0.050000000000000003, 0.050000000000000003, 0.0},
            {-0.050000000000000003, -0.050000000000000003, 0.0},
            {0.050000000000000003, -0.050000000000000003, 0.0},
            {0.050000000000000003, 0.050000000000000003, 0.0},
        };
        const std::vector<int3> floorSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto floorPrepared = comparison::prepareWall(
            floorShape, {-0.055, -0.055, -0.055}, 0.0025000000000000001, {45, 45, 45},
            0.070710678118654766, false, floorSurfaceNodes, floorSurfaceFaces);
        const int floorGeometry = simulation.addGeometry(*floorPrepared);

        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(0.0025000000000000001);
        sphereBody.setMaterial(simulation.materials(), particleMaterial);
        sphereBody.setPosition({0.0, 0.0, 0.0025500000000000002});
        sphereBody.setVelocity({0.0, 0.0, -3.8999999999999999});
        sphereBody.setAngularVelocity({0.0, 0.0, 0.0});
        sphereBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int sphereIndex = simulation.addSphere(sphereBody);
        bodies.push_back({"sphere", true, sphereIndex});

        LSParticle floorBody;
        floorBody.setMaterial(simulation.materials(), wallMaterial);
        floorBody.setGeometry(simulation.geometries(), floorGeometry);
        floorBody.setInfiniteMass();
        floorBody.setPosition({0.0, 0.0, 0.0});
        floorBody.setVelocity({0.0, 0.0, 0.0});
        floorBody.setAngularVelocity({0.0, 0.0, 0.0});
        floorBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int floorIndex = simulation.addLSParticle(floorBody);
        bodies.push_back({"floor", false, floorIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{particleMaterial, particleMaterial};
        contact0.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(false);
        contact0.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        materialPair contact1{particleMaterial, wallMaterial};
        contact1.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact1.setLinearStiffnessEnabled(false);
        contact1.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact1.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact1.setDampingEnabled(false);
        simulation.setMaterialPair(contact1);

        materialPair contact2{wallMaterial, wallMaterial};
        contact2.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact2.setLinearStiffnessEnabled(false);
        contact2.setFrictionCoefficients(0.0, 0.0, 0.0);
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
        std::cerr << "dem03-legacy-hertzhistory-wall: " << error.what() << '\n';
        return 1;
    }
}
