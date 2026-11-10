/** @file dem05-sliptime.cpp
 * @brief Direct native-core LAMMPS comparison: dem05-sliptime.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem05-sliptime-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem05-sliptime");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem05-sliptime", "LAMMPS", "PASS",
            1.0000000000000001e-05, 100000, {0.0, 0.0, -1.0},
            {
                {40000, "vx", "sphere", "relative", 1.2399999999998963, 0.01},
                {40000, "wy", "sphere", "relative", 160.00000000010371, 0.01},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int particleMaterial = simulation.addMaterial(material{10000.0, 0.25, 1.0, 2700.0});
        const int wallMaterial = simulation.addMaterial(material{10000.0, 0.25, 1.0, 2700.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::PlaneWall floorShape{{0.0, 0.0, 1.0}, 4.0};
        const std::vector<math::Vec3> floorSurfaceNodes{
            {-2.0, 2.0, 0.0},
            {-2.0, -2.0, 0.0},
            {2.0, -2.0, 0.0},
            {2.0, 2.0, 0.0},
        };
        const std::vector<int3> floorSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto floorPrepared = comparison::prepareWall(
            floorShape, {-2.1000000000000001, -2.1000000000000001, -2.1000000000000001}, 0.050000000000000003, {85, 85, 85},
            2.8284271247461903, false, floorSurfaceNodes, floorSurfaceFaces);
        const int floorGeometry = simulation.addGeometry(*floorPrepared);

        levelset::PlaneWall roofShape{{0.0, 0.0, -1.0}, 4.0};
        const std::vector<math::Vec3> roofSurfaceNodes{
            {-2.0, -2.0, 0.0},
            {-2.0, 2.0, 0.0},
            {2.0, 2.0, 0.0},
            {2.0, -2.0, 0.0},
        };
        const std::vector<int3> roofSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto roofPrepared = comparison::prepareWall(
            roofShape, {-2.1000000000000001, -2.1000000000000001, -2.1000000000000001}, 0.050000000000000003, {85, 85, 85},
            2.8284271247461903, false, roofSurfaceNodes, roofSurfaceFaces);
        const int roofGeometry = simulation.addGeometry(*roofPrepared);

        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(0.0025000000000000001);
        sphereBody.setMaterial(simulation.materials(), particleMaterial);
        sphereBody.setPosition({0.0, 0.0, 0.002499982});
        sphereBody.setVelocity({1.3999999999999999, 0.0, 0.0});
        sphereBody.setAngularVelocity({0.0, 0.0, 0.0});
        sphereBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int sphereIndex = simulation.addSphere(sphereBody);
        bodies.push_back({"sphere", true, sphereIndex});

        LSParticle floorBody;
        floorBody.setMaterial(simulation.materials(), wallMaterial);
        floorBody.setGeometry(simulation.geometries(), floorGeometry);
        floorBody.setInfiniteMass();
        floorBody.setPosition({0.69999999999999996, 0.0, 0.0});
        floorBody.setVelocity({0.0, 0.0, 0.0});
        floorBody.setAngularVelocity({0.0, 0.0, 0.0});
        floorBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int floorIndex = simulation.addLSParticle(floorBody);
        bodies.push_back({"floor", false, floorIndex});

        LSParticle roofBody;
        roofBody.setMaterial(simulation.materials(), wallMaterial);
        roofBody.setGeometry(simulation.geometries(), roofGeometry);
        roofBody.setInfiniteMass();
        roofBody.setPosition({0.69999999999999996, 0.0, 1.0});
        roofBody.setVelocity({0.0, 0.0, 0.0});
        roofBody.setAngularVelocity({0.0, 0.0, 0.0});
        roofBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int roofIndex = simulation.addLSParticle(roofBody);
        bodies.push_back({"roof", false, roofIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{particleMaterial, particleMaterial};
        contact0.setLinearStiffness(10000.0, 2800.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(true);
        contact0.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        materialPair contact1{particleMaterial, wallMaterial};
        contact1.setLinearStiffness(10000.0, 2800.0, 0.0, 0.0);
        contact1.setLinearStiffnessEnabled(true);
        contact1.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact1.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact1.setDampingEnabled(false);
        simulation.setMaterialPair(contact1);

        materialPair contact2{wallMaterial, wallMaterial};
        contact2.setLinearStiffness(10000.0, 2800.0, 0.0, 0.0);
        contact2.setLinearStiffnessEnabled(true);
        contact2.setFrictionCoefficients(0.40000000000000002, 0.0, 0.0);
        contact2.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact2.setDampingEnabled(false);
        simulation.setMaterialPair(contact2);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-4.9568542494923804, -5.6568542494923806, -6.1568542494923806}, {6.3568542494923808, 5.6568542494923806, 7.1568542494923797});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem05-sliptime: " << error.what() << '\n';
        return 1;
    }
}
