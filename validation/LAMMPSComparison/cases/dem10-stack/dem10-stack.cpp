/** @file dem10-stack.cpp
 * @brief Direct native-core LAMMPS comparison: dem10-stack.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem10-stack-3d-si.yaml
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
        const auto options = comparison::parseOptions(argc, argv, "dem10-stack");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem10-stack", "LAMMPS", "PASS",
            1.9999999999999999e-06, 10000, {0.0, 0.0, -9.8100000000000005},
            {
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int lowerMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 20000.0});
        const int upperMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 10000.0});
        const int wallMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 20000.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::PlaneWall floorShape{{0.0, 0.0, 1.0}, 0.01};
        const std::vector<math::Vec3> floorSurfaceNodes{
            {-0.0050000000000000001, 0.0050000000000000001, 0.0},
            {-0.0050000000000000001, -0.0050000000000000001, 0.0},
            {0.0050000000000000001, -0.0050000000000000001, 0.0},
            {0.0050000000000000001, 0.0050000000000000001, 0.0},
        };
        const std::vector<int3> floorSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto floorPrepared = comparison::prepareWall(
            floorShape, {-0.0060000000000000001, -0.0060000000000000001, -0.0060000000000000001}, 0.00050000000000000001, {25, 25, 25},
            0.0070710678118654762, false, floorSurfaceNodes, floorSurfaceFaces);
        const int floorGeometry = simulation.addGeometry(*floorPrepared);

        levelset::PlaneWall roofShape{{0.0, 0.0, -1.0}, 0.01};
        const std::vector<math::Vec3> roofSurfaceNodes{
            {-0.0050000000000000001, -0.0050000000000000001, 0.0},
            {-0.0050000000000000001, 0.0050000000000000001, 0.0},
            {0.0050000000000000001, 0.0050000000000000001, 0.0},
            {0.0050000000000000001, -0.0050000000000000001, 0.0},
        };
        const std::vector<int3> roofSurfaceFaces{
            {0, 1, 2},
            {0, 2, 3},
        };
        const auto roofPrepared = comparison::prepareWall(
            roofShape, {-0.0060000000000000001, -0.0060000000000000001, -0.0060000000000000001}, 0.00050000000000000001, {25, 25, 25},
            0.0070710678118654762, false, roofSurfaceNodes, roofSurfaceFaces);
        const int roofGeometry = simulation.addGeometry(*roofPrepared);

        // Initial body states, in the original stable insertion order.
        particle lowerBody;
        lowerBody.setRadius(0.00050000000000000001);
        lowerBody.setMaterial(simulation.materials(), lowerMaterial);
        lowerBody.setPosition({0.0, 0.0, 0.00044999999999999999});
        lowerBody.setVelocity({0.0, 0.0, 0.0});
        lowerBody.setAngularVelocity({0.0, 0.0, 0.0});
        lowerBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int lowerIndex = simulation.addSphere(lowerBody);
        bodies.push_back({"lower", true, lowerIndex});

        particle upperBody;
        upperBody.setRadius(0.00050000000000000001);
        upperBody.setMaterial(simulation.materials(), upperMaterial);
        upperBody.setPosition({0.0, 0.0, 0.0013500000000000001});
        upperBody.setVelocity({0.0, 0.0, 0.0});
        upperBody.setAngularVelocity({0.0, 0.0, 0.0});
        upperBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int upperIndex = simulation.addSphere(upperBody);
        bodies.push_back({"upper", true, upperIndex});

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

        LSParticle roofBody;
        roofBody.setMaterial(simulation.materials(), wallMaterial);
        roofBody.setGeometry(simulation.geometries(), roofGeometry);
        roofBody.setInfiniteMass();
        roofBody.setPosition({0.0, 0.0, 0.0018});
        roofBody.setVelocity({0.0, 0.0, 0.0});
        roofBody.setAngularVelocity({0.0, 0.0, 0.0});
        roofBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int roofIndex = simulation.addLSParticle(roofBody);
        bodies.push_back({"roof", false, roofIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{lowerMaterial, lowerMaterial};
        contact0.setLinearStiffness(1000.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(true);
        contact0.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        materialPair contact1{lowerMaterial, upperMaterial};
        contact1.setLinearStiffness(1000.0, 0.0, 0.0, 0.0);
        contact1.setLinearStiffnessEnabled(true);
        contact1.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact1.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact1.setDampingEnabled(false);
        simulation.setMaterialPair(contact1);

        materialPair contact2{lowerMaterial, wallMaterial};
        contact2.setLinearStiffness(1000.0, 0.0, 0.0, 0.0);
        contact2.setLinearStiffnessEnabled(true);
        contact2.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact2.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact2.setDampingEnabled(false);
        simulation.setMaterialPair(contact2);

        materialPair contact3{upperMaterial, upperMaterial};
        contact3.setLinearStiffness(1000.0, 0.0, 0.0, 0.0);
        contact3.setLinearStiffnessEnabled(true);
        contact3.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact3.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact3.setDampingEnabled(false);
        simulation.setMaterialPair(contact3);

        materialPair contact4{upperMaterial, wallMaterial};
        contact4.setLinearStiffness(1000.0, 0.0, 0.0, 0.0);
        contact4.setLinearStiffnessEnabled(true);
        contact4.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact4.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact4.setDampingEnabled(false);
        simulation.setMaterialPair(contact4);

        materialPair contact5{wallMaterial, wallMaterial};
        contact5.setLinearStiffness(1000.0, 0.0, 0.0, 0.0);
        contact5.setLinearStiffnessEnabled(true);
        contact5.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact5.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact5.setDampingEnabled(false);
        simulation.setMaterialPair(contact5);

        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.014142135623730952, -0.014142135623730952, -0.016104135623730954}, {0.014142135623730952, 0.014142135623730952, 0.01790413562373095});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem10-stack: " << error.what() << '\n';
        return 1;
    }
}
