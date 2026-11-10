/** @file dem11-regionwall.cpp
 * @brief Direct native-core LAMMPS comparison: dem11-regionwall.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem11-regionwall-3d-si.yaml
 * Historical PASS is provenance; current acceptance uses only the declared checks.
 */
#include "comparisonCommon.h"
#include "dem11-regionwall-cylinder-geometry.h"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;
        const auto options = comparison::parseOptions(argc, argv, "dem11-regionwall");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        SphereDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem11-regionwall", "LAMMPS", "PASS",
            1e-08, 20500, {0.0, 0.0, 0.0},
            {
                {20500, "vx", "sphere", "relative", -1.0, 0.01},
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int sphereMaterial = simulation.addMaterial(material{106904480061.14214, 0.25, 0.5, 7000.0});
        const int wallMaterial = simulation.addMaterial(material{106904480061.14214, 0.25, 0.5, 7000.0});

        // Preserve the prepared local grid, surface and sampling phase.
        levelset::CylinderWall cylinderShape{{0.0, 0.0, -0.050000000000000003}, {0.0, 0.0, 0.050000000000000003}, 0.01};
        const auto cylinderPrepared = comparison::prepareWall(
            cylinderShape, {-0.012500000000000001, -0.012625000000000001, -0.052499999999999998}, 0.00025000000000000001, {101, 102, 421},
            0.050990195135927854, true, comparison::dem11_regionwall_cylinder::surfaceNodes, comparison::dem11_regionwall_cylinder::surfaceFaces);
        const int cylinderGeometry = simulation.addGeometry(*cylinderPrepared);

        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(0.0025000000000000001);
        sphereBody.setMaterial(simulation.materials(), sphereMaterial);
        sphereBody.setPosition({0.0071999999999999998, 0.0, 0.0});
        sphereBody.setVelocity({2.0, 0.0, 0.0});
        sphereBody.setAngularVelocity({0.0, 0.0, 0.0});
        sphereBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int sphereIndex = simulation.addSphere(sphereBody);
        bodies.push_back({"sphere", true, sphereIndex});

        LSParticle cylinderBody;
        cylinderBody.setMaterial(simulation.materials(), wallMaterial);
        cylinderBody.setGeometry(simulation.geometries(), cylinderGeometry);
        cylinderBody.setInfiniteMass();
        cylinderBody.setPosition({0.0, 0.0, 0.0});
        cylinderBody.setVelocity({0.0, 0.0, 0.0});
        cylinderBody.setAngularVelocity({0.0, 0.0, 0.0});
        cylinderBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int cylinderIndex = simulation.addLSParticle(cylinderBody);
        bodies.push_back({"cylinder", false, cylinderIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        materialPair contact0{sphereMaterial, sphereMaterial};
        contact0.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
        contact0.setLinearStiffnessEnabled(false);
        contact0.setFrictionCoefficients(0.0, 0.0, 0.0);
        contact0.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        contact0.setDampingEnabled(false);
        simulation.setMaterialPair(contact0);

        materialPair contact1{sphereMaterial, wallMaterial};
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
        simulation.setBoundary({-0.10198039027185571, -0.10198039027185571, -0.10198039027185571}, {0.10198039027185571, 0.10198039027185571, 0.10198039027185571});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem11-regionwall: " << error.what() << '\n';
        return 1;
    }
}
