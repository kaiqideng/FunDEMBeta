/**
 * @brief PFC official linear-restitution Md=0 comparison using direct native-core construction.
 * @see https://docs.itascacg.com/pfc700/pfc/docproject/source/manual/examples/verification/linear_contact_model/cmlinear_restitution.html
 */
#include "comparisonCommon.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

int main(int argc, char** argv)
{
    using namespace fundem;
    using namespace fundem::comparison;
    using Real = math::Real;
    using Vec3 = math::Vec3;
    constexpr const char* id = "pfc-linear-restitution-md0-beta06";
    try
    {
        const RunOptions options = parseOptions(argc, argv, id);
        if (options.help)
        {
            printUsage(id);
            return 0;
        }
        SphereDEM simulation(options.mode);
        const CaseSettings settings{id, "PFC", "PASS", 1.0000000000000001e-05, 100000, {0.0, 0.0, -10.0}, {}};
        // Historical PASS is provenance; no fresh acceptance gate is configured.
        const Real young = 10000.0, poisson = 0.33000000000000002, density = 2650.0;
        const Real restitution = 0.094780224842154856;
        // Identical material coefficients retain the harmonic contact restitution.
        const int sphereMaterial = simulation.addMaterial(material{young, poisson, restitution, density});
        const int wallMaterial = simulation.addMaterial(material{young, poisson, restitution, density});

        // Retain the accepted analytic plane and fixed grid phase.
        const levelset::PlaneWall plane({0.0, 0.0, 1.0}, 2.0);
        const Real wallRadius = 1.4142135623730951;
        const std::vector<Vec3> nodes{{-1.0, 1.0, 0.0}, {-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, {1.0, 1.0, 0.0}};
        const std::vector<int3> faces{{0, 1, 2}, {0, 2, 3}};
        const auto geometry = prepareWall(plane, {-1.05, -1.05, -1.05}, 0.025000000000000001, {85, 85, 85}, wallRadius, false, nodes, faces);
        const int wallGeometry = simulation.addGeometry(*geometry);

        // Search bounds include every body radius and ballistic travel over the full duration.
        Vec3 lower{std::numeric_limits<Real>::max()}, upper{-std::numeric_limits<Real>::max()};
        Real maximumRadius = 0.0;
        const Real duration = settings.timeStep * settings.steps;
        const auto includeInSearchBounds = [&](const Vec3& position, const Vec3& velocity, Real radius)
        {
            maximumRadius = std::max(maximumRadius, radius);
            for (int axis = 0; axis < 3; ++axis)
            {
                const Real travel = std::abs(velocity[axis]) * duration + 0.5 * std::abs(settings.gravity[axis]) * duration * duration;
                lower[axis] = std::min(lower[axis], position[axis] - radius - travel);
                upper[axis] = std::max(upper[axis], position[axis] + radius + travel);
            }
        };
        std::vector<Body> bodies;
        particle sphere({1.1499999999999999, 0.050000000000000003, 1}, math::Quaternion{1, 0, 0, 0}, {0, 0, 0}, {0, 0, 0});
        sphere.setRadius(0.025000000000000001);
        sphere.setMaterial(simulation.materials(), sphereMaterial);
        bodies.push_back({"sphere", true, simulation.addSphere(sphere)});
        includeInSearchBounds(sphere.position(), sphere.velocity(), sphere.radius());

        LSParticle wall({1.1499999999999999, 0.050000000000000003, 0}, math::Quaternion{1, 0, 0, 0}, {0, 0, 0}, {0, 0, 0});
        wall.setMaterial(simulation.materials(), wallMaterial);
        wall.setGeometry(simulation.geometries(), wallGeometry);
        wall.setInfiniteMass();
        bodies.push_back({"wall", false, simulation.addLSParticle(wall)});
        includeInSearchBounds(wall.position(), wall.velocity(), wall.boundingRadius());

        materialPair pair{sphereMaterial, wallMaterial};
        pair.setLinearStiffness(50000.0, 0.0, 0.0, 0.0);
        pair.setLinearStiffnessEnabled(true);
        pair.setFrictionCoefficients(0.0, 0.0, 0.0);
        pair.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        pair.setDampingEnabled(false); // Restitution-derived automatic damping.
        simulation.setMaterialPair(pair);

        simulation.setBoundary(lower - Vec3{maximumRadius}, upper + Vec3{maximumRadius});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << id << ": ERROR: " << error.what() << '\n';
        return 1;
    }
}
