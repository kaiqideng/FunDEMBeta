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
    constexpr const char* id = "pfc-linear-restitution-md0-sweep";
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
        const int wallMaterial = simulation.addMaterial(material{young, poisson, 1.0, density});
        // Source beta values: 0, .01, .02, .03, .04, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.
        // e_ball = e_contact/(2-e_contact), so harmonic(e_ball, e_wall=1) = e_contact.
        struct BallSettings {
            const char* id;
            Real restitution, x;
        };
        const std::vector<BallSettings> balls{{"ball-01", 1, 0.050000000000000003},
                                              {"ball-02", 0.9399976271057352, 0.15000000000000002},
                                              {"ball-03", 0.88517327192459394, 0.25},
                                              {"ball-04", 0.83489372754441116, 0.35000000000000003},
                                              {"ball-05", 0.78862498034990658, 0.45000000000000001},
                                              {"ball-06", 0.74591352594170146, 0.55000000000000004},
                                              {"ball-07", 0.57387074184313736, 0.65000000000000013},
                                              {"ball-08", 0.35742362021006774, 0.75000000000000011},
                                              {"ball-09", 0.22874735907083726, 0.85000000000000009},
                                              {"ball-10", 0.14536170331729087, 0.95000000000000007},
                                              {"ball-11", 0.088751503041589533, 1.05},
                                              {"ball-12", 0.049747659602316682, 1.1500000000000001},
                                              {"ball-13", 0.023535120638067897, 1.2500000000000002},
                                              {"ball-14", 0.0076402406044936916, 1.3500000000000001},
                                              {"ball-15", 0.00076245881078512495, 1.4500000000000002},
                                              {"ball-16", 0, 1.55}};
        std::vector<int> ballMaterials;
        std::vector<int> materialIndices{wallMaterial};
        for (const auto& ball : balls)
        {
            const int index = simulation.addMaterial(material{young, poisson, ball.restitution, density});
            ballMaterials.push_back(index);
            materialIndices.push_back(index);
        }

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
        LSParticle wall({0.80000000000000004, 0.10000000000000001, 0}, math::Quaternion{1, 0, 0, 0}, {0, 0, 0}, {0, 0, 0});
        wall.setMaterial(simulation.materials(), wallMaterial);
        wall.setGeometry(simulation.geometries(), wallGeometry);
        wall.setInfiniteMass();
        bodies.push_back({"wall", false, simulation.addLSParticle(wall)});
        includeInSearchBounds(wall.position(), wall.velocity(), wall.boundingRadius());

        for (std::size_t i = 0; i < ballMaterials.size(); ++i)
        {
            particle sphere({balls[i].x, 0.050000000000000003, 1.0}, math::Quaternion{1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
            sphere.setRadius(0.025000000000000001);
            sphere.setMaterial(simulation.materials(), ballMaterials[i]);
            bodies.push_back({balls[i].id, true, simulation.addSphere(sphere)});
            includeInSearchBounds(sphere.position(), sphere.velocity(), sphere.radius());
        }

        // All 153 unordered material pairs, including self-pairs, share these explicit coefficients.
        for (std::size_t first = 0; first < materialIndices.size(); ++first)
            for (std::size_t second = first; second < materialIndices.size(); ++second)
            {
                materialPair pair{materialIndices[first], materialIndices[second]};
                pair.setLinearStiffness(50000.0, 0.0, 0.0, 0.0);
                pair.setLinearStiffnessEnabled(true);
                pair.setFrictionCoefficients(0.0, 0.0, 0.0);
                pair.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
                pair.setDampingEnabled(false); // Restitution-derived automatic damping.
                simulation.setMaterialPair(pair);
            }

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
