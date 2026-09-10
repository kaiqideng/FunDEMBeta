/** Verifies a damped linear two-sphere collision against its analytical restitution and duration. */
#include "execution/contactFunctions.h"
#include "material/material.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"
#include "validationCommon.h"

#include <algorithm>
#include <exception>

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;
        using Real = math::Real;
        using Vec3 = math::Vec3;

        constexpr Real radius = 0.05;
        constexpr Real density = 1000.0;
        constexpr Real particleNormalStiffness = 2.0e5;
        constexpr Real pairNormalStiffness = 1.0e5;
        constexpr Real targetRestitution = 0.8;
        constexpr Real impactSpeed = 1.0;

        SphereDEM simulation;
        const int materialIndex = simulation.addMaterial(material{particleNormalStiffness, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, targetRestitution, density});
        for (int particleIndex = 0; particleIndex < 2; ++particleIndex)
        {
            particle sphere;
            sphere.setPosition({particleIndex == 0 ? -radius : radius, 0.0, 0.0});
            sphere.setVelocity({particleIndex == 0 ? 0.5 * impactSpeed : -0.5 * impactSpeed, 0.0, 0.0});
            sphere.setRadius(radius);
            sphere.setMaterial(simulation.materials(), materialIndex);
            simulation.addSphere(sphere);
        }

        const Real particleMass = 1.0 / simulation.spheres().host().front().inverseMass();
        const Real effectiveMass = 0.5 * particleMass;
        const Real dissipation = execution::contactDissipationFactor(targetRestitution);
        const Real undampedFrequency = std::sqrt(pairNormalStiffness / effectiveMass);
        const Real frequencyRatio = std::sqrt(1.0 - dissipation * dissipation);
        const Real dampedFrequency = undampedFrequency * frequencyRatio;
        const Real analyticalContactDuration = math::pi / dampedFrequency;
        const Real releaseAngle = math::pi + std::atan2(-2.0 * dissipation * frequencyRatio, 1.0 - 2.0 * dissipation * dissipation);
        const Real unilateralRestitution = -std::exp(-dissipation * releaseAngle / frequencyRatio) * (std::cos(releaseAngle) - dissipation * std::sin(releaseAngle) / frequencyRatio);
        const Real timeStep = analyticalContactDuration / 1000.0;

        simulation.setBoundary({-0.2, -0.1, -0.1}, {0.2, 0.1, 0.1});
        simulation.setGravity(Vec3::zero());
        simulation.setTimeStep(timeStep);

        auto output = validation::openDataFile(validation::outputDirectory(argc, argv, "normalCollision"), "normalCollision.dat");
        output << "# time overlap normalForce\n";

        bool contactStarted = false;
        bool contactFinished = false;
        Real numericalContactDuration = 0.0;
        for (int stepIndex = 0; stepIndex < 2200; ++stepIndex)
        {
            simulation.step();
            const auto& spheres = simulation.spheres().host();
            const Real overlap = 2.0 * radius - math::distance(spheres[0].position(), spheres[1].position());
            Real normalForce = 0.0;
            for (const contact& value : simulation.sphereInteractions().contacts().host())
            {
                normalForce += value.normalForceMagnitude();
            }
            output << simulation.time() << ' ' << std::max(overlap, Real{0.0}) << ' ' << normalForce << '\n';

            if (overlap > 0.0)
            {
                contactStarted = true;
            }
            const Real separationVelocity = math::dot(spheres[1].velocity() - spheres[0].velocity(), math::normalizedOrZero(spheres[1].position() - spheres[0].position()));
            if (contactStarted && overlap <= 0.0 && separationVelocity > 0.0)
            {
                numericalContactDuration = simulation.time();
                contactFinished = true;
                break;
            }
        }

        validation::require(contactFinished, "The two-sphere collision did not finish.");
        const auto& spheres = simulation.spheres().host();
        const Real reboundSpeed = spheres[1].velocity().x - spheres[0].velocity().x;
        const Real numericalRestitution = reboundSpeed / impactSpeed;
        const Real restitutionError = std::abs(numericalRestitution - unilateralRestitution) / unilateralRestitution;
        const Real durationError = std::abs(numericalContactDuration - analyticalContactDuration) / analyticalContactDuration;
        const Real momentumError = std::abs(spheres[0].velocity().x + spheres[1].velocity().x) / impactSpeed;
        const Real maximumError = std::max({restitutionError, durationError, momentumError});
        constexpr Real tolerance = 5.0e-3;

        std::cout << "numerical restitution=" << numericalRestitution << ", unilateral analytical restitution=" << unilateralRestitution << ", configured restitution=" << targetRestitution
                  << ", numerical duration=" << numericalContactDuration << " s, analytical duration=" << analyticalContactDuration << " s\n";
        validation::require(maximumError <= tolerance, "Normal-collision response exceeded tolerance.");
        validation::printPass("normalCollision", maximumError, tolerance);
    }
    catch (const std::exception& error)
    {
        std::cerr << "normalCollision: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
