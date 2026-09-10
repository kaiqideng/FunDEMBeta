/** Verifies tangential history integration, Coulomb return, and history reset. */
#include "execution/contactFunctions.h"
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

        constexpr Real overlap = 1.0e-3;
        constexpr Real timeStep = 1.0e-4;
        constexpr Real normalStiffness = 1.0e5;
        constexpr Real slidingStiffness = 4.0e4;
        constexpr Real frictionCoefficient = 0.3;
        constexpr Real normalForceReference = normalStiffness * overlap;
        constexpr Real frictionLimit = frictionCoefficient * normalForceReference;
        const Vec3 normal = Vec3::unitZ();
        const Vec3 relativeVelocity{1.0, 0.0, 0.0};

        Vec3 springDeformation = Vec3::zero();
        Real normalForce = 0.0;
        Real maximumDiskViolation = 0.0;
        auto output = validation::openDataFile(validation::outputDirectory(argc, argv, "frictionLaw"), "frictionLaw.dat");
        output << "# step deformationX tangentialForceX normalForce frictionLimit\n";

        Vec3 force = Vec3::zero();
        for (int stepIndex = 0; stepIndex < 20; ++stepIndex)
        {
            force =
                execution::calculateLinearContactForce(normalForce, springDeformation, relativeVelocity, normal, overlap, timeStep, normalStiffness, slidingStiffness, frictionCoefficient, 1.0, 1.0);
            const Vec3 tangentialForce = force - math::dot(force, normal) * normal;
            maximumDiskViolation = std::max(maximumDiskViolation, std::max(Real{0.0}, math::norm(tangentialForce) - frictionLimit));
            output << stepIndex + 1 << ' ' << springDeformation.x << ' ' << tangentialForce.x << ' ' << normalForce << ' ' << frictionLimit << '\n';
        }

        const Vec3 finalTangentialForce = force - math::dot(force, normal) * normal;
        const Real normalError = validation::relativeError(normalForce, normalForceReference);
        const Real limitError = std::abs(math::norm(finalTangentialForce) - frictionLimit) / frictionLimit;
        const Real deformationError = std::abs(math::norm(springDeformation) - frictionLimit / slidingStiffness) / (frictionLimit / slidingStiffness);

        execution::calculateLinearContactForce(normalForce, springDeformation, relativeVelocity, normal, 0.0, timeStep, normalStiffness, slidingStiffness, frictionCoefficient, 1.0, 1.0);
        validation::require(springDeformation == Vec3::zero(), "Tangential history was not cleared when contact opened.");

        const Real maximumError = std::max({normalError, limitError, deformationError, maximumDiskViolation});
        constexpr Real tolerance = 1.0e-12;
        validation::require(maximumError <= tolerance, "Coulomb friction invariants exceeded tolerance.");
        validation::printPass("frictionLaw", maximumError, tolerance);
    }
    catch (const std::exception& error)
    {
        std::cerr << "frictionLaw: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
