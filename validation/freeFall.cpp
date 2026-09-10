/** Verifies velocity-Verlet translation and constant-axis rigid rotation. */
#include "material/material.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"
#include "validationCommon.h"

#include <algorithm>
#include <exception>

namespace
{

using namespace fundem;
using Real = math::Real;
using Vec3 = math::Vec3;

/** Applies one constant host torque after built-in force assembly. */
class ConstantTorqueSolver final : public SphereDEM
{
public:
    void setAppliedTorque(const Vec3& value) noexcept { appliedTorque_ = value; }

protected:
    void addSphereExternalForceAndTorque(particleContainer::host_container_type& particles) override
    {
        if (!particles.empty())
        {
            particles.front().addTorque(appliedTorque_);
        }
    }

private:
    Vec3 appliedTorque_{Vec3::zero()};
};

} // namespace

int main(int argc, char** argv)
{
    try
    {
        constexpr Real radius = 0.1;
        constexpr Real density = 1250.0;
        constexpr Real timeStep = 1.0e-3;
        constexpr int stepCount = 1000;
        constexpr Real duration = stepCount * timeStep;
        const Vec3 initialPosition{0.2, -0.1, 2.0};
        const Vec3 initialVelocity{0.3, -0.2, 1.0};
        const Vec3 gravity{0.0, 0.0, -9.81};

        ConstantTorqueSolver simulation;
        const int materialIndex = simulation.addMaterial(material{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, density});
        particle sphere;
        sphere.setPosition(initialPosition);
        sphere.setVelocity(initialVelocity);
        sphere.setRadius(radius);
        sphere.setMaterial(simulation.materials(), materialIndex);
        simulation.addSphere(sphere);

        const Real axialInertia = simulation.spheres().host().front().inertiaTensor()(2, 2);
        const Real angularAcceleration = 1.0;
        simulation.setAppliedTorque({0.0, 0.0, axialInertia * angularAcceleration});
        simulation.setBoundary({-2.0, -2.0, -5.0}, {2.0, 2.0, 4.0});
        simulation.setGravity(gravity);
        simulation.setTimeStep(timeStep);
        simulation.solve(stepCount);

        const particle& result = simulation.spheres().host().front();
        const Vec3 expectedPosition = initialPosition + duration * initialVelocity + 0.5 * duration * duration * gravity;
        const Vec3 expectedVelocity = initialVelocity + duration * gravity;
        const Vec3 expectedAngularVelocity{0.0, 0.0, angularAcceleration * duration};
        const math::Quaternion expectedOrientation = math::Quaternion::fromUnitAxisAngle(Vec3::unitZ(), 0.5 * angularAcceleration * duration * duration);
        const Vec3 actualAxis = math::rotateUnit(result.orientation(), Vec3::unitX());
        const Vec3 expectedAxis = math::rotateUnit(expectedOrientation, Vec3::unitX());

        const Real positionError = validation::relativeError(result.position(), expectedPosition);
        const Real velocityError = validation::relativeError(result.velocity(), expectedVelocity);
        const Real angularVelocityError = validation::relativeError(result.angularVelocity(), expectedAngularVelocity);
        const Real orientationError = validation::relativeError(actualAxis, expectedAxis);
        const Real maximumError = std::max({positionError, velocityError, angularVelocityError, orientationError});
        constexpr Real tolerance = 1.0e-6;

        auto output = validation::openDataFile(validation::outputDirectory(argc, argv, "freeFall"), "freeFall.dat");
        output << "# time x y z vx vy vz omegaZ expectedX expectedY expectedZ expectedVx expectedVy expectedVz expectedOmegaZ\n";
        output << simulation.time() << ' ' << result.position().x << ' ' << result.position().y << ' ' << result.position().z << ' ' << result.velocity().x << ' ' << result.velocity().y << ' '
               << result.velocity().z << ' ' << result.angularVelocity().z << ' ' << expectedPosition.x << ' ' << expectedPosition.y << ' ' << expectedPosition.z << ' ' << expectedVelocity.x << ' '
               << expectedVelocity.y << ' ' << expectedVelocity.z << ' ' << expectedAngularVelocity.z << '\n';

        validation::require(maximumError <= tolerance, "Free-fall or constant-torque integration exceeded tolerance.");
        validation::printPass("freeFall", maximumError, tolerance);
    }
    catch (const std::exception& error)
    {
        std::cerr << "freeFall: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
