#include "rigidBody.h"

#include "execution/motionIntegration.h"

namespace fundem
{

void rigidBody::updatePositionAndOrientation(Real timeStep) noexcept
{
    if (!isInContainer())
    {
        return;
    }
    pointMass::updatePosition(timeStep);
    execution::integrateQuaternion(orientation_, angularVelocity_, timeStep);
}

void rigidBody::updateVelocityAndAngularVelocity(const Vec3& gravity, Real timeStep) noexcept
{
    if (!isInContainer() || inverseMass() == 0.0)
    {
        return;
    }

    pointMass::updateVelocity(gravity, timeStep);
    execution::integrateAngularVelocity(angularVelocity_, torque_, orientation_, inertiaTensor_, inverseInertiaTensor_, timeStep);
}

} // namespace fundem
