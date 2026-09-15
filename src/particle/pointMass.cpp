#include "pointMass.h"

#include "execution/motionIntegration.h"

namespace fundem
{

void pointMass::updatePosition(Real timeStep) noexcept
{
    if (!isInContainer())
    {
        return;
    }
    execution::integratePosition(position_, velocity_, timeStep);
}

void pointMass::updateVelocity(const Vec3& gravity, Real timeStep) noexcept
{
    if (!isInContainer() || inverseMass_ == 0.0)
    {
        return;
    }
    execution::integrateVelocity(velocity_, force_, inverseMass_, gravity, timeStep);
}

} // namespace fundem
