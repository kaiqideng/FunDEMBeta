#include "interaction/bond.h"
#include "material/material.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace
{

using namespace fundem;
using Real = math::Real;
using Vec3 = math::Vec3;

class FiberTestSolver : public SphereDEM
{
public:
    inline static constexpr int particleCount = 11;
    inline static constexpr int bondCount = particleCount - 1;
    inline static constexpr int loadedParticleIndex = particleCount - 1;
    inline static constexpr Real loadRampDuration = 1.0;
    inline static constexpr Real verticalLoad = 100.0e3;
    inline static constexpr Real dampingCoefficient = 0.001;

    Real maximumParticleSpeed() const noexcept
    {
        Real result = 0.0;
        for (const particle& value : spheres().host())
        {
            result = std::max(result, math::norm(value.velocity()));
        }
        return result;
    }

protected:
    void addSphereExternalForceAndTorque(particleContainer::host_container_type& particles) override
    {
        const Real loadScale = std::min(time() / loadRampDuration, 1.0);
        particles[loadedParticleIndex].addForce({0.0, 0.0, loadScale * verticalLoad});

        for (particle& value : particles)
        {
            const Vec3 dampingForce = dampingCoefficient * math::norm(value.force()) * math::normalizedOrZero(value.velocity());
            const Vec3 dampingTorque = dampingCoefficient * math::norm(value.torque()) * math::normalizedOrZero(value.angularVelocity());
            value.addForce(-dampingForce);
            value.addTorque(-dampingTorque);
        }
    }
};

} // namespace

int main()
{
    try
    {
        constexpr Real radius = 0.2;
        constexpr Real initialLength = 2.0 * radius;
        constexpr Real density = 7800.0;
        constexpr Real youngsModulus = 200.0e9;
        constexpr Real poissonRatio = 0.3;
        constexpr Real crossSectionArea = math::pi * radius * radius;
        constexpr Real areaMomentOfInertia = 0.25 * crossSectionArea * radius * radius;
        constexpr Real polarMomentOfInertia = 2.0 * areaMomentOfInertia;
        constexpr Real shearModulus = youngsModulus / (2.0 * (1.0 + poissonRatio));
        constexpr Real normalStiffness = youngsModulus * crossSectionArea / initialLength;
        constexpr Real shearStiffness = 12.0 * youngsModulus * areaMomentOfInertia / (initialLength * initialLength * initialLength);
        constexpr Real bendingStiffness = youngsModulus * areaMomentOfInertia / initialLength;
        constexpr Real torsionalStiffness = shearModulus * polarMomentOfInertia / initialLength;
        constexpr Real timeStep = 4.0e-5;
        constexpr Real speedTolerance = 1.0e-9;
        constexpr Real maximumSimulationTime = 400.0;
        constexpr int maximumStepCount = static_cast<int>(maximumSimulationTime / timeStep);

        FiberTestSolver simulation;
        const int materialIndex = simulation.addMaterial(material{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, density});

        for (int particleIndex = 0; particleIndex < FiberTestSolver::particleCount; ++particleIndex)
        {
            particle value;
            value.setPosition({particleIndex * initialLength, 0.0, 0.0});
            value.setRadius(radius);
            value.setMaterial(simulation.materials(), materialIndex);
            if (particleIndex == 0)
            {
                value.setInfiniteMass();
            }
            simulation.addSphere(value);
        }

        for (int bondIndex = 0; bondIndex < FiberTestSolver::bondCount; ++bondIndex)
        {
            bond connection{initialLength};
            if (!connection.setConnection(simulation.spheres(), bondIndex, bondIndex + 1, -Vec3::unitX()) ||
                !connection.setStiffness(normalStiffness, shearStiffness, bendingStiffness, torsionalStiffness))
            {
                throw std::runtime_error("Failed to create the bonded fiber.");
            }
            simulation.addBond(connection);
        }

        simulation.setBoundary({-radius, -radius, -radius}, {FiberTestSolver::bondCount * initialLength + radius, radius, radius});
        simulation.setGravity(Vec3::zero());
        simulation.setTimeStep(timeStep);

        bool converged = false;
        Real maximumSpeed = 0.0;
        for (int stepIndex = 0; stepIndex < maximumStepCount; ++stepIndex)
        {
            simulation.step();
            if (simulation.time() < FiberTestSolver::loadRampDuration)
            {
                continue;
            }
            maximumSpeed = simulation.maximumParticleSpeed();
            if (maximumSpeed < speedTolerance)
            {
                converged = true;
                break;
            }
        }

        if (!converged)
        {
            std::cerr << "Final time=" << simulation.time() << " s, maximum particle speed=" << maximumSpeed << " m/s.\n";
            throw std::runtime_error("The bonded fiber did not bring every particle below the speed tolerance.");
        }
        if (simulation.spheres().host().front().position() != Vec3::zero())
        {
            throw std::runtime_error("The fixed fiber endpoint moved.");
        }
        if (simulation.spheres().host().back().position().z <= 0.0)
        {
            throw std::runtime_error("The loaded fiber endpoint did not deflect upward.");
        }

        constexpr Real fiberLength = FiberTestSolver::bondCount * initialLength;
        constexpr Real analyticalTipDisplacement = FiberTestSolver::verticalLoad * fiberLength * fiberLength * fiberLength / (3.0 * youngsModulus * areaMomentOfInertia);
        const Real simulatedTipDisplacement = simulation.spheres().host().back().position().z;
        const Real absoluteDisplacementError = std::abs(simulatedTipDisplacement - analyticalTipDisplacement);
        const Real relativeDisplacementError = absoluteDisplacementError / analyticalTipDisplacement;

        std::cout << std::setprecision(12) << "Fiber converged at t=" << simulation.time() << " s with maximum particle speed " << maximumSpeed << " m/s.\n"
                  << "Tip displacement: simulated=" << simulatedTipDisplacement << " m, analytical=" << analyticalTipDisplacement << " m, absolute error=" << absoluteDisplacementError
                  << " m, relative error=" << 100.0 * relativeDisplacementError << "%.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
