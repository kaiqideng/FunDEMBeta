#include "data/myLSObject.h"
#include "solver/SphereDEM.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using namespace fundem;
using math::Real;
using math::Vec3;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

Real selectedEnergy(const particle& value, bool linear, bool angular)
{
    const Vec3 bodySpin = math::inverseRotateUnit(value.orientation(), value.angularVelocity());
    return (linear ? 0.5 * math::normSquared(value.velocity()) / value.inverseMass() : 0.0) +
           (angular ? 0.5 * math::dot(bodySpin, value.inertiaTensor() * bodySpin) : 0.0);
}

void checkState(const particle& value, const particle& before, bool linear, bool angular)
{
    require(math::nearlyEqual(value.velocity(), linear ? Vec3::zero() : before.velocity()), "Reset changed the wrong linear-velocity component.");
    require(math::nearlyEqual(value.angularVelocity(), angular ? Vec3::zero() : before.angularVelocity()), "Reset changed the wrong angular-velocity component.");
    require(value.position() == before.position() && value.orientation() == before.orientation(), "Reset changed particle pose.");
    require(value.materialIndex() == before.materialIndex() && value.inverseMass() == before.inverseMass(), "Reset changed material or mass.");
}

void checkReset(executionMode mode, bool linear, bool angular)
{
    SphereDEM simulation;
    simulation.setExecutionMode(mode);
    simulation.addMaterial(material{1.e9, 0.25, 0.5, 1000.0});
    materialPair pair{0, 0};
    pair.setLinearStiffness(5.e4, 2.5e4);
    pair.setLevelSetStiffness(5.e6, 2.5e6);
    simulation.setMaterialPair(pair);
    levelset::Sphere geometry{0.1};
    geometry.buildSurfaceNode(1);
    geometry.buildLSGrid(0.05, 2, false);
    simulation.addGeometry(geometry);
    for (int index = 0; index < 2; ++index)
    {
        particle sphere;
        sphere.setRadius(0.1);
        sphere.setMaterial(simulation.materials(), 0);
        sphere.setPosition({-1.0, Real(index), 0.0});
        sphere.setVelocity({0.1, 0.2, 0.3});
        sphere.setAngularVelocity({0.2, 0.3, 0.4});
        simulation.addSphere(sphere);
        LSParticle levelSet;
        levelSet.setMaterial(simulation.materials(), 0);
        levelSet.setGeometry(simulation.geometries(), 0);
        levelSet.setPosition({1.0, Real(index), 0.0});
        levelSet.setVelocity({0.3, 0.2, 0.1});
        levelSet.setAngularVelocity({0.4, 0.3, 0.2});
        simulation.addLSParticle(levelSet);
    }
    simulation.setBoundary({-5.0, -5.0, -5.0}, {5.0, 5.0, 5.0});
    simulation.setGravity(Vec3::zero());
    simulation.setTimeStep(1.e-5);
    simulation.solve(1);
    simulation.copyDeviceToHost();
    const particle sphereBefore = simulation.spheres().host()[0];
    const particle otherSphere = simulation.spheres().host()[1];
    const LSParticle levelSetBefore = simulation.LSParticles().host()[0];
    const LSParticle otherLevelSet = simulation.LSParticles().host()[1];
    const Real timeBefore = simulation.time();
    const Real sphereEnergy = simulation.resetSphereVelocities({0}, linear, angular);
    const Real levelSetEnergy = simulation.resetLSParticleVelocities({0}, linear, angular);
    const auto close = [](Real actual, Real expected) { return std::abs(actual - expected) < 1.e-12 * (1.0 + std::abs(expected)); };
    require(close(sphereEnergy, selectedEnergy(sphereBefore, linear, angular)), "Sphere reset returned the wrong removed kinetic energy.");
    require(close(levelSetEnergy, selectedEnergy(levelSetBefore, linear, angular)), "LS reset returned the wrong removed kinetic energy.");
    require(simulation.time() == timeBefore, "Reset advanced solver time.");
    checkState(simulation.spheres().host()[0], sphereBefore, linear, angular);
    checkState(simulation.LSParticles().host()[0], levelSetBefore, linear, angular);
    checkState(simulation.spheres().host()[1], otherSphere, false, false);
    checkState(simulation.LSParticles().host()[1], otherLevelSet, false, false);
    simulation.solve(1);
    simulation.copyDeviceToHost();
    require(math::nearlyEqual(simulation.spheres().host()[0].velocity(), linear ? Vec3::zero() : sphereBefore.velocity()), "Sphere continuation lost the selected reset.");
    require(math::nearlyEqual(simulation.spheres().host()[0].angularVelocity(), angular ? Vec3::zero() : sphereBefore.angularVelocity()), "Sphere continuation changed the unselected angular velocity.");
    require(math::nearlyEqual(simulation.LSParticles().host()[0].velocity(), linear ? Vec3::zero() : levelSetBefore.velocity()), "LS continuation changed the unselected linear velocity.");
    require(math::nearlyEqual(simulation.LSParticles().host()[0].angularVelocity(), angular ? Vec3::zero() : levelSetBefore.angularVelocity()), "LS continuation lost the selected reset.");
}
} // namespace

int main()
{
    try
    {
        for (const bool linear : {false, true})
            for (const bool angular : {false, true})
            {
                checkReset(executionMode::CPU, linear, angular);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
                checkReset(executionMode::GPU, linear, angular);
                checkReset(executionMode::Hybrid, linear, angular);
#endif
            }
        std::cout << "Runtime particle mutation checks passed.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
