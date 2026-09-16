#include "TestSupport.h"
#include "execution/cuda/contactForceKernel.cuh"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
using namespace fundem;
using math::Real;
using math::Vec3;

constexpr Real timeStep = 1.0e-4;
constexpr Real overlap = 0.02;
constexpr Real effectiveRadius = 0.1;

void requireNear(Real actual, Real expected, const char* message)
{
    if (!math::isFinite(actual) || std::abs(actual - expected) > 1.0e-11 * (1.0 + std::abs(expected)))
        throw std::runtime_error(message);
}

void requireNear(const Vec3& actual, const Vec3& expected, const char* message)
{
    if (!math::isFinite(actual) || math::norm(actual - expected) > 1.0e-11 * (1.0 + math::norm(expected)))
        throw std::runtime_error(message);
}

particle makeSphere(const materialContainer& materials, int materialIndex, const Vec3& position)
{
    particle value;
    value.setPosition(position);
    value.setRadius(0.1);
    value.setMaterial(materials, materialIndex);
    FUNDEM_TEST_REQUIRE(value.inverseMass() > 0.0);
    return value;
}

template <class MasterStorage, class SlaveStorage>
contact makeContact(const MasterStorage& masters, int masterIndex, const SlaveStorage& slaves, int slaveIndex)
{
    contact value;
    FUNDEM_TEST_REQUIRE(value.setMasterSlaveParticle(masters, masterIndex, slaves, slaveIndex));
    value.setPoint(Vec3::zero());
    value.setNormal(Vec3::unitX());
    value.setOverlap(overlap);
    value.setArea(3.0); // Sphere-LS stiffness must not be multiplied by a nodal area.
    value.setEffectiveMass(1.0);
    value.setEffectiveRadius(effectiveRadius);
    // Stationary particles and restitution 1 remove damping; each spring is
    // beyond its Coulomb limit, so the returned force and torque are analytic.
    value.setSlidingSpringDeformation(Vec3::unitY());
    value.setRollingSpringDeformation(Vec3::unitZ());
    value.setTorsionalSpringDeformation(Vec3::unitX());
    return value;
}

void checkResponse(contactContainer& gpuContacts, contactContainer& cpuContacts, Real normalStiffness,
                   Real rollingStiffness, Real torsionalStiffness, Real rollingFriction, Real torsionalFriction)
{
    FUNDEM_TEST_REQUIRE(cpuContacts.host()[0].calculateForce(timeStep));
    FUNDEM_TEST_REQUIRE_CUDA(cudaDeviceSynchronize());
    gpuContacts.copyDeviceToHost();
    const auto& gpu = gpuContacts.host()[0];
    const auto& cpu = cpuContacts.host()[0];
    const Real normalForce = normalStiffness * overlap;
    // Sliding coefficients are 0.4 and 0.8 in both fixtures: harmonic mean 8/15.
    const Vec3 expectedForce{normalForce, -(8.0 / 15.0) * normalForce, 0.0};
    const Vec3 expectedTorque{-2.0 * effectiveRadius * torsionalFriction * normalForce,
                              effectiveRadius * rollingFriction * normalForce, 0.0};
    requireNear(gpu.normalForceMagnitude(), normalForce, "CUDA contact selected the wrong normal stiffness.");
    requireNear(gpu.force(), expectedForce, "CUDA sliding friction must retain the pair harmonic mean.");
    requireNear(gpu.torque(), expectedTorque, "CUDA rolling/torsional friction selected the wrong material coefficients.");
    requireNear(gpu.rollingSpringDeformation(), {0.0, 0.0, rollingFriction * normalForce / rollingStiffness},
                "CUDA rolling spring did not return to the expected friction limit.");
    requireNear(gpu.torsionalSpringDeformation(), {torsionalFriction * normalForce / torsionalStiffness, 0.0, 0.0},
                "CUDA torsional spring did not return to the expected friction limit.");
    requireNear(cpu.force(), expectedForce, "CPU contact force disagrees with the analytical material rule.");
    requireNear(cpu.torque(), expectedTorque, "CPU contact torque disagrees with the analytical material rule.");
    requireNear(gpu.slidingSpringDeformation(), cpu.slidingSpringDeformation(), "CPU/CUDA sliding histories differ.");
    requireNear(gpu.rollingSpringDeformation(), cpu.rollingSpringDeformation(), "CPU/CUDA rolling histories differ.");
    requireNear(gpu.torsionalSpringDeformation(), cpu.torsionalSpringDeformation(), "CPU/CUDA torsional histories differ.");
    requireNear(gpu.rollingElasticEnergy(), cpu.rollingElasticEnergy(), "CPU/CUDA rolling energies differ.");
    requireNear(gpu.torsionalElasticEnergy(), cpu.torsionalElasticEnergy(), "CPU/CUDA torsional energies differ.");
}

void testSphereLevelSet(bool infiniteMass)
{
    materialContainer materials;
    materials.host().push_back(material{1000.0, 500.0, 200.0, 100.0, 0.4, 0.3, 0.2, 1.0, 1000.0});
    materials.host().push_back(LSMaterial{9000.0, 8000.0, 0.8, 1.0, 1200.0});
    FUNDEM_TEST_REQUIRE(materials.host()[1].rollingFrictionCoefficient() == 0.0);
    FUNDEM_TEST_REQUIRE(materials.host()[1].torsionalFrictionCoefficient() == 0.0);
    particleContainer spheres;
    spheres.host().push_back(makeSphere(materials, 0, Vec3::zero()));
    LSParticleContainer levelSets;
    LSParticle levelSet;
    levelSet.setPosition({-0.19, 0.0, 0.0});
    levelSet.setMaterial(materials, 1);
    // Geometry search is outside this seeded-contact test. Supply the rigid
    // mass state directly; the force launcher consumes no level-set grid.
    static_cast<rigidBody&>(levelSet).setMassAndInertia(2.0, math::Mat3::identity());
    if (infiniteMass)
        levelSet.setInfiniteMass();
    FUNDEM_TEST_REQUIRE((levelSet.inverseMass() == 0.0) == infiniteMass);
    levelSets.host().push_back(levelSet);
    const auto initial = makeContact(spheres, 0, levelSets, 0);
    contactContainer cpuContacts;
    contactContainer gpuContacts;
    cpuContacts.host().push_back(initial);
    gpuContacts.host().push_back(initial);
    materials.copyHostToDevice();
    spheres.copyHostToDevice();
    levelSets.copyHostToDevice();
    gpuContacts.copyHostToDevice();
    cuda::launchContactForce(gpuContacts, spheres, levelSets, materials, timeStep);
    checkResponse(gpuContacts, cpuContacts, 1000.0, 200.0, 100.0, 0.3, 0.2);
    spheres.copyDeviceToHost();
    levelSets.copyDeviceToHost();
    requireNear(spheres.host()[0].force(), gpuContacts.host()[0].force(), "CUDA contact did not apply the master sphere force.");
    requireNear(levelSets.host()[0].force(), -gpuContacts.host()[0].force(), "CUDA contact did not apply the LS reaction force.");
}

void testSphereSphereHarmonicFriction()
{
    materialContainer materials;
    materials.host().push_back(material{1000.0, 500.0, 200.0, 100.0, 0.4, 0.3, 0.2, 1.0, 1000.0});
    materials.host().push_back(material{3000.0, 700.0, 600.0, 300.0, 0.8, 0.9, 0.6, 1.0, 1200.0});
    particleContainer spheres;
    spheres.host().push_back(makeSphere(materials, 0, Vec3::zero()));
    spheres.host().push_back(makeSphere(materials, 1, {-0.19, 0.0, 0.0}));
    const auto initial = makeContact(spheres, 0, spheres, 1);
    contactContainer cpuContacts;
    contactContainer gpuContacts;
    cpuContacts.host().push_back(initial);
    gpuContacts.host().push_back(initial);
    materials.copyHostToDevice();
    spheres.copyHostToDevice();
    gpuContacts.copyHostToDevice();
    cuda::launchContactForce(gpuContacts, spheres, materials, timeStep);
    // Finite sphere pairs retain series stiffness and harmonic friction:
    // kn=750, kr=150, kt=75, mu_r=0.45, mu_t=0.3.
    checkResponse(gpuContacts, cpuContacts, 750.0, 150.0, 75.0, 0.45, 0.3);
}
} // namespace

int main()
{
    try
    {
        testSphereLevelSet(false);
        testSphereLevelSet(true);
        testSphereSphereHarmonicFriction();
        std::cout << "PASS CPU/CUDA mixed-material rolling and torsional friction for finite/fixed LS particles and sphere pairs\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
