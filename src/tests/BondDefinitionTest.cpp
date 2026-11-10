// Focused API checks for fracture selection and bounding-sphere reference geometry.
#include "solver/solvers.h"
#include "execution/bondFunctions.h"
#include "execution/cpu/particleFunctions.h"
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "execution/cuda/bondForceKernel.cuh"
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{
using namespace fundem;
using math::Real;
using math::Vec3;

void require(bool value, const char* message)
{
    if (!value)
        throw std::runtime_error(message);
}

bool close(Real a, Real b) { return std::abs(a - b) <= 1.e-12 * std::max(1., std::max(std::abs(a), std::abs(b))); }

bool close(const Vec3& a, const Vec3& b) { return math::norm(a - b) <= 1.e-11 * std::max(1., std::max(math::norm(a), math::norm(b))); }

template <class T, class = void> struct hasMasterTorqueField : std::false_type {};
template <class T> struct hasMasterTorqueField<T, std::void_t<typename T::masterTorqueField>> : std::true_type {};
template <class T, class = void> struct hasSlaveTorqueField : std::false_type {};
template <class T> struct hasSlaveTorqueField<T, std::void_t<typename T::slaveTorqueField>> : std::true_type {};

static_assert(!hasMasterTorqueField<bond>::value && !hasSlaveTorqueField<bond>::value, "Endpoint torques must not occupy persistent Bond or device-layout fields");
static_assert(std::is_same_v<decltype(std::declval<const bond&>().masterTorque()), Vec3> &&
                  std::is_same_v<decltype(std::declval<const bond&>().slaveTorque()), Vec3>,
              "Derived endpoint torques must be returned by value");

struct ExpectedLoads {
    Vec3 force;
    Vec3 torque;
    Vec3 masterTorque;
    Vec3 slaveTorque;
};

// Use the established full EVM response, not the new endpoint-torque getter/helper.
ExpectedLoads expectedLoads(const bond& value, const rigidBody& master, const rigidBody& slave)
{
    const auto rotated = [](const rigidBody& body, const Vec3& local) { return math::rotateUnit(body.orientation(), local); };
    const Vec3 masterOffset = rotated(master, value.masterEndpointLocalPosition());
    const Vec3 slaveOffset = rotated(slave, value.slaveEndpointLocalPosition());
    ExpectedLoads result;
    Vec3 normal;
    Real axialDisplacement, normalEnergy, shearEnergy, bendingEnergy, torsionalEnergy;
    execution::calculateBondForceAndTorque(result.force,
                                           result.torque,
                                           result.masterTorque,
                                           result.slaveTorque,
                                           normal,
                                           axialDisplacement,
                                           normalEnergy,
                                           shearEnergy,
                                           bendingEnergy,
                                           torsionalEnergy,
                                           master.position() + masterOffset,
                                           slave.position() + slaveOffset,
                                           rotated(master, value.masterEndpointLocalNormal()),
                                           rotated(master, value.masterEndpointLocalTangent1()),
                                           rotated(master, value.masterEndpointLocalTangent2()),
                                           rotated(slave, value.slaveEndpointLocalNormal()),
                                           rotated(slave, value.slaveEndpointLocalTangent1()),
                                           rotated(slave, value.slaveEndpointLocalTangent2()),
                                           value.coefficientB1(),
                                           value.coefficientB2(),
                                           value.coefficientB3(),
                                           value.coefficientB4(),
                                           value.bendingStiffness(),
                                           value.equivalentLength());
    execution::applyBondDamageToResponse(result.force,
                                         result.torque,
                                         result.masterTorque,
                                         result.slaveTorque,
                                         normalEnergy,
                                         shearEnergy,
                                         bendingEnergy,
                                         torsionalEnergy,
                                         normal,
                                         axialDisplacement,
                                         value.normalStiffness(),
                                         value.damageFactor());
    result.masterTorque += math::cross(masterOffset, result.force);
    result.slaveTorque += math::cross(slaveOffset, -result.force);
    return result;
}

struct Fixture {
    materialContainer materials;
    particleContainer particles;
    bondContainer bonds;

    Fixture(fractureMode mode, bool strengths)
    {
        materials.host().push_back(material{1.e6, .25, 1., 1000.});
        for (int i = 0; i < 2; ++i)
        {
            particle p;
            p.setRadius(i == 0 ? .1 : .3);
            p.setPosition({Real(i), 0., 0.});
            p.setMaterial(materials, 0);
            particles.host().push_back(p);
        }
        bond b{1.};
        require(b.setConnection(particles, 0, 1), "Connect spheres using the existing bounding-sphere reference geometry");
        require(b.setStiffness(1000., 400., 200., 100.) && b.setCrossSectionArea(1.), "Configure custom stiffness and area");
        require(b.setModeICriticalEnergy(.005) && b.setModeIICriticalEnergy(.005), "Retain positive stored fracture energies");
        require(b.setFractureMode(mode), "Configure explicit fracture mode");
        if (strengths)
            require(b.setTensileStrength(1.) && b.setCohesion(1.), "Configure fracture strengths");
        bonds.host().push_back(b);
    }

    bond& value() { return bonds.host().front(); }
    void pose(Real opening) { particles.host()[1].setPosition({1. + opening, 0., 0.}); }
    void evaluate() { require(value().updateMasterSlaveParticle(particles) && value().calculateForce(), "Evaluate native CPU response"); }
};

void deform(Fixture& fixture, bool compression)
{
    fixture.particles.host()[0].setPosition({.03, -.02, .01});
    fixture.particles.host()[1].setPosition({compression ? .82 : 1.04, .13, -.04});
    fixture.particles.host()[0].setOrientation(math::Quaternion::fromUnitAxisAngle(math::normalizedOrZero(Vec3{1., 2., -1.}), .19));
    fixture.particles.host()[1].setOrientation(math::Quaternion::fromUnitAxisAngle(math::normalizedOrZero(Vec3{-2., 1., 3.}), -.27));
}

void requireLoads(const bond& value, const ExpectedLoads& expected)
{
    require(close(value.force(), expected.force) && close(value.torque(), expected.torque) && close(value.masterTorque(), expected.masterTorque) &&
                close(value.slaveTorque(), expected.slaveTorque),
            "Derived endpoint loads must preserve the full EVM endpoint couples and final-force lever arms");
}

void testDerivedEndpointTorques()
{
    for (bool compression : {false, true})
        for (Real damage : {0., .35, 1.})
        {
            Fixture fixture(fractureMode::unbreakable, false);
            require(fixture.value().setDamageFactor(damage), "Set a prescribed damage level without advancing fracture history");
            deform(fixture, compression);
            const auto expected = expectedLoads(fixture.value(), fixture.particles.host()[0], fixture.particles.host()[1]);
            fixture.evaluate();
            requireLoads(fixture.value(), expected);
            if (damage < 1.)
                require(math::norm(expected.masterTorque) > 1. && math::norm(expected.slaveTorque) > 1., "The regression must exercise nonzero endpoint couples");
        }

    Fixture coincident(fractureMode::unbreakable, false);
    coincident.particles.host()[1].setPosition(coincident.particles.host()[0].position());
    const auto expected = expectedLoads(coincident.value(), coincident.particles.host()[0], coincident.particles.host()[1]);
    coincident.evaluate();
    requireLoads(coincident.value(), expected);
    require(math::isFinite(math::norm(coincident.value().masterTorque())) && math::isFinite(math::norm(coincident.value().slaveTorque())),
            "Coincident endpoints must retain the original finite zero-length direction fallback");

    Fixture cleared(fractureMode::unbreakable, false);
    deform(cleared, false);
    cleared.evaluate();
    require(cleared.value().updateMasterSlaveParticle(cleared.particles), "Refresh snapshots before the next force evaluation");
    require(math::norm(cleared.value().masterTorque()) == 0. && math::norm(cleared.value().slaveTorque()) == 0., "Snapshot refresh must clear derived torques with the other loads");
    require(cleared.value().calculateForce(), "Evaluate the refreshed snapshots");
    require(!cleared.value().calculateForce(), "Reject evaluation without refreshed snapshots");
    require(math::norm(cleared.value().masterTorque()) == 0. && math::norm(cleared.value().slaveTorque()) == 0., "Rejected evaluation must not leave derived endpoint loads");
    const bond unconnected;
    require(math::norm(unconnected.masterTorque()) == 0. && math::norm(unconnected.slaveTorque()) == 0., "Unconnected Bonds have no endpoint loads");
}

void testCPUScatter()
{
    for (int count : {1, cpu::parallelInteractionThreshold})
    {
        Fixture fixture(fractureMode::unbreakable, false);
        require(fixture.value().setDamageFactor(.35), "Configure partial damage for assembly checks");
        deform(fixture, true);
        const auto expected = expectedLoads(fixture.value(), fixture.particles.host()[0], fixture.particles.host()[1]);
        const bond prototype = fixture.value();
        fixture.bonds.host().resize(count, prototype);
        cpu::addBondForceAndTorque(fixture.bonds, fixture.particles);
        require(close(fixture.particles.host()[0].force(), Real(count) * expected.force) && close(fixture.particles.host()[1].force(), -Real(count) * expected.force) &&
                    close(fixture.particles.host()[0].torque(), Real(count) * expected.masterTorque) && close(fixture.particles.host()[1].torque(), Real(count) * expected.slaveTorque),
                "Serial and staged-parallel CPU scatter must assemble the established endpoint loads");
        for (const bond& value : fixture.bonds.host())
            requireLoads(value, expected);
    }
}

void testLevelSetEndpointTorques()
{
    // Only body kinematics are needed here; this check does not perform LS contact detection.
    particleContainer spheres;
    LSParticleContainer levelSets;
    particle sphere;
    sphere.setRadius(.2);
    spheres.host().push_back(sphere);
    for (int index = 0; index < 2; ++index)
    {
        LSParticle value;
        value.setPosition({Real(index + 1), 0., 0.});
        levelSets.host().push_back(value);
    }
    bond mixed{.8}, same{.9};
    require(mixed.setConnection(spheres, 0, levelSets, 0) && same.setConnection(levelSets, 0, 1), "Connect Sphere-LS and LS-LS owner-container pairings");
    for (bond* value : {&mixed, &same})
        require(value->setStiffness(1000., 400., 200., 100.) && value->setCrossSectionArea(0.) && value->setDamageFactor(.25), "Configure mixed-family response");
    bondContainer mixedBonds, sameBonds;
    mixedBonds.host().push_back(mixed);
    sameBonds.host().push_back(same);
    spheres.host()[0].setOrientation(math::Quaternion::fromUnitAxisAngle(Vec3::unitY(), .15));
    levelSets.host()[0].setOrientation(math::Quaternion::fromUnitAxisAngle(Vec3::unitZ(), -.2));
    levelSets.host()[1].setOrientation(math::Quaternion::fromUnitAxisAngle(Vec3::unitX(), .17));
    levelSets.host()[0].setPosition({.9, .12, -.03});
    levelSets.host()[1].setPosition({1.95, -.08, .04});
    const auto mixedExpected = expectedLoads(mixedBonds.host().front(), spheres.host()[0], levelSets.host()[0]);
    const auto sameExpected = expectedLoads(sameBonds.host().front(), levelSets.host()[0], levelSets.host()[1]);
    cpu::addBondForceAndTorque(mixedBonds, spheres, levelSets);
    cpu::addBondForceAndTorque(sameBonds, levelSets);
    requireLoads(mixedBonds.host().front(), mixedExpected);
    requireLoads(sameBonds.host().front(), sameExpected);
    require(close(spheres.host()[0].torque(), mixedExpected.masterTorque) &&
                close(levelSets.host()[0].torque(), mixedExpected.slaveTorque + sameExpected.masterTorque) && close(levelSets.host()[1].torque(), sameExpected.slaveTorque),
            "Derived loads must use the correct owning containers for Sphere-LS and LS-LS Bonds");
}

void testModes()
{
    static_assert(static_cast<int>(fractureMode::unbreakable) == 0 && static_cast<int>(fractureMode::brittle) == 1 && static_cast<int>(fractureMode::damage) == 2,
                  "Renaming a fracture mode must preserve its stored device values");
    bond defaults;
    require(defaults.fractureMode() == fractureMode::damage, "Default fracture selection must preserve existing caller behaviour");
    require(!defaults.setFractureMode(static_cast<fractureMode>(-1)) && defaults.fractureMode() == fractureMode::damage, "Reject an invalid selection without modifying state");

    Fixture unbreakable(fractureMode::unbreakable, false);
    require(unbreakable.value().isValid(), "Positive-area unbreakable bonds do not require fracture strengths");
    unbreakable.pose(.01);
    unbreakable.evaluate();
    require(unbreakable.value().damageFactor() == 0. && !unbreakable.value().getDamageStarted() && close(unbreakable.value().force().x, 10.),
            "Explicit unbreakable selection must skip fracture despite positive area and energies");

    Fixture brittle(fractureMode::brittle, true);
    require(brittle.value().isValid(), "Explicit brittle bond must be valid");
    brittle.pose(.0014);
    brittle.evaluate();
    require(brittle.value().damageFactor() == 1. && brittle.value().getDamageStarted() && math::norm(brittle.value().force()) == 0.,
            "Explicit brittle selection must fail immediately even with positive stored energies");

    Fixture damage(fractureMode::damage, true);
    require(damage.value().isValid(), "Explicit progressive bond must be valid");
    damage.pose(.0014);
    damage.evaluate();
    require(damage.value().getDamageStarted() && damage.value().damageFactor() == 0., "Progressive first crossing must retain the established onset convention");
    damage.pose(.003);
    damage.evaluate();
    require(damage.value().damageFactor() > 0. && damage.value().damageFactor() < 1., "Progressive selection must retain a softening interval");

    Fixture noArea(fractureMode::brittle, false);
    require(noArea.value().setCrossSectionArea(0.) && noArea.value().isValid(), "Legacy zero area must disable fracture in every selection");
    noArea.pose(.01);
    noArea.evaluate();
    require(noArea.value().damageFactor() == 0. && close(noArea.value().force().x, 10.), "Zero area must preserve the elastic response");

    Fixture pureBrittle(fractureMode::brittle, true);
    require(pureBrittle.value().setStiffness(0., 0., 100., 100.) && pureBrittle.value().isValid(), "Explicit brittle selection does not need progressive normal/shear stiffness");
    require(pureBrittle.value().setFractureMode(fractureMode::damage) && !pureBrittle.value().isValid(), "Positive-energy progressive selection still requires both equivalent-component stiffnesses");
}

void testSphereReferenceGeometry()
{
    Fixture f(fractureMode::unbreakable, false);
    require(f.value().equivalentLength() == 1. && close(f.value().point().x, .4), "Separated unequal spheres use the midpoint of their facing surfaces, not their centres");
    require(close(f.value().masterEndpointLocalPosition().x, -.1) && close(f.value().slaveEndpointLocalPosition().x, -.1),
            "Reference anchors retain the supplied equivalent length about the contact point");
    bond reversed{.25};
    require(reversed.setConnection(f.particles, 1, 0) && close(reversed.point().x, .4) && reversed.normal() == Vec3::unitX() && reversed.equivalentLength() == .25,
            "Swapping spheres preserves the bond point and reverses its normal without replacing the equivalent length");
    f.particles.host()[1].setPosition({.3, 0., 0.});
    bond overlapping{.25};
    require(overlapping.setConnection(f.particles, 0, 1) && close(overlapping.point().x, .05), "Overlapping spheres use the midpoint of the penetration segment");
    f.particles.host()[1].setPosition({1., 0., 0.});
    f.particles.host()[1].setRadius(.1);
    bond equal{.25};
    require(equal.setConnection(f.particles, 0, 1) && close(equal.point().x, .5), "Equal spheres use the centre midpoint");
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
void testDeviceEndpointTorques()
{
    for (Real damage : {0., .35, 1.})
    {
        Fixture fixture(fractureMode::unbreakable, false);
        require(fixture.value().setDamageFactor(damage), "Configure a retained device damage level");
        deform(fixture, true);
        require(math::norm(fixture.value().masterParticle().position() - fixture.particles.host()[0].position()) > .01,
                "Device torque readback must be tested with deliberately stale creation-time Bond snapshots");
        const auto expected = expectedLoads(fixture.value(), fixture.particles.host()[0], fixture.particles.host()[1]);
        fixture.particles.setUseDevice(true);
        fixture.bonds.setUseDevice(true);
        fixture.particles.copyHostToDevice();
        fixture.bonds.copyHostToDevice();
        cuda::launchBondForce(fixture.bonds, fixture.particles, nullptr);
        fixture.particles.copyDeviceToHost();
        fixture.bonds.copyDeviceToHost();
        requireLoads(fixture.value(), expected);
        require(close(fixture.particles.host()[0].force(), expected.force) && close(fixture.particles.host()[1].force(), -expected.force) &&
                    close(fixture.particles.host()[0].torque(), expected.masterTorque) && close(fixture.particles.host()[1].torque(), expected.slaveTorque),
                "CUDA scatter and host-derived endpoint torques must agree after current particles and Bonds are downloaded");
        require(math::norm(fixture.value().masterParticle().position() - fixture.particles.host()[0].position()) > .01,
                "Host torque getters must read live owning particles, not silently rely on updated CPU-only snapshots");
    }
}

void testDeviceModes()
{
    for (fractureMode mode : {fractureMode::unbreakable, fractureMode::brittle, fractureMode::damage})
    {
        Fixture cpu(mode, mode != fractureMode::unbreakable), gpu(mode, mode != fractureMode::unbreakable);
        cpu.pose(.0014);
        gpu.pose(.0014);
        cpu.evaluate();
        gpu.particles.setUseDevice(true);
        gpu.bonds.setUseDevice(true);
        gpu.particles.copyHostToDevice();
        gpu.bonds.copyHostToDevice();
        cuda::launchBondForce(gpu.bonds, gpu.particles, nullptr);
        gpu.bonds.copyDeviceToHost();
        require(gpu.value().fractureMode() == mode, "The fracture selector must survive the device-layout round trip");
        require(gpu.value().damageFactor() == cpu.value().damageFactor() && gpu.value().getDamageStarted() == cpu.value().getDamageStarted() &&
                    math::norm(gpu.value().force() - cpu.value().force()) < 1.e-10,
                "All three fracture selections must agree between native CPU and GPU");
    }
}
#endif
} // namespace

int main()
{
    try
    {
        testModes();
        testSphereReferenceGeometry();
        testDerivedEndpointTorques();
        testCPUScatter();
        testLevelSetEndpointTorques();
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testDeviceModes();
        testDeviceEndpointTorques();
#endif
        std::cout << "Bond definition: all focused API checks passed.\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
