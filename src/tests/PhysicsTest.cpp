#include "execution/bondFunctions.h"
#include "execution/sphFunctions.h"
#include "interaction/contact.h"
#include "interaction/contactSearch.h"
#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "material/material.h"
#include "solver/LSDEM.h"
#include "solver/SPHDEM.h"
#include "solver/SphereDEM.h"
#include "execution/cpu/particleFunctions.h"
#include "execution/contactFunctions.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{

using namespace fundem;

template <class T, class = void> struct hasRollingFrictionSetter : std::false_type {};
template <class T> struct hasRollingFrictionSetter<T, std::void_t<decltype(std::declval<T&>().setRollingFrictionCoefficient(0.0))>> : std::true_type {};
template <class T, class = void> struct hasTorsionalFrictionSetter : std::false_type {};
template <class T> struct hasTorsionalFrictionSetter<T, std::void_t<decltype(std::declval<T&>().setTorsionalFrictionCoefficient(0.0))>> : std::true_type {};

static_assert(hasRollingFrictionSetter<material>::value && hasTorsionalFrictionSetter<material>::value, "Sphere materials must retain rotational friction settings.");
static_assert(!hasRollingFrictionSetter<LSMaterial>::value && !hasTorsionalFrictionSetter<LSMaterial>::value, "LSMaterial must not expose rotational friction settings.");

bool nearlyEqual(math::Real first, math::Real second, math::Real tolerance = 1.0e-12) noexcept { return std::abs(first - second) <= tolerance * (1.0 + std::max(std::abs(first), std::abs(second))); }

void require(bool condition, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

material testMaterial() { return {1.0e5, 5.0e4, 0.0, 0.0, 0.4, 0.0, 0.0, 0.5, 1000.0}; }

void testSphereDEMExecutionModes()
{
    SphereDEM simulation;
    require(simulation.mode() == executionMode::CPU && !simulation.spheres().usesDevice() && !simulation.LSParticles().usesDevice(), "SphereDEM must use the pure CPU mode by default.");

    bool unsupportedHybridRejected = false;
    try
    {
        LSDEM LSOnly;
        LSOnly.setExecutionMode(executionMode::Hybrid);
    }
    catch (const std::logic_error&)
    {
        unsupportedHybridRejected = true;
    }
    require(unsupportedHybridRejected, "LSDEM must reject Hybrid mode because it contains only one particle type.");

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    simulation.setExecutionMode(executionMode::Hybrid);
    require(simulation.mode() == executionMode::Hybrid && simulation.spheres().usesDevice() && !simulation.LSParticles().usesDevice(),
            "The Sphere GPU and LSParticle CPU mode must set only the sphere backend to GPU.");

    simulation.setExecutionMode(executionMode::GPU);
    require(simulation.mode() == executionMode::GPU && simulation.spheres().usesDevice() && simulation.LSParticles().usesDevice(), "The pure GPU mode must set both particle containers to GPU.");
#endif

    simulation.setExecutionMode(executionMode::CPU);
    require(simulation.mode() == executionMode::CPU && !simulation.spheres().usesDevice() && !simulation.LSParticles().usesDevice(), "Selecting the pure CPU mode must clear both GPU backend tags.");
}

void testSPHDEMExecutionModes()
{
    SPHDEM simulation;
    require(simulation.mode() == executionMode::CPU && !simulation.SPHParticles().usesDevice() && !simulation.LSParticles().usesDevice(), "SPHDEM must use the pure CPU mode by default.");

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    simulation.setExecutionMode(executionMode::Hybrid);
    require(simulation.SPHParticles().usesDevice() && !simulation.LSParticles().usesDevice(), "Hybrid SPHDEM must place only SPHParticles on the GPU.");

    simulation.setExecutionMode(executionMode::GPU);
    require(simulation.SPHParticles().usesDevice() && simulation.LSParticles().usesDevice(), "GPU SPHDEM must place both particle containers on the GPU.");
#endif

    simulation.setExecutionMode(executionMode::CPU);
    require(!simulation.SPHParticles().usesDevice() && !simulation.LSParticles().usesDevice(), "CPU SPHDEM must keep both particle containers on the host.");
}

void testContainerOwnership()
{
    SphereDEM firstSolver;
    SphereDEM secondSolver;
    firstSolver.addMaterial(testMaterial());
    secondSolver.addMaterial(testMaterial());

    particle sphere;
    sphere.setPosition({0.0, 0.0, 0.0});
    sphere.setRadius(0.1);
    sphere.setMaterial(firstSolver.materials(), 0);
    require(firstSolver.addSphere(sphere) == 0, "A sphere referencing the solver material must be accepted.");

    bool rejected = false;
    try
    {
        secondSolver.addSphere(sphere);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "A sphere referencing another solver material container must be rejected.");

    SphereDEM materialTypeSolver;
    materialTypeSolver.addMaterial(LSMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0});
    particle sphereWithLevelSetMaterial;
    sphereWithLevelSetMaterial.setPosition({0.0, 0.0, 0.0});
    sphereWithLevelSetMaterial.setRadius(0.1);
    sphereWithLevelSetMaterial.setMaterial(materialTypeSolver.materials(), 0);
    rejected = false;
    try
    {
        materialTypeSolver.addSphere(sphereWithLevelSetMaterial);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "A sphere must reject LSMaterial so the material marker cannot misclassify a sphere-sphere contact as nodal contact.");
}

void configureContinuationTest(SphereDEM& simulation)
{
    simulation.addMaterial(testMaterial());
    for (int index = 0; index < 2; ++index)
    {
        particle value;
        value.setPosition({0.15 * index, 0.0, 0.0});
        value.setRadius(0.1);
        value.setMaterial(simulation.materials(), 0);
        simulation.addSphere(value);
    }
    simulation.setBoundary({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5});
    simulation.setGravity(math::Vec3::zero());
    simulation.setTimeStep(1.0e-5);
}

void testSplitSolveContinuation()
{
    SphereDEM uninterrupted;
    SphereDEM split;
    configureContinuationTest(uninterrupted);
    configureContinuationTest(split);
    uninterrupted.solve(30);
    split.solve(11);
    split.solve(19);

    for (int particleIndex = 0; particleIndex < 2; ++particleIndex)
    {
        const particle& expected = uninterrupted.spheres().host()[particleIndex];
        const particle& actual = split.spheres().host()[particleIndex];
        require(math::nearlyEqual(actual.position(), expected.position(), 1.0e-13, 1.0e-13), "Splitting solve() must not change a particle position.");
        require(math::nearlyEqual(actual.velocity(), expected.velocity(), 1.0e-13, 1.0e-13), "Splitting solve() must not change a particle velocity.");
    }
}

void testAppendAfterSolve()
{
    SphereDEM simulation;
    configureContinuationTest(simulation);
    simulation.solve(1);

    const int addedMaterialIndex = simulation.addMaterial(testMaterial());
    const int addedLSMaterialIndex = simulation.addMaterial(LSMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0});

    particle addedSphere;
    addedSphere.setPosition({0.4, 0.0, 0.0});
    addedSphere.setRadius(0.05);
    addedSphere.setMaterial(simulation.materials(), addedMaterialIndex);
    require(simulation.addSphere(addedSphere) == 2, "A sphere must be appendable after solve() completes.");

    levelset::PlaneWall wall{{0.0, 0.0, 1.0}, 1.0};
    wall.buildLSGrid(0.1, 2);
    const int geometryIndex = simulation.addGeometry(wall, true);
    LSParticle fixedWall;
    fixedWall.setPosition({0.0, 0.0, -0.4});
    fixedWall.setMaterial(simulation.materials(), addedLSMaterialIndex);
    fixedWall.setGeometry(simulation.geometries(), geometryIndex);
    require(simulation.addLSParticle(fixedWall) == 0, "An LSParticle must be appendable after solve() completes.");

    bond addedBond;
    require(addedBond.setEquivalentLength(0.15), "An appended bond requires a positive equivalent length.");
    require(addedBond.setConnection(simulation.spheres(), 0, 1, math::Vec3::unitX()), "An appended bond must reference existing solver particles.");
    require(addedBond.setStiffness(1.0e5, 5.0e4, 1.0e3, 1.0e3), "An appended bond requires valid stiffness.");
    require(simulation.addBond(addedBond) == 0, "A bond must be appendable after solve() completes.");

    simulation.solve(1);
    require(simulation.materials().hostSize() == 3, "Reinitialization after an append must retain every material.");
    require(simulation.geometries().hostDescriptors().size() == 1, "Reinitialization after an append must retain every geometry.");
    require(simulation.spheres().hostSize() == 3, "Reinitialization after an append must retain every sphere.");
    require(simulation.LSParticles().hostSize() == 1, "Reinitialization after an append must retain every LSParticle.");
    require(simulation.sphereInteractions().bonds().hostSize() == 1, "Reinitialization after an append must retain every bond.");
    require(math::isFinite(simulation.spheres().host().back().position()), "An appended sphere must remain valid after the next solve().");
}

void testLevelSetContainerOwnership()
{
    LSDEM firstSolver;
    LSDEM secondSolver;
    firstSolver.addMaterial(LSMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0});
    secondSolver.addMaterial(LSMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0});
    levelset::PlaneWall wall{{0.0, 0.0, 1.0}, 1.0};
    wall.buildLSGrid(0.1, 2);
    firstSolver.addGeometry(wall, true);
    secondSolver.addGeometry(wall, true);

    LSParticle particle;
    particle.setMaterial(firstSolver.materials(), 0);
    particle.setGeometry(firstSolver.geometries(), 0);
    require(firstSolver.addLSParticle(particle) == 0, "An LSParticle referencing the solver containers must be accepted.");

    bool rejected = false;
    try
    {
        secondSolver.addLSParticle(particle);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "An LSParticle referencing another solver material or geometry container must be rejected.");
}

void testLSParticleMotionReset(executionMode mode)
{
    LSDEM simulation;
    simulation.setExecutionMode(mode);
    simulation.addMaterial(LSMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0});
    levelset::PlaneWall wall{{0.0, 0.0, 1.0}, 1.0};
    wall.buildLSGrid(0.1, 2);
    simulation.addGeometry(wall, true);

    LSParticle particle;
    particle.setMaterial(simulation.materials(), 0);
    particle.setGeometry(simulation.geometries(), 0);
    const int particleIndex = simulation.addLSParticle(particle);

    simulation.setBoundary({-2.0, -2.0, -2.0}, {2.0, 2.0, 2.0});
    simulation.setGravity(math::Vec3::zero());
    simulation.setTimeStep(0.1);
    simulation.setLSParticlePosition(particleIndex, {0.1, 0.2, 0.3});
    simulation.setLSParticleVelocity(particleIndex, {0.5, 0.0, 0.0});
    simulation.setLSParticleAngularVelocity(particleIndex, {0.0, 0.0, 1.0});
    simulation.solve(1);
    require(math::nearlyEqual(simulation.LSParticles().host()[particleIndex].position(), math::Vec3{0.15, 0.2, 0.3}), "The pre-solve LSParticle motion state must drive integration.");

    simulation.setLSParticlePosition(particleIndex, {-0.4, 0.0, 0.0});
    simulation.setLSParticleVelocity(particleIndex, {0.0, 0.25, 0.0});
    simulation.setLSParticleAngularVelocity(particleIndex, {0.0, 2.0, 0.0});
    simulation.solve(1);
    const LSParticle& updatedParticle = simulation.LSParticles().host()[particleIndex];
    require(math::nearlyEqual(updatedParticle.position(), math::Vec3{-0.4, 0.025, 0.0}), "The post-solve LSParticle position and velocity reset must survive reinitialization.");
    require(math::nearlyEqual(updatedParticle.velocity(), math::Vec3{0.0, 0.25, 0.0}), "The post-solve LSParticle velocity reset must survive reinitialization.");
    require(math::nearlyEqual(updatedParticle.angularVelocity(), math::Vec3{0.0, 2.0, 0.0}), "The post-solve LSParticle angular-velocity reset must survive reinitialization.");

    bool rejected = false;
    try
    {
        simulation.setLSParticlePosition(1, math::Vec3::zero());
    }
    catch (const std::out_of_range&)
    {
        rejected = true;
    }
    require(rejected, "An LSParticle motion setter must reject an invalid particle index.");
}

void testContactCalculationGate()
{
    materialContainer materials;
    materials.host().push_back(testMaterial());
    particleContainer particles;
    for (int index = 0; index < 2; ++index)
    {
        particle value;
        value.setPosition({0.15 * index, 0.0, 0.0});
        value.setRadius(0.1);
        value.setMaterial(materials, 0);
        particles.host().push_back(value);
    }

    contact value;
    require(value.setMasterSlaveParticle(particles, 0, 1), "Valid contact particle copies must be accepted.");
    value.setPoint({0.075, 0.0, 0.0});
    value.setNormal({-1.0, 0.0, 0.0});
    value.setOverlap(0.05);
    value.setArea(2.0);
    value.setEffectiveMass(1.0);
    value.setEffectiveRadius(0.0);
    require(!value.calculateForce(1.0e-4), "A contact outside a container must not calculate force.");

    require(value.setMasterSlaveParticle(particles, 0, 1), "Contact particle copies must be refreshable.");
    contactContainer contacts;
    contacts.host().push_back(value);
    require(contacts.host()[0].calculateForce(1.0e-4), "A bound contact with refreshed particle copies must calculate force.");
    require(nearlyEqual(contacts.host()[0].normalForceMagnitude(), 2500.0), "A standard material must select the spherical stiffness model independently of effective radius.");
    require(!contacts.host()[0].calculateForce(1.0e-4), "A second contact calculation without refreshing particle copies must be rejected.");
}

void testParallelContactForceAssembly()
{
    constexpr int contactCount = 512;
    materialContainer materials;
    materials.host().push_back(testMaterial());
    particleContainer particles;
    for (int index = 0; index < 2; ++index)
    {
        particle value;
        value.setPosition({0.15 * index, 0.0, 0.0});
        value.setRadius(0.1);
        value.setMaterial(materials, 0);
        particles.host().push_back(value);
    }

    contactContainer contacts;
    contacts.host().reserve(contactCount);
    for (int contactIndex = 0; contactIndex < contactCount; ++contactIndex)
    {
        contact value;
        require(value.setMasterSlaveParticle(particles, 0, 1), "Parallel contact assembly requires valid particle copies.");
        value.setPoint({0.075, 0.0, 0.0});
        value.setNormal({-1.0, 0.0, 0.0});
        value.setOverlap(0.05);
        value.setArea(2.0);
        value.setEffectiveMass(1.0);
        value.setEffectiveRadius(0.0);
        contacts.host().push_back(value);
    }

    cpu::addContactForceAndTorque(contacts, particles, particles, 1.0e-4);
    const math::Vec3 expectedMasterForce = math::Real(contactCount) * contacts.host().front().force();
    require(math::norm(particles.host()[0].force() - expectedMasterForce) <= 1.0e-10 * math::norm(expectedMasterForce), "Parallel contact force gathering must retain every master contribution.");
    require(math::norm(particles.host()[1].force() + expectedMasterForce) <= 1.0e-10 * math::norm(expectedMasterForce), "Parallel contact force gathering must retain every slave contribution.");
}

void testParallelGridContactSearch()
{
    materialContainer materials;
    materials.host().push_back(testMaterial());
    particleContainer particles;
    for (int z = 0; z < 5; ++z)
    {
        for (int y = 0; y < 6; ++y)
        {
            for (int x = 0; x < 10; ++x)
            {
                particle value;
                value.setPosition({0.16 * x, 0.16 * y, 0.16 * z});
                value.setRadius(0.1);
                value.setMaterial(materials, 0);
                particles.host().push_back(value);
            }
        }
    }

    contactContainer directContacts;
    contactContainer gridContacts;
    contactSearch directSearch;
    contactSearch gridSearch;
    directSearch.setDirectSearch(true);
    require(directSearch.findContacts(directContacts, particles) == gridSearch.findContacts(gridContacts, particles),
            "The parallel grid search must find the same number of contacts as direct search.");

    const auto contactKeys = [](const contactContainer& contacts)
    {
        std::vector<long long> keys(contacts.hostSize());
        for (int contactIndex = 0; contactIndex < static_cast<int>(contacts.hostSize()); ++contactIndex)
        {
            const contact& value = contacts.host()[contactIndex];
            keys[contactIndex] = (static_cast<long long>(value.masterParticleIndex()) << 32) | static_cast<unsigned int>(value.slaveParticleIndex());
        }
        std::sort(keys.begin(), keys.end());
        return keys;
    };
    require(contactKeys(directContacts) == contactKeys(gridContacts), "The parallel grid search and direct search must return identical particle pairs.");
    require(gridSearch.occupiedCellCount() > 0, "The CPU spatial grid must contain occupied cells.");
}

void testMaterialTypeAndLevelSetAssignment()
{
    const material standardMaterial = testMaterial();
    LSMaterial levelSetMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0};
    require(standardMaterial.type() == materialType::standard, "material must use the standard material type.");
    require(levelSetMaterial.type() == materialType::levelSet, "LSMaterial must preserve its level-set marker when viewed through material.");
    const material& unifiedMaterial = levelSetMaterial;
    require(nearlyEqual(levelSetMaterial.slidingFrictionCoefficient(), 0.4) && unifiedMaterial.rollingFrictionCoefficient() == 0.0 && unifiedMaterial.torsionalFrictionCoefficient() == 0.0 &&
                unifiedMaterial.rollingStiffness() == 0.0 && unifiedMaterial.torsionalStiffness() == 0.0,
            "LSMaterial must expose sliding friction and keep unused rotational properties zero in unified storage.");

    materialContainer materials;
    materials.host().push_back(standardMaterial);
    materials.host().push_back(levelSetMaterial);

    LSParticle particle;
    particle.setMaterial(materials, 0);
    require(particle.materialIndex() == -1, "LSParticle must reject a non-level-set material.");
    particle.setMaterial(materials, 1);
    require(particle.materialIndex() == 1 && particle.getMaterial().isLevelSet(), "LSParticle must accept a level-set material.");
}

void testRotationalContactFriction()
{
    enum class pairKind { sphereLevelSet, sphereSphere, levelSetLevelSet, frictionlessSphereLevelSet };
    for (const pairKind kind : {pairKind::sphereLevelSet, pairKind::sphereSphere, pairKind::levelSetLevelSet, pairKind::frictionlessSphereLevelSet})
        for (const bool infiniteMass : {false, true})
        {
            const bool nodalContact = kind == pairKind::levelSetLevelSet;
            const bool mixedContact = kind == pairKind::sphereLevelSet || kind == pairKind::frictionlessSphereLevelSet;
            const bool frictionlessSphere = kind == pairKind::frictionlessSphereLevelSet;
            materialContainer materials;
            const LSMaterial levelSetMaterial{2.0e5, 8.0e4, 0.6, 1.0, 1200.0};
            materials.host().push_back(nodalContact ? static_cast<material>(levelSetMaterial)
                                                    : material{1.0e5, 5.0e4, 2.0e4, 1.0e4, 0.4, frictionlessSphere ? 0.0 : 0.3, frictionlessSphere ? 0.0 : 0.2, 1.0, 1000.0});
            materials.host().push_back(kind == pairKind::sphereSphere ? material{1.0e5, 5.0e4, 2.0e4, 1.0e4, 0.4, 0.9, 0.6, 1.0, 1000.0}
                                                                     : static_cast<material>(levelSetMaterial));
            particle master;
            master.setPosition({0.0, 0.0, 0.0});
            master.setAngularVelocity({100.0, 0.0, 100.0});
            master.setRadius(0.1);
            master.setMaterial(materials, 0);
            particle slave;
            slave.setPosition({0.15, 0.0, 0.0});
            slave.setRadius(0.1);
            slave.setMaterial(materials, 1);
            if (infiniteMass)
                slave.setInfiniteMass();
            particleContainer masters;
            particleContainer slaves;
            masters.host().push_back(master);
            slaves.host().push_back(slave);

            contact value;
            require(value.setMasterSlaveParticle(masters, 0, slaves, 0), "A contact must accept separate master/slave containers.");
            value.setPoint({0.075, 0.0, 0.0});
            value.setNormal({-1.0, 0.0, 0.0});
            value.setOverlap(0.05);
            value.setArea(2.0);
            value.setEffectiveMass(1.0);
            value.setEffectiveRadius(0.1);
            value.setRollingSpringDeformation({0.0, 0.1, 0.0});
            value.setTorsionalSpringDeformation({0.1, 0.0, 0.0});
            contactContainer contacts;
            contacts.host().push_back(value);
            require(contacts.host()[0].calculateForce(1.0), "A configured contact must calculate force.");
            const contact& evaluated = contacts.host()[0];
            const math::Real scale = nodalContact ? 2.0 : 1.0;
            const math::Real expectedStiffness = mixedContact ? materials.host()[0].normalStiffness()
                : execution::effectiveStiffness(scale * materials.host()[0].normalStiffness(), scale * materials.host()[1].normalStiffness(), false, infiniteMass);
            const math::Real normalForce = 0.05 * expectedStiffness;
            require(nearlyEqual(evaluated.normalForceMagnitude(), normalForce), "Normal contact stiffness changed with rotational friction selection.");
            const math::Real rollingFriction = nodalContact || frictionlessSphere ? 0.0 : mixedContact ? 0.3 : execution::harmonicMean(0.3, 0.9);
            const math::Real torsionalFriction = nodalContact || frictionlessSphere ? 0.0 : mixedContact ? 0.2 : execution::harmonicMean(0.2, 0.6);
            require(nearlyEqual(evaluated.torque().x, -0.2 * torsionalFriction * normalForce) && nearlyEqual(evaluated.torque().y, 0.0) &&
                        nearlyEqual(evaluated.torque().z, -0.1 * rollingFriction * normalForce),
                    "Sphere-LS rotational friction must use only the sphere; sphere-sphere must retain harmonic combination; LS-LS must have no direct rotational torque.");
            if (nodalContact || frictionlessSphere)
                require(evaluated.rollingSpringDeformation() == math::Vec3::zero() && evaluated.torsionalSpringDeformation() == math::Vec3::zero(),
                        "An inactive rotational friction law must clear both spring histories.");
        }
}

void testBondOwnershipAndValidity()
{
    SphereDEM solver;
    solver.addMaterial(testMaterial());
    for (int index = 0; index < 2; ++index)
    {
        particle value;
        value.setPosition({0.2 * index, 0.0, 0.0});
        value.setRadius(0.1);
        value.setMaterial(solver.materials(), 0);
        solver.addSphere(value);
    }

    bond value;
    require(value.setEquivalentLength(0.2), "A positive bond equivalent length must be accepted.");
    require(value.setConnection(solver.spheres(), 0, 1, {1.0, 0.0, 0.0}), "A valid sphere bond connection must be accepted.");
    require(!value.isValid(), "A bond without stiffness must not be valid.");
    require(value.setStiffness(1.0e5, 5.0e4, 1.0e3, 1.0e3), "Finite non-negative bond stiffnesses must be accepted.");
    require(value.isValid(), "A fully configured bond must be valid.");
    require(value.setEquivalentLength(0.25), "A changed positive bond equivalent length must be accepted.");
    require(!value.isValid(), "Changing the equivalent length must invalidate dependent bond geometry and stiffness.");
    require(value.coefficientB1() == 0.0 && value.coefficientB2() == 0.0 && value.coefficientB3() == 0.0 && value.coefficientB4() == 0.0,
            "Changing the equivalent length must clear the old bond coefficients.");
    require(value.setConnection(solver.spheres(), 0, 1, {1.0, 0.0, 0.0}), "Bond geometry must be configurable again after changing equivalent length.");
    require(value.setStiffness(1.0e5, 5.0e4, 1.0e3, 1.0e3), "Bond stiffness must be configurable again after changing equivalent length.");
    require(value.isValid(), "Reconfiguring length-dependent bond values must restore validity.");
    require(solver.addBond(value) == 0, "A bond referencing the solver sphere container must be accepted.");
}

void testBondCrossSectionArea()
{
    bond connection;
    require(connection.crossSectionArea() == 0.0, "The default bond cross-sectional area must disable fracture.");
    require(connection.setCrossSectionArea(0.02) && nearlyEqual(connection.crossSectionArea(), 0.02), "A positive finite cross-sectional area must be accepted.");
    require(!connection.setCrossSectionArea(-1.0) && !connection.setCrossSectionArea(std::numeric_limits<math::Real>::infinity()) && nearlyEqual(connection.crossSectionArea(), 0.02),
            "Invalid cross-sectional areas must not overwrite the last valid value.");
    require(connection.setCrossSectionArea(0.0), "A zero cross-sectional area must remain supported.");
    math::Real damage = 0.25;
    math::Real maximumRatio = 0.7;
    require(!execution::updateBKDamage(damage, maximumRatio, 1.0, 10.0, 10.0, 10.0, 10.0, connection.crossSectionArea(), 1.0, 1.0, 1.75, 0.9) &&
                damage == 0.25 && maximumRatio == 0.7,
            "Zero cross-sectional area must disable further fracture without erasing existing damage history.");
}

void testModeICompressionDoesNotDamageBond()
{
    math::Real compressionDamage = 0.0;
    math::Real compressionMaximumRatio = 0.0;
    const bool compressionFailed = execution::updateBKDamage(compressionDamage, compressionMaximumRatio, -0.1, 2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    require(!compressionFailed && compressionDamage == 0.0 && compressionMaximumRatio == 0.0, "Pure axial compression must not drive mode-I bond damage.");

    math::Real tensionDamage = 0.0;
    math::Real tensionMaximumRatio = 0.0;
    const bool tensionFailed = execution::updateBKDamage(tensionDamage, tensionMaximumRatio, 0.1, 2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    require(tensionFailed && tensionDamage == 1.0, "Axial opening above GIc must completely break the bond.");
}

void testMaterialDensityValidation()
{
    material value;
    bool rejected = false;
    try
    {
        value.setDensity(0.0);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "Zero material density must be rejected.");
    require(!value.isValid(), "Rejecting zero density must not make the material valid.");
}

void testRotatingFixedWallAcceleration()
{
    math::Vec3 velocity = math::Vec3::zero();
    math::Vec3 acceleration = math::Vec3::zero();
    execution::accumulateVirtualParticleVelocityAndAcceleration(velocity, acceleration, math::Vec3::unitX(), math::Quaternion::identity(), math::Vec3::zero(), {0.0, 0.0, 2.0});
    require(nearlyEqual(velocity.x, 0.0) && nearlyEqual(velocity.y, 2.0) && nearlyEqual(velocity.z, 0.0), "A rotating fixed wall must provide the virtual-particle tangential velocity.");
    require(nearlyEqual(acceleration.x, -4.0) && nearlyEqual(acceleration.y, 0.0) && nearlyEqual(acceleration.z, 0.0),
            "A rotating fixed wall must provide the virtual-particle centripetal acceleration.");
}

void testSystemSPHProperties()
{
    constexpr math::Real spacing = 0.02;
    constexpr math::Real smoothingLength = 0.026;
    constexpr math::Real referenceDensity = 1000.0;
    SPHDEM solver;
    solver.setSPHProperties(spacing, smoothingLength, referenceDensity, 1.0e-3);
    const SPHParticle particle{{0.0, 0.0, 0.0}};
    require(solver.addSPHParticle(particle) == 0, "A valid SPH particle must be accepted after system properties are set.");
    require(nearlyEqual(solver.SPHParticles().host()[0].mass(), referenceDensity * spacing * spacing * spacing), "SPH mass must be derived from the system spacing and reference density.");
    require(nearlyEqual(solver.SPHParticles().host()[0].density(), referenceDensity), "SPH density must be initialized from the system reference density.");
    require(nearlyEqual(solver.SPHSmoothingLength(), smoothingLength), "The smoothing length must be stored once by SPHDEM.");
    require(execution::calculateLatticeKernelSum3D(spacing, smoothingLength) > 0.0, "The system reference kernel sum must be finite and positive.");
}

void testSPHWallReaction(executionMode mode)
{
    constexpr math::Real spacing = 0.02;
    SPHDEM simulation{mode};
    simulation.setBoundary({-0.2, -0.2, -0.1}, {0.2, 0.2, 0.2});
    simulation.setGravity({0.0, 0.0, -9.81});
    simulation.setTimeStep(1.0e-5);
    simulation.setSPHProperties(spacing, 1.3 * spacing, 1000.0, 1.0e-3);
    simulation.addSPHParticle(SPHParticle{{0.0, 0.0, 0.015}});

    const int materialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, 0.0, 1.0, 1.0});
    levelset::PlaneWall plane{math::Vec3::unitZ(), 0.2};
    plane.buildLSGrid(spacing, 3);
    const int geometryIndex = simulation.addGeometry(plane, true);
    LSParticle wall;
    wall.setMaterial(simulation.materials(), materialIndex);
    wall.setGeometry(simulation.geometries(), geometryIndex);
    simulation.addLSParticle(wall);

    simulation.solve(0);
    require(simulation.virtualParticles().hostSize() > 0, "A level-set SPH wall must generate virtual particles.");
    const math::Vec3 fluidForce = simulation.SPHParticles().host()[0].force();
    const math::Vec3 wallForce = simulation.LSParticles().host()[0].force();
    const math::Real forceScale = math::norm(fluidForce) + math::norm(wallForce);
    require(forceScale > math::defaultTolerance, "Gravity-induced wall pressure must generate a non-zero SPH wall reaction.");
    require(math::norm(fluidForce + wallForce) <= 1.0e-11 * forceScale, "SPH wall force must be equal and opposite to the fluid force contribution.");
}

void testSPHKernelAndMomentumConservation()
{
    constexpr math::Real smoothingLength = 0.03;
    constexpr math::Real epsilon = 1.0e-6 * smoothingLength;
    const math::Vec3 separation{0.012, -0.006, 0.003};
    const math::Vec3 analyticalGradient = execution::wendlandKernelGradient3D(separation, smoothingLength);
    const auto numericalGradient = [&](const math::Vec3& direction)
    {
        return (execution::wendlandKernel3D(math::norm(separation + epsilon * direction), smoothingLength) -
                execution::wendlandKernel3D(math::norm(separation - epsilon * direction), smoothingLength)) /
               (2.0 * epsilon);
    };
    require(nearlyEqual(analyticalGradient.x, numericalGradient(math::Vec3::unitX()), 1.0e-6) && nearlyEqual(analyticalGradient.y, numericalGradient(math::Vec3::unitY()), 1.0e-6) &&
                nearlyEqual(analyticalGradient.z, numericalGradient(math::Vec3::unitZ()), 1.0e-6),
            "The Wendland kernel gradient must match the derivative of the kernel value.");

    constexpr math::Real particleMass = 0.008;
    constexpr math::Real dynamicViscosity = 1.0e-3;
    constexpr math::Real soundSpeed = 10.0;
    const math::Vec3 firstPosition{0.0, 0.0, 0.0};
    const math::Vec3 secondPosition{0.02, 0.0, 0.0};
    const math::Vec3 firstVelocity{0.1, 0.02, 0.0};
    const math::Vec3 secondVelocity{-0.1, -0.01, 0.0};
    constexpr math::Real firstDensity = 1010.0;
    constexpr math::Real secondDensity = 990.0;
    const math::Real firstPressure = execution::computePressureFromDensity(firstDensity, 1000.0, soundSpeed);
    const math::Real secondPressure = execution::computePressureFromDensity(secondDensity, 1000.0, soundSpeed);
    const math::Vec3 firstAcceleration =
        execution::calculateFluidPressureAcceleration(firstPosition,
                                                      firstVelocity,
                                                      firstDensity,
                                                      firstPressure,
                                                      secondPosition,
                                                      secondVelocity,
                                                      secondDensity,
                                                      secondPressure,
                                                      particleMass,
                                                      1000.0,
                                                      smoothingLength,
                                                      soundSpeed) +
        execution::calculateFluidViscousAcceleration(firstPosition, firstVelocity, firstDensity, secondPosition, secondVelocity, secondDensity, particleMass, dynamicViscosity, smoothingLength);
    const math::Vec3 secondAcceleration =
        execution::calculateFluidPressureAcceleration(secondPosition,
                                                      secondVelocity,
                                                      secondDensity,
                                                      secondPressure,
                                                      firstPosition,
                                                      firstVelocity,
                                                      firstDensity,
                                                      firstPressure,
                                                      particleMass,
                                                      1000.0,
                                                      smoothingLength,
                                                      soundSpeed) +
        execution::calculateFluidViscousAcceleration(secondPosition, secondVelocity, secondDensity, firstPosition, firstVelocity, firstDensity, particleMass, dynamicViscosity, smoothingLength);
    const math::Vec3 forceSum = particleMass * (firstAcceleration + secondAcceleration);
    const math::Real forceScale = particleMass * (math::norm(firstAcceleration) + math::norm(secondAcceleration));
    require(math::norm(forceSum) <= 1.0e-10 * (1.0 + forceScale), "A fluid-particle pair must conserve linear momentum.");

    math::Real density = 1000.0;
    math::Real pressure = 0.0;
    execution::integrateDensityAndPressure(density, pressure, -std::numeric_limits<math::Real>::infinity(), 1000.0, soundSpeed, 1.0);
    require(nearlyEqual(density, execution::numericalSPHDensityFloorRatio * 1000.0), "A failed density update must use the reference-density-scaled numerical floor.");
}

} // namespace

int main()
{
    try
    {
        testSphereDEMExecutionModes();
        testSPHDEMExecutionModes();
        testContainerOwnership();
        testSplitSolveContinuation();
        testAppendAfterSolve();
        testLevelSetContainerOwnership();
        testLSParticleMotionReset(executionMode::CPU);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testLSParticleMotionReset(executionMode::GPU);
#endif
        testContactCalculationGate();
        testParallelContactForceAssembly();
        testParallelGridContactSearch();
        testMaterialTypeAndLevelSetAssignment();
        testRotationalContactFriction();
        testBondOwnershipAndValidity();
        testBondCrossSectionArea();
        testModeICompressionDoesNotDamageBond();
        testMaterialDensityValidation();
        testRotatingFixedWallAcceleration();
        testSystemSPHProperties();
        testSPHWallReaction(executionMode::CPU);
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        testSPHWallReaction(executionMode::Hybrid);
        testSPHWallReaction(executionMode::GPU);
#endif
        testSPHKernelAndMomentumConservation();
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
