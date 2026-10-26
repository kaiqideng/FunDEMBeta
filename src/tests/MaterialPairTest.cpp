/**
 * @file MaterialPairTest.cpp
 * @brief Checks bulk material validation, explicit coefficient switches, and stable symmetric pair storage.
 */
#include "material/materialPair.h"
#include "data/myLSObject.h"
#include "solver/SphereDEM.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{

using namespace fundem;

void require(bool condition, const char* message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

template <class Exception = std::invalid_argument, class Function> void requireRejected(Function function, const char* message)
{
    try
    {
        function();
    }
    catch (const Exception&)
    {
        return;
    }
    throw std::runtime_error(message);
}

template <class Layout> struct fieldCount;
template <class... Fields> struct fieldCount<deviceLayout<Fields...>> : std::integral_constant<int, sizeof...(Fields)> {};
static_assert(fieldCount<material::DeviceLayout>::value == 4, "Bulk material device storage must contain exactly four fields.");
static_assert(fieldCount<materialPair::DeviceLayout>::value == 17, "Pair device storage must contain thirteen coefficients, two switches, and two material IDs.");
static_assert(std::is_trivial_v<materialPairCoefficients> && std::is_standard_layout_v<materialPairCoefficients>, "Coefficient transfer data must remain a plain aggregate.");

void testMaterial()
{
    material value;
    require(!value.isValid(), "A default material must be incomplete.");
    value.setYoungModulus(2.0e7);
    require(!value.isValid(), "A material without density must remain incomplete.");
    value.setDensity(2500.0);
    require(value.isValid(), "Assigning modulus and density must complete the material.");
    value.setProperties(3.0e7, -0.2, 0.9, 2600.0);
    require(value.youngModulus() == 3.0e7 && value.poissonRatio() == -0.2 && value.restitutionCoefficient() == 0.9 && value.density() == 2600.0, "All four material properties must be retained.");
    requireRejected([&] { value.setProperties(9.0e7, 0.5, 0.5, 3000.0); }, "Incompressible Poisson ratio must be rejected.");
    require(value.youngModulus() == 3.0e7 && value.density() == 2600.0, "Rejected bulk replacement must not partially update the material.");
    requireRejected([&] { value.setPoissonRatio(-1.0); }, "Poisson ratio -1 must be rejected.");
    requireRejected([&] { value.setYoungModulus(0.0); }, "Zero Young's modulus must be rejected.");
    requireRejected([&] { value.setDensity(-1.0); }, "Negative density must be rejected.");
    requireRejected([&] { value.setRestitutionCoefficient(1.1); }, "Restitution above unity must be rejected.");
    requireRejected([&] { value.setRestitutionCoefficient(std::numeric_limits<math::Real>::quiet_NaN()); }, "NaN restitution must be rejected.");
    value.setRestitutionCoefficient(0.0);
    value.setRestitutionCoefficient(1.0);
}

void testCoefficientSwitches()
{
    materialPair pair(3, 1);
    require(pair.masterMaterialIndex() == 1 && pair.slaveMaterialIndex() == 3, "Pair indices must be canonical and symmetric.");
    require(!pair.linearStiffnessEnabled() && !pair.dampingEnabled() && !pair.hasLevelSetStiffness(),
            "Exactly two contact-law switches must initially be off; LS coefficients must initially be unassigned.");
    require(pair.coefficients().normalStiffness_ == 0.0 && pair.coefficients().normalDampingCoefficient_ == 0.0, "Ordinary stiffness and damping values must default to zero.");
    require(pair.normalStiffnessPerUnitArea() == -1.0 && pair.slidingStiffnessPerUnitArea() == -1.0, "Unassigned LS coefficients must retain their internal sentinel.");

    pair.setLevelSetStiffness(15.0, 16.0);
    require(pair.hasLevelSetStiffness() && !pair.linearStiffnessEnabled() && !pair.dampingEnabled(), "Assigning LS stiffness must not enable either contact-law switch.");
    pair.setNormalStiffnessPerUnitArea(0.0);
    pair.setSlidingStiffnessPerUnitArea(0.0);
    require(pair.hasLevelSetStiffness() && !pair.linearStiffnessEnabled() && !pair.dampingEnabled(), "Explicit zero LS stiffness must be valid without a third switch.");

    pair.setLinearStiffness(11.0, 12.0, 13.0, 14.0);
    require(pair.linearStiffnessEnabled() && !pair.dampingEnabled(), "Setting linear stiffness must enable only its group.");
    pair.setLinearStiffnessEnabled(false);
    require(!pair.linearStiffnessEnabled() && pair.normalStiffness() == 11.0 && pair.rollingStiffness() == 13.0, "Disabling the linear model must retain every stored stiffness.");
    pair.setLinearStiffnessEnabled(true);
    require(pair.normalStiffness() == 11.0 && pair.slidingStiffness() == 12.0 && pair.rollingStiffness() == 13.0 && pair.torsionalStiffness() == 14.0,
            "Re-enabling the linear model must restore the configured group values.");
    pair.setLinearStiffnessEnabled(false);
    pair.setRollingStiffness(7.0);
    require(pair.linearStiffnessEnabled() && pair.normalStiffness() == 11.0 && pair.rollingStiffness() == 7.0,
            "An individual ordinary stiffness setter must enable the complete linear group, not a per-direction mode.");
    pair.setLinearStiffness(0.0, 0.0, 0.0, 0.0);
    require(pair.linearStiffnessEnabled() && pair.normalStiffness() == 0.0, "Enabled zero linear stiffness must differ from Hertz-Mindlin fallback.");

    pair.setDampingCoefficients(1.0, 2.0, 3.0, 4.0);
    require(pair.dampingEnabled(), "Setting explicit damping must enable its entire group.");
    pair.setDampingEnabled(false);
    require(!pair.dampingEnabled() && pair.normalDampingCoefficient() == 1.0 && pair.torsionalDampingCoefficient() == 4.0, "Disabling explicit damping must retain all four values.");
    pair.setDampingEnabled(true);
    require(pair.slidingDampingCoefficient() == 2.0 && pair.rollingDampingCoefficient() == 3.0, "Re-enabling explicit damping must preserve the stored group.");
    pair.setDampingEnabled(false);
    pair.setNormalDampingCoefficient(0.0);
    require(pair.dampingEnabled() && pair.normalDampingCoefficient() == 0.0 && pair.slidingDampingCoefficient() == 2.0,
            "An individual damping setter must enable the complete explicit damping group without changing other values.");
    pair.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
    require(pair.dampingEnabled(), "Enabled zero dashpots must explicitly suppress restitution-based damping.");
    pair.setFrictionCoefficients(0.1, 0.2, 0.3);
    const auto saved = pair.coefficients();
    requireRejected([&] { pair.setLinearStiffness(100.0, -1.0); }, "Negative linear stiffness must be rejected.");
    requireRejected([&] { pair.setDampingCoefficients(100.0, 200.0, 300.0, -1.0); }, "Negative damping must be rejected.");
    requireRejected([&] { pair.setFrictionCoefficients(100.0, std::numeric_limits<math::Real>::infinity()); }, "Infinite friction must be rejected.");
    requireRejected([&] { pair.setLevelSetStiffness(100.0, -1.0); }, "Users must not assign the internal LS sentinel.");
    require(pair.normalStiffness() == saved.normalStiffness_ && pair.normalDampingCoefficient() == saved.normalDampingCoefficient_ &&
                pair.slidingFrictionCoefficient() == saved.slidingFrictionCoefficient_ && pair.normalStiffnessPerUnitArea() == saved.normalStiffnessPerUnitArea_ &&
                pair.linearStiffnessEnabled() == saved.linearStiffnessEnabled_ && pair.dampingEnabled() == saved.dampingEnabled_,
            "Rejected group setters must leave values and both switches unchanged.");

    pair.setLinearStiffnessEnabled(false);
    pair.setDampingEnabled(false);
    pair.setFrictionCoefficients(0.4);
    pair.setNormalStiffnessPerUnitArea(21.0);
    require(!pair.linearStiffnessEnabled() && !pair.dampingEnabled(), "Friction and LS coefficient changes must never activate either model switch.");
}

void testStableStorage()
{
    requireRejected([] { materialPair invalid(-1, 0); }, "Negative material IDs must be rejected.");
    materialPairContainer pairs(2);
    require(pairs.materialCount() == 2 && pairs.hostSize() == 3, "Two materials require three pairs including self-pairs.");
    require(&pairs.at(0, 1) == &pairs.at(1, 0), "Both pair orderings must identify the same record.");
    pairs.at(1, 0).setLinearStiffness(100.0, 200.0);
    pairs.at(1, 1).setDampingCoefficients(1.0, 2.0, 3.0, 4.0);
    pairs.resizeForMaterials(5);
    require(pairs.hostSize() == 15 && pairs.at(0, 1).normalStiffness() == 100.0 && pairs.at(1, 1).torsionalDampingCoefficient() == 4.0,
            "Growing the table must preserve all previously configured pairs.");
    int expected = 0;
    const materialPairContainer& constPairs = pairs;
    for (int larger = 0; larger < 5; ++larger)
    {
        for (int smaller = 0; smaller <= larger; ++smaller)
        {
            require(materialPairContainer::index(larger, smaller) == expected && materialPairContainer::index(smaller, larger) == expected,
                    "Triangular indices must be symmetric, contiguous, and stable.");
            const auto& pair = constPairs.at(larger, smaller);
            require(pair.masterMaterialIndex() == smaller && pair.slaveMaterialIndex() == larger, "Every triangular slot must retain the correct material IDs.");
            ++expected;
        }
    }
    pairs.resizeForMaterials(5);
    require(pairs.hostSize() == 15, "An unchanged material count must not append duplicate records.");
    requireRejected([&] { pairs.resizeForMaterials(4); }, "Append-only pair storage must reject material removal.");
    requireRejected<std::out_of_range>([&] { pairs.at(5, 0); }, "Out-of-range material IDs must be rejected.");
    requireRejected<std::out_of_range>([&] { constPairs.at(0, -1); }, "Negative lookup IDs must be rejected.");
}

void testMovedStorage()
{
    materialPairContainer source(2);
    source.at(0, 1).setLinearStiffness(101.0, 202.0);
    materialPairContainer target(std::move(source));
    require(target.materialCount() == 2 && target.hostSize() == 3 && target.at(1, 0).normalStiffness() == 101.0, "Move construction must transfer the complete pair table.");
    require(source.materialCount() == 0 && source.hostSize() == 0, "Move construction must leave a consistent empty source.");
    requireRejected<std::out_of_range>([&] { source.at(0, 0); }, "A moved-from table must reject lookups instead of accessing empty storage.");
    source.resizeForMaterials(1);
    source.at(0, 0).setNormalDampingCoefficient(3.0);
    require(source.materialCount() == 1 && source.at(0, 0).normalDampingCoefficient() == 3.0, "A moved-from table must support reuse.");

    materialPairContainer assigned(3);
    assigned = std::move(target);
    require(assigned.materialCount() == 2 && assigned.hostSize() == 3 && assigned.at(0, 1).slidingStiffness() == 202.0, "Move assignment must replace the destination table.");
    require(target.materialCount() == 0 && target.hostSize() == 0, "Move assignment must leave a consistent empty source.");
    target.resizeForMaterials(1);
    require(target.hostSize() == 1 && !target.at(0, 0).linearStiffnessEnabled() && !target.at(0, 0).dampingEnabled(), "A move-assigned source must support fresh initialization.");
    materialPairContainer& sameTable = assigned;
    assigned = std::move(sameTable);
    require(assigned.materialCount() == 2 && assigned.at(0, 1).normalStiffness() == 101.0, "Self-move assignment must preserve a usable table.");
}

void testMaterialReplacement()
{
    SphereDEM simulation;
    simulation.setBoundary({-4.0, -4.0, -4.0}, {4.0, 4.0, 4.0});
    simulation.setGravity(math::Vec3::zero());
    simulation.setTimeStep(1.0e-6);
    simulation.addMaterial(material{1.0e7, 0.2, 0.7, 1000.0});
    simulation.addMaterial(material{1.0e7, 0.2, 0.7, 1000.0});
    materialPair pair(0, 0);
    pair.setLevelSetStiffness(1000.0, 500.0);
    pair.setLinearStiffness(4000.0, 2000.0, 100.0, 50.0);
    pair.setFrictionCoefficients(10.0, 10.0, 10.0);
    pair.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
    simulation.setMaterialPair(pair);
    for (int index = 0; index < 3; ++index)
    {
        particle value;
        value.setRadius(0.1);
        value.setMaterial(simulation.materials(), index == 2 ? 1 : 0);
        value.setPosition({index == 2 ? 3.0 : 0.15 * index, 0.0, 0.0});
        if (index == 0)
            value.setVelocity({0.0, 0.01, 0.0});
        if (index == 1)
            value.setInfiniteMass();
        simulation.addSphere(value);
    }
    levelset::Sphere geometry(0.1);
    geometry.buildSurfaceNode(0);
    geometry.buildLSGrid(0.05, 2, false);
    simulation.addGeometry(geometry);
    for (int index = 0; index < 2; ++index)
    {
        LSParticle value;
        value.setMaterial(simulation.materials(), 0);
        value.setGeometry(simulation.geometries(), 0);
        value.setPosition({-1.0 - index, 0.0, 0.0});
        if (index == 1)
            value.setInfiniteMass();
        simulation.addLSParticle(value);
    }
    bond connection;
    require(connection.setEquivalentLength(0.15) && connection.setConnection(simulation.spheres(), 0, 1) &&
                connection.setStiffness(100.0, 50.0, 10.0, 5.0),
            "The material replacement fixture requires a connected bond.");
    simulation.addBond(connection);
    simulation.solve(1);
    require(simulation.sphereInteractions().contacts().hostSize() == 1, "The material replacement fixture requires one active contact.");

    const auto spheresBefore = simulation.spheres().host();
    const auto levelSetsBefore = simulation.LSParticles().host();
    require(levelSetsBefore.front().inverseMass() > 0.0 && levelSetsBefore.back().inverseMass() == 0.0,
            "The material replacement fixture must include both finite and fixed LS particles.");
    const contact contactBefore = simulation.sphereInteractions().contacts().host().front();
    const bond bondBefore = simulation.sphereInteractions().bonds().host().front();
    const auto timeBefore = simulation.time();
    const int stepBefore = simulation.stepCount();
    const material replacement{3.0e7, 0.3, 0.5, 2000.0};
    requireRejected<std::out_of_range>([&] { simulation.setMaterial(-1, replacement); }, "Material replacement must reject negative indices.");
    requireRejected<std::out_of_range>([&] { simulation.setMaterial(2, replacement); }, "Material replacement must reject missing indices.");
    requireRejected([&] { simulation.setMaterial(0, material{}); }, "Material replacement must reject incomplete definitions.");
    require(simulation.materials().host()[0].density() == 1000.0, "Rejected replacements must preserve the prior definition.");
    simulation.initialize();
    requireRejected<std::logic_error>([&] { simulation.setMaterial(0, replacement); }, "Initialized solvers must reject material mutation.");
    simulation.solve(0);
    simulation.setMaterial(0, replacement);

    const auto near = [](math::Real actual, math::Real expected, const char* message)
    {
        require(std::isfinite(actual) && std::abs(actual - expected) <= 1.0e-11 * (1.0 + std::abs(expected)), message);
    };
    const auto unchangedMotion = [](const particle& actual, const particle& before)
    {
        require(actual.position() == before.position() && actual.orientation() == before.orientation() &&
                    actual.velocity() == before.velocity() && actual.angularVelocity() == before.angularVelocity(),
                "Replacing a material must preserve all rigid-body kinematics.");
    };
    const auto refreshed = [&](const particle& actual, const particle& before)
    {
        unchangedMotion(actual, before);
        require(actual.getMaterial().youngModulus() == replacement.youngModulus() && actual.getMaterial().poissonRatio() == replacement.poissonRatio() &&
                    actual.getMaterial().restitutionCoefficient() == replacement.restitutionCoefficient() && actual.getMaterial().density() == replacement.density(),
                "Every particle using the replaced material must refresh its host material snapshot.");
        near(actual.inverseMass(), before.inverseMass() / 2.0, "Doubling density must double finite particle mass and preserve fixed particles.");
        for (int component = 0; component < 9; ++component)
        {
            near(actual.inertiaTensor()[component], before.inertiaTensor()[component] * 2.0, "Density replacement must rescale the full inertia tensor.");
            near(actual.inverseInertiaTensor()[component], before.inverseInertiaTensor()[component] / 2.0, "Density replacement must rescale inverse inertia.");
        }
    };
    for (int index = 0; index < 2; ++index)
    {
        refreshed(simulation.spheres().host()[index], spheresBefore[index]);
        refreshed(simulation.LSParticles().host()[index], levelSetsBefore[index]);
    }
    unchangedMotion(simulation.spheres().host()[2], spheresBefore[2]);
    require(simulation.spheres().host()[2].getMaterial().density() == 1000.0 &&
                simulation.spheres().host()[2].inverseMass() == spheresBefore[2].inverseMass(),
            "Replacing one material must not change particles using another material.");
    const contact& retained = simulation.sphereInteractions().contacts().host().front();
    require(retained.slidingSpringDeformation() == contactBefore.slidingSpringDeformation() &&
                retained.rollingSpringDeformation() == contactBefore.rollingSpringDeformation() &&
                retained.torsionalSpringDeformation() == contactBefore.torsionalSpringDeformation(),
            "Replacing a material must retain contact spring history.");
    const bond& retainedBond = simulation.sphereInteractions().bonds().host().front();
    require(retainedBond.damageFactor() == bondBefore.damageFactor() && retainedBond.maximumEnergyReleaseRatio() == bondBefore.maximumEnergyReleaseRatio() &&
                retainedBond.masterEndpointLocalPosition() == bondBefore.masterEndpointLocalPosition() &&
                retainedBond.slaveEndpointLocalPosition() == bondBefore.slaveEndpointLocalPosition(),
            "Replacing a material must retain bond damage and reference history.");
    pair.setNormalStiffness(8000.0);
    simulation.setMaterialPair(pair);
    simulation.initialize();
    require(simulation.stepCount() == stepBefore && simulation.time() == timeBefore, "Material update initialization must not advance or restart simulation time.");
    const contact& updated = simulation.sphereInteractions().contacts().host().front();
    near(updated.effectiveMass(), contactBefore.effectiveMass() * 2.0, "Reinitialization must refresh contact effective mass after density changes.");
    near(updated.normalForceMagnitude(), contactBefore.normalForceMagnitude() * 2.0, "Reinitialization must use the updated pair stiffness immediately.");
    require(math::nearlyEqual(updated.slidingSpringDeformation(), contactBefore.slidingSpringDeformation(), 1.0e-12, 1.0e-12),
            "Force refresh must preserve the existing contact spring when it remains below the friction limit.");
    simulation.solve(1);
    require(simulation.stepCount() == stepBefore + 1 && simulation.time() > timeBefore, "Continuation must advance from the existing solver state.");
}

} // namespace

int main()
{
    try
    {
        testMaterial();
        testCoefficientSwitches();
        testStableStorage();
        testMovedStorage();
        testMaterialReplacement();
        std::cout << "Material and material-pair tests passed.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
