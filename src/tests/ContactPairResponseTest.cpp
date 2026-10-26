/** Verifies material-pair defaults, activation, Hertz energy and LS nodal damping. */
#include "data/myLSObject.h"
#include "execution/contactMaterialFunctions.h"
#include "solver/LSDEM.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using namespace fundem;
using execution::contactType;
using math::Real;

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void near(Real actual, Real expected, const char* message) { require(std::isfinite(actual) && std::abs(actual - expected) <= 1.0e-11 * (1.0 + std::abs(expected)), message); }

execution::contactResponse response(const materialPair& pair, contactType type = contactType::sphereSphere, Real penetration = 0.002, Real area = 0.4, Real weight = 1.0)
{
    return execution::resolveContactResponse(pair.coefficients(), type, 1.2e7, 0.2, 0.7, 8.0e6, 0.3, 0.5, 1.25 * weight, 0.15, penetration, area, weight);
}

void testHertzMindlin()
{
    materialPair pair(0, 0);
    const Real effectiveYoung = 1.0 / ((1.0 - 0.2 * 0.2) / 1.2e7 + (1.0 - 0.3 * 0.3) / 8.0e6);
    const Real firstShear = 1.2e7 / (2.0 * 1.2);
    const Real secondShear = 8.0e6 / (2.0 * 1.3);
    const Real effectiveShear = 1.0 / (1.8 / firstShear + 1.7 / secondShear);
    const Real a = std::sqrt(0.15 * 0.002);
    const Real force = (4.0 / 3.0) * effectiveYoung * a * 0.002;
    const Real e = 2.0 / (1.0 / 0.7 + 1.0 / 0.5);
    const Real logarithm = std::log(e);
    const Real beta = -logarithm / std::sqrt(logarithm * logarithm + math::pi * math::pi);
    for (contactType type : {contactType::sphereSphere, contactType::sphereLevelSet})
    {
        const auto result = response(pair, type);
        near(result.normalStiffness_ * 0.002, force, "Hertz normal force does not match both elastic compliances.");
        near(result.slidingStiffness_, 8.0 * effectiveShear * a, "Mindlin tangential stiffness is incorrect.");
        near(result.rollingStiffness_, result.slidingStiffness_, "Default Hertz-Mindlin rolling stiffness must equal sliding stiffness.");
        near(result.torsionalStiffness_, result.slidingStiffness_, "Default Hertz-Mindlin torsional stiffness must equal sliding stiffness.");
        near(result.normalElasticEnergy_, 0.4 * force * 0.002, "Hertz energy must be 2/5 Fn*overlap, not the linear 1/2 expression.");
        near(result.normalDampingCoefficient_, 2.0 * std::sqrt(5.0 / 6.0) * beta * std::sqrt(1.25 * 2.0 * effectiveYoung * a), "Hertz damping must use the normal tangent stiffness.");
        near(result.slidingDampingCoefficient_, 2.0 * std::sqrt(5.0 / 6.0) * beta * std::sqrt(1.25 * 8.0 * effectiveShear * a), "Mindlin restitution damping is incorrect.");
        const Real change = 1.0e-7;
        const Real derivative = (response(pair, type, 0.002 + change).normalElasticEnergy_ - response(pair, type, 0.002 - change).normalElasticEnergy_) / (2.0 * change);
        require(std::abs(derivative - force) < 1.0e-8 * force, "Normal energy derivative must equal the elastic force.");
    }
}

void testGroupActivation()
{
    for (contactType type : {contactType::sphereSphere, contactType::sphereLevelSet})
    {
        materialPair pair(0, 1);
        const auto defaults = response(pair, type);
        pair.setLinearStiffness(70000.0, 25000.0, 100.0, 50.0);
        auto result = response(pair, type);
        near(result.normalStiffness_, 70000.0, "Explicit pair normal stiffness must not be mixed a second time.");
        near(result.slidingStiffness_, 25000.0, "The linear group must use the configured sliding stiffness.");
        near(result.rollingStiffness_, 100.0, "The linear group must use the configured rolling stiffness.");
        near(result.torsionalStiffness_, 50.0, "The linear group must use the configured torsional stiffness.");
        near(result.normalElasticEnergy_, 0.5 * 70000.0 * 0.002 * 0.002, "Explicit linear normal energy changed.");
        pair.setLinearStiffnessEnabled(false);
        result = response(pair, type);
        near(result.normalStiffness_, defaults.normalStiffness_, "Disabling the linear group must restore Hertz normal stiffness.");
        near(result.slidingStiffness_, defaults.slidingStiffness_, "Disabling the linear group must restore Mindlin sliding stiffness.");
        near(result.rollingStiffness_, defaults.slidingStiffness_, "Default rolling stiffness must follow Mindlin sliding stiffness.");
        near(result.torsionalStiffness_, defaults.slidingStiffness_, "Default torsional stiffness must follow Mindlin sliding stiffness.");
        near(pair.normalStiffness(), 70000.0, "Disabling the linear group erased its configured value.");
        near(pair.rollingStiffness(), 100.0, "Disabling the linear group erased its configured rolling value.");
        pair.setDampingCoefficients(10.0, 20.0, 30.0, 40.0);
        result = response(pair, type);
        near(result.normalDampingCoefficient_, 10.0, "Explicit normal damping must apply alongside Hertz stiffness.");
        near(result.slidingDampingCoefficient_, 20.0, "Explicit sliding damping must apply alongside Mindlin stiffness.");
        near(result.rollingDampingCoefficient_, 30.0, "The damping group must activate rolling damping.");
        near(result.torsionalDampingCoefficient_, 40.0, "The damping group must activate torsional damping.");
        pair.setDampingEnabled(false);
        result = response(pair, type);
        near(result.normalDampingCoefficient_, defaults.normalDampingCoefficient_, "The damping switch must restore the normal restitution law.");
        near(result.slidingDampingCoefficient_, defaults.slidingDampingCoefficient_, "The damping switch must restore the sliding restitution law.");
        near(result.rollingDampingCoefficient_, defaults.rollingDampingCoefficient_, "The damping switch must restore the rolling restitution law.");
        near(result.torsionalDampingCoefficient_, defaults.torsionalDampingCoefficient_, "The damping switch must restore the torsional restitution law.");
        near(pair.torsionalDampingCoefficient(), 40.0, "Disabling damping must preserve every configured coefficient.");
        pair.setLinearStiffnessEnabled(true);
        pair.setDampingCoefficients(0.0, 0.0, 0.0, 0.0);
        result = response(pair, type);
        near(result.normalDampingCoefficient_, 0.0, "Enabled zero normal damping must override restitution.");
        near(result.slidingDampingCoefficient_, 0.0, "Enabled zero sliding damping must override restitution.");
        near(result.rollingDampingCoefficient_, 0.0, "Enabled zero rolling damping must override restitution.");
        near(result.torsionalDampingCoefficient_, 0.0, "Enabled zero torsional damping must override restitution.");
        pair.setDampingEnabled(false);
        result = response(pair, type);
        require(result.normalDampingCoefficient_ > 0.0 && result.slidingDampingCoefficient_ > 0.0 && result.rollingDampingCoefficient_ > 0.0 && result.torsionalDampingCoefficient_ > 0.0,
                "One damping switch must restore all four restitution-derived coefficients.");
    }
}

void testNodalDistribution()
{
    materialPair pair(0, 0);
    pair.setLevelSetStiffness(200000.0, 70000.0);
    pair.setDampingCoefficients(1000.0, 500.0, 200.0, 100.0);
    const auto first = response(pair, contactType::levelSetLevelSet, 0.002, 0.1, 0.25);
    const auto second = response(pair, contactType::levelSetLevelSet, 0.002, 0.3, 0.75);
    near(first.normalStiffness_ + second.normalStiffness_, 200000.0 * 0.4, "LS stiffness must scale with nodal area.");
    near(first.normalDampingCoefficient_, 250.0, "Explicit LS normal damping is not area weighted.");
    near(second.normalDampingCoefficient_, 750.0, "Explicit LS normal damping is not area weighted.");
    near(first.slidingDampingCoefficient_ + second.slidingDampingCoefficient_, 500.0, "Surface refinement must not multiply the total explicit damping.");
    near(execution::nodalDampingWeight(1.25 * 0.25, 0.3, 0.5), 0.25, "Existing nodal mass weight must recover area fraction.");
    near(first.rollingStiffness_, 0.0, "LS nodal contacts must not add a separate rolling spring.");
    near(first.torsionalDampingCoefficient_, 0.0, "LS nodal contacts must not add a separate torsional damper.");
    pair.setLinearStiffness(1.0, 2.0, 3.0, 4.0);
    const auto linearOn = response(pair, contactType::levelSetLevelSet, 0.002, 0.1, 0.25);
    pair.setLinearStiffnessEnabled(false);
    const auto linearOff = response(pair, contactType::levelSetLevelSet, 0.002, 0.1, 0.25);
    near(linearOn.normalStiffness_, first.normalStiffness_, "The sphere linear switch must not alter LS per-area stiffness.");
    near(linearOff.slidingStiffness_, first.slidingStiffness_, "Disabling sphere linear stiffness must not disable LS stiffness.");
    pair.setDampingEnabled(false);
    const auto automatic = response(pair, contactType::levelSetLevelSet, 0.002, 0.1, 0.25);
    require(automatic.normalDampingCoefficient_ > 0.0 && automatic.slidingDampingCoefficient_ > 0.0, "The damping group must also restore restitution damping for LS nodes.");
    materialPair missing(0, 0);
    missing.setNormalStiffnessPerUnitArea(1000.0);
    missing.setDampingCoefficients(10.0, 20.0, 30.0, 40.0);
    const auto incomplete = response(missing, contactType::levelSetLevelSet);
    near(incomplete.normalStiffness_, 0.0, "An incomplete LS pair must not produce partial stiffness response.");
    near(incomplete.normalDampingCoefficient_, 0.0, "An incomplete LS pair must not produce partial dashpot response.");
}

void testDampingWithoutSprings()
{
    math::Vec3 sliding{1.0, 2.0, 3.0};
    const auto slidingForce = execution::integrateTangentialSpringForce(sliding, {2.0, 0.0, 0.0}, math::Vec3::unitZ(), 100.0, 1.0, 0.0, 3.0, 0.1);
    near(slidingForce.x, -6.0, "An active dashpot must work independently of a zero sliding stiffness.");
    require(sliding == math::Vec3::zero(), "A disabled spring must clear its elastic history.");
    math::Vec3 torsional{1.0, 2.0, 3.0};
    const auto torsionalForce = execution::integrateTorsionalSpringForce(torsional, {4.0, 0.0, 2.0}, math::Vec3::unitZ(), 100.0, 1.0, 0.0, 3.0, 0.1);
    near(torsionalForce.z, -6.0, "An active torsional dashpot must work without its spring.");
    near(torsionalForce.x, 0.0, "Torsional damping must stay on the normal axis.");
    require(torsional == math::Vec3::zero(), "A disabled torsional spring retained elastic history.");
}

void testSignedNormalResponse()
{
    using math::Vec3;
    Real magnitude = 0.0;
    Vec3 sliding = Vec3::zero();
    const auto separating = execution::calculateLinearContactForce(magnitude, sliding, {100.0, 0.0, 10.0}, Vec3::unitZ(), 0.1, 0.1, 100.0, 50.0, 0.25, 10.0, 0.0);
    near(separating.z, -90.0, "A positive-overlap separating contact must retain its signed normal dashpot response.");
    near(magnitude, 90.0, "Normal-force magnitude must remain the absolute value of the signed response.");
    near(math::norm(Vec3{separating.x, separating.y, 0.0}), 0.25 * magnitude, "Sliding Coulomb return must use the absolute normal load.");
    require(sliding.x > 0.0, "Signed normal force must not disable tangential spring history.");
    for (Real penetration : {0.0, -0.1})
    {
        magnitude = 123.0;
        sliding = {1.0, 2.0, 3.0};
        const auto absent = execution::calculateLinearContactForce(magnitude, sliding, {100.0, 0.0, 10.0}, Vec3::unitZ(), penetration, 0.1, 100.0, 50.0, 0.25, 10.0, 0.0);
        require(absent == Vec3::zero() && sliding == Vec3::zero() && magnitude == 0.0, "Nonpositive overlap must still return zero and clear sliding history.");
    }
    sliding = Vec3::zero();
    const auto compressing = execution::calculateLinearContactForce(magnitude, sliding, {0.0, 0.0, -2.0}, Vec3::unitZ(), 0.1, 0.1, 100.0, 50.0, 0.25, 10.0, 0.0);
    near(compressing.z, 30.0, "Disabling the tensile-force clamp changed a positive compressive response.");
    near(magnitude, 30.0, "A positive normal force changed its stored magnitude.");

    materialPair pair(0, 0);
    pair.setDampingCoefficients(1000.0, 0.0, 0.0, 0.0);
    for (contactType type : {contactType::sphereSphere, contactType::sphereLevelSet})
    {
        const auto hertz = response(pair, type);
        const Real elasticForce = hertz.normalStiffness_ * 0.002;
        const Real separatingSpeed = (2.0 * elasticForce + 1.0) / hertz.normalDampingCoefficient_;
        sliding = Vec3::zero();
        const auto force = execution::calculateLinearContactForce(magnitude,
                                                                  sliding,
                                                                  {0.0, 0.0, separatingSpeed},
                                                                  Vec3::unitZ(),
                                                                  0.002,
                                                                  0.1,
                                                                  hertz.normalStiffness_,
                                                                  hertz.slidingStiffness_,
                                                                  0.0,
                                                                  hertz.normalDampingCoefficient_,
                                                                  hertz.slidingDampingCoefficient_);
        near(force.z, -elasticForce - 1.0, "Resolved Hertz-Mindlin contacts must use the same signed normal response as explicit linear contacts.");
        near(magnitude, elasticForce + 1.0, "Resolved Hertz-Mindlin normal magnitude lost its absolute-value contract.");
    }
}

void testMissingLSConfiguration()
{
    LSDEM simulation;
    const int materialIndex = simulation.addMaterial(material{1.0e7, 0.25, 0.8, 1000.0});
    levelset::Sphere sphere(0.5);
    sphere.buildSurfaceNode(0);
    sphere.buildLSGrid(0.25, 2);
    const int geometryIndex = simulation.addGeometry(sphere);
    for (int i = 0; i < 2; ++i)
    {
        LSParticle particle;
        particle.setMaterial(simulation.materials(), materialIndex);
        particle.setGeometry(simulation.geometries(), geometryIndex);
        particle.setPosition({3.0 * i, 0.0, 0.0});
        simulation.addLSParticle(particle);
    }
    bool rejected = false;
    try
    {
        simulation.initialize();
    }
    catch (const std::invalid_argument& error)
    {
        rejected = std::string(error.what()).find("LS-LS") != std::string::npos;
    }
    require(rejected, "Initialization must report missing same-material LS stiffness before force evaluation.");
    materialPair pair(materialIndex, materialIndex);
    pair.setLevelSetStiffness(0.0, 0.0);
    simulation.setMaterialPair(pair);
    simulation.initialize(); // Explicit zero LS stiffness is configured, not missing.
}

void testMaterialClone()
{
    LSDEM simulation;
    simulation.addMaterial(material{1.0e7, 0.25, 0.8, 1000.0});
    for (int i = 0; i < 20; ++i)
        simulation.addMaterial(simulation.materials().host()[0]);
    require(simulation.materials().hostSize() == 21 && simulation.materialPairs().hostSize() == 231, "Cloning existing materials must append a complete pair table.");
    for (const material& value : simulation.materials().host())
        near(value.youngModulus(), 1.0e7, "Appending an aliased material must survive vector reallocation.");
}
} // namespace

int main()
{
    try
    {
        testHertzMindlin();
        testGroupActivation();
        testNodalDistribution();
        testDampingWithoutSprings();
        testSignedNormalResponse();
        testMissingLSConfiguration();
        testMaterialClone();
        std::cout << "PASS Hertz-Mindlin, explicit pair overrides, activation, nodal damping and missing LS configuration\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
