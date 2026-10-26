/**
 * @file ContactMaterialTest.cu
 * @brief Checks material-pair response and particle scattering against CPU and analytical results.
 */
#include "TestSupport.h"
#include "execution/cuda/contactForceKernel.cuh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace
{
using namespace fundem;
using math::Real;
using math::Vec3;

constexpr Real timeStep = 1.0e-4;
constexpr Real overlap = 0.02;

class testStream
{
public:
    testStream() { FUNDEM_TEST_REQUIRE_CUDA(cudaStreamCreateWithFlags(&value_, cudaStreamNonBlocking)); }
    ~testStream() { cudaStreamDestroy(value_); }
    operator cudaStream_t() const noexcept { return value_; }

private:
    cudaStream_t value_{nullptr};
};

void requireNear(Real actual, Real expected, const char* message)
{
    if (!math::isFinite(actual) || std::abs(actual - expected) > 1.0e-10 * (1.0 + std::abs(expected)))
        throw std::runtime_error(message);
}

void requireNear(const Vec3& actual, const Vec3& expected, const char* message)
{
    if (!math::isFinite(actual) || math::norm(actual - expected) > 1.0e-10 * (1.0 + math::norm(expected)))
        throw std::runtime_error(message);
}

void fillMaterials(materialContainer& materials)
{
    materials.host().push_back(material{1.0e7, 0.22, 0.6, 1000.0});
    materials.host().push_back(material{1.7e7, 0.31, 0.8, 1200.0});
}

template <class Particle> Particle makeParticle(const materialContainer& materials, int materialIndex, const Vec3& position, bool infiniteMass = false)
{
    Particle value;
    value.setPosition(position);
    if constexpr (std::is_same_v<Particle, particle>)
        value.setRadius(0.1);
    value.setMaterial(materials, materialIndex);
    if constexpr (std::is_base_of_v<LSParticle, Particle>)
        static_cast<rigidBody&>(value).setMassAndInertia(materialIndex == 0 ? 2.0 : 3.0, math::Mat3::identity());
    if (infiniteMass)
        value.setInfiniteMass();
    return value;
}

template <class MasterStorage, class SlaveStorage> contact makeContact(const MasterStorage& masters, int masterIndex, const SlaveStorage& slaves, int slaveIndex, Real effectiveRadius)
{
    contact value;
    FUNDEM_TEST_REQUIRE(value.setMasterSlaveParticle(masters, masterIndex, slaves, slaveIndex));
    value.setPoint({0.01, 0.0, 0.0});
    value.setNormal(Vec3::unitX());
    value.setOverlap(overlap);
    value.setArea(3.0); // Sphere/LS response must not multiply stiffness by nodal area.
    const Real inverseMassSum = masters.host()[masterIndex].inverseMass() + slaves.host()[slaveIndex].inverseMass();
    value.setEffectiveMass(1.0 / inverseMassSum);
    value.setEffectiveRadius(effectiveRadius);
    value.setSlidingSpringDeformation({0.0, 1.0e-4, 2.0e-4});
    value.setRollingSpringDeformation({0.0, 0.0, 1.0e-4});
    value.setTorsionalSpringDeformation({2.0e-4, 0.0, 0.0});
    return value;
}

void compareContacts(const contact& actual, const contact& expected)
{
    requireNear(actual.force(), expected.force(), "CPU/CUDA contact forces differ.");
    requireNear(actual.torque(), expected.torque(), "CPU/CUDA contact torques differ.");
    requireNear(actual.normalForceMagnitude(), expected.normalForceMagnitude(), "CPU/CUDA normal forces differ.");
    requireNear(actual.slidingSpringDeformation(), expected.slidingSpringDeformation(), "CPU/CUDA sliding histories differ.");
    requireNear(actual.rollingSpringDeformation(), expected.rollingSpringDeformation(), "CPU/CUDA rolling histories differ.");
    requireNear(actual.torsionalSpringDeformation(), expected.torsionalSpringDeformation(), "CPU/CUDA torsional histories differ.");
    requireNear(actual.normalElasticEnergy(), expected.normalElasticEnergy(), "CPU/CUDA normal energies differ.");
    requireNear(actual.slidingElasticEnergy(), expected.slidingElasticEnergy(), "CPU/CUDA sliding energies differ.");
    requireNear(actual.rollingElasticEnergy(), expected.rollingElasticEnergy(), "CPU/CUDA rolling energies differ.");
    requireNear(actual.torsionalElasticEnergy(), expected.torsionalElasticEnergy(), "CPU/CUDA torsional energies differ.");
}

template <class MasterStorage, class SlaveStorage>
std::vector<contact> evaluate(MasterStorage& masters, SlaveStorage& slaves, materialContainer& materials, materialPairContainer& pairs, const std::vector<contact>& initial)
{
    bool sameContainer = false;
    if constexpr (std::is_same_v<MasterStorage, SlaveStorage>)
        sameContainer = &masters == &slaves;
    constexpr auto type = execution::contactTypeFor<typename MasterStorage::host_type, typename SlaveStorage::host_type>();
    contactContainer cpuContacts, gpuContacts;
    for (const auto& value : initial)
    {
        cpuContacts.host().push_back(value);
        gpuContacts.host().push_back(value);
    }
    std::vector<Vec3> masterForces(masters.hostSize(), Vec3::zero()), masterTorques(masters.hostSize(), Vec3::zero());
    std::vector<Vec3> slaveForces(slaves.hostSize(), Vec3::zero()), slaveTorques(slaves.hostSize(), Vec3::zero());
    for (auto& value : cpuContacts.host())
    {
        FUNDEM_TEST_REQUIRE(value.calculateForce(timeStep, pairs, type));
        const int master = value.masterParticleIndex(), slave = value.slaveParticleIndex();
        masterForces[master] += value.force();
        masterTorques[master] += math::cross(value.point() - masters.host()[master].position(), value.force()) + value.torque();
        auto& forces = sameContainer ? masterForces : slaveForces;
        auto& torques = sameContainer ? masterTorques : slaveTorques;
        forces[slave] -= value.force();
        torques[slave] += math::cross(value.point() - slaves.host()[slave].position(), -value.force()) - value.torque();
    }
    testStream stream;
    materials.copyHostToDeviceAsync(stream);
    pairs.copyHostToDeviceAsync(stream);
    masters.copyHostToDeviceAsync(stream);
    if (!sameContainer)
        slaves.copyHostToDeviceAsync(stream);
    gpuContacts.copyHostToDeviceAsync(stream);
    if constexpr (std::is_same_v<MasterStorage, SlaveStorage>)
    {
        if (sameContainer)
            cuda::launchContactForce(gpuContacts, masters, materials, pairs, timeStep, stream);
        else
            cuda::launchContactForce(gpuContacts, masters, slaves, materials, pairs, timeStep, stream);
    }
    else
        cuda::launchContactForce(gpuContacts, masters, slaves, materials, pairs, timeStep, stream);
    gpuContacts.copyDeviceToHostAsync(stream);
    masters.copyDeviceToHostAsync(stream);
    if (!sameContainer)
        slaves.copyDeviceToHostAsync(stream);
    FUNDEM_TEST_REQUIRE_CUDA(cudaStreamSynchronize(stream));
    for (int i = 0; i < static_cast<int>(initial.size()); ++i)
        compareContacts(gpuContacts.host()[i], cpuContacts.host()[i]);
    for (int i = 0; i < static_cast<int>(masters.hostSize()); ++i)
    {
        requireNear(masters.host()[i].force(), masterForces[i], "CUDA master force scattering is incorrect.");
        requireNear(masters.host()[i].torque(), masterTorques[i], "CUDA master torque scattering is incorrect.");
    }
    if (!sameContainer)
        for (int i = 0; i < static_cast<int>(slaves.hostSize()); ++i)
        {
            requireNear(slaves.host()[i].force(), slaveForces[i], "CUDA slave reaction scattering is incorrect.");
            requireNear(slaves.host()[i].torque(), slaveTorques[i], "CUDA slave reaction torque scattering is incorrect.");
        }
    return {gpuContacts.host().begin(), gpuContacts.host().end()};
}

Real dissipationFactor(const particle& master, const particle& slave)
{
    const Real first = master.getMaterial().restitutionCoefficient();
    const Real second = slave.getMaterial().restitutionCoefficient();
    const Real logarithm = std::log(2.0 * first * second / (first + second));
    return -logarithm / std::sqrt(logarithm * logarithm + math::pi * math::pi);
}

void checkOrdinaryResponse(const contact& result, const contact& initial, const particle& master, const particle& slave, const materialPair& pair)
{
    const auto& first = master.getMaterial();
    const auto& second = slave.getMaterial();
    const Real effectiveE = 1.0 / ((1.0 - first.poissonRatio() * first.poissonRatio()) / first.youngModulus() + (1.0 - second.poissonRatio() * second.poissonRatio()) / second.youngModulus());
    const Real firstG = first.youngModulus() / (2.0 * (1.0 + first.poissonRatio()));
    const Real secondG = second.youngModulus() / (2.0 * (1.0 + second.poissonRatio()));
    const Real effectiveG = 1.0 / ((2.0 - first.poissonRatio()) / firstG + (2.0 - second.poissonRatio()) / secondG);
    const Real radius = std::sqrt(initial.effectiveRadius() * overlap);
    const bool linear = pair.linearStiffnessEnabled();
    const Real kn = linear ? pair.normalStiffness() : (4.0 / 3.0) * effectiveE * radius;
    const Real ks = linear ? pair.slidingStiffness() : 8.0 * effectiveG * radius;
    const Real kr = linear ? pair.rollingStiffness() : ks;
    const Real kt = linear ? pair.torsionalStiffness() : ks;
    const Real beta = dissipationFactor(master, slave);
    const Real dampingScale = linear ? 1.0 : std::sqrt(5.0 / 6.0);
    const Real cn = pair.dampingEnabled() ? pair.normalDampingCoefficient() : 2.0 * dampingScale * beta * std::sqrt(initial.effectiveMass() * (linear ? 1.0 : 1.5) * kn);
    const Real cs = pair.dampingEnabled() ? pair.slidingDampingCoefficient() : 2.0 * dampingScale * beta * std::sqrt(initial.effectiveMass() * ks);
    const Real cr = pair.dampingEnabled() ? pair.rollingDampingCoefficient() : 2.0 * beta * std::sqrt(initial.effectiveMass() * kr);
    const Real ct = pair.dampingEnabled() ? pair.torsionalDampingCoefficient() : 2.0 * beta * std::sqrt(initial.effectiveMass() * kt);
    const Vec3 relativeVelocity =
        execution::relativeVelocityAtContactPoint(master.position(), master.velocity(), master.angularVelocity(), slave.position(), slave.velocity(), slave.angularVelocity(), initial.point());
    const Real vn = math::dot(relativeVelocity, initial.normal());
    const Vec3 vt = relativeVelocity - vn * initial.normal();
    const Real signedNormalForce = kn * overlap - cn * vn;
    const Real fn = std::abs(signedNormalForce);
    const Vec3 slidingSpring = ks > 0.0 && fn > 0.0 ? initial.slidingSpringDeformation() + timeStep * vt : Vec3::zero();
    const Vec3 slidingForce = fn > 0.0 ? -ks * slidingSpring - cs * vt : Vec3::zero();
    const Vec3 relativeOmega = master.angularVelocity() - slave.angularVelocity();
    const Real momentRadius = initial.effectiveRadius();
    const Vec3 rollingVelocity = -momentRadius * math::cross(initial.normal(), relativeOmega);
    const Vec3 torsionalVelocity = 2.0 * momentRadius * math::dot(initial.normal(), relativeOmega) * initial.normal();
    const Vec3 rollingSpring = kr > 0.0 && fn > 0.0 ? initial.rollingSpringDeformation() + timeStep * rollingVelocity : Vec3::zero();
    const Vec3 torsionalSpring = kt > 0.0 && fn > 0.0 ? initial.torsionalSpringDeformation() + timeStep * torsionalVelocity : Vec3::zero();
    const Vec3 rollingForce = fn > 0.0 ? -kr * rollingSpring - cr * rollingVelocity : Vec3::zero();
    const Vec3 torsionalForce = fn > 0.0 ? -kt * torsionalSpring - ct * torsionalVelocity : Vec3::zero();
    FUNDEM_TEST_REQUIRE(math::norm(slidingForce) <= 10.0 * fn && math::norm(rollingForce) <= 10.0 * fn && math::norm(torsionalForce) <= 10.0 * fn);
    requireNear(result.force(), signedNormalForce * initial.normal() + slidingForce, "Grouped stiffness/damping or Hertz-Mindlin force is incorrect.");
    requireNear(result.normalForceMagnitude(), fn, "The recorded normal force must remain an absolute magnitude.");
    requireNear(result.torque(), momentRadius * math::cross(initial.normal(), rollingForce) + 2.0 * momentRadius * torsionalForce, "Grouped rotational response is incorrect.");
    requireNear(result.normalElasticEnergy(), (linear ? 0.5 : 0.4) * kn * overlap * overlap, "The selected normal potential is incorrect.");
    requireNear(result.slidingSpringDeformation(), slidingSpring, "The selected sliding spring state is incorrect.");
    requireNear(result.rollingSpringDeformation(), rollingSpring, "The selected rolling spring state is incorrect.");
    requireNear(result.torsionalSpringDeformation(), torsionalSpring, "The selected torsional spring state is incorrect.");
    requireNear(result.slidingElasticEnergy(), 0.5 * ks * math::normSquared(slidingSpring), "The selected sliding energy is incorrect.");
    requireNear(result.rollingElasticEnergy(), 0.5 * kr * math::normSquared(rollingSpring), "The selected rolling energy is incorrect.");
    requireNear(result.torsionalElasticEnergy(), 0.5 * kt * math::normSquared(torsionalSpring), "The selected torsional energy is incorrect.");
}

enum class pairMode
{
    hertz,
    explicitLinear,
    activeZero,
    linearDisabled,
    dampingDisabled,
    explicitDamping,
    zeroDamping,
    damperOnly
};

template <class Slave> void testSpherePair(pairMode mode, bool infiniteMass = false, bool separating = false)
{
    materialContainer materials;
    fillMaterials(materials);
    materialPairContainer pairs(2);
    auto& pair = pairs.at(0, 1);
    FUNDEM_TEST_REQUIRE(!pair.linearStiffnessEnabled() && !pair.dampingEnabled());
    pair.setFrictionCoefficients(10.0, 10.0, 10.0);
    if (mode != pairMode::hertz)
    {
        const bool zeroStiffness = mode == pairMode::activeZero;
        const bool zeroDamping = zeroStiffness || mode == pairMode::zeroDamping;
        const bool zeroTangentialStiffness = zeroStiffness || mode == pairMode::damperOnly;
        pair.setLinearStiffness(zeroStiffness ? 0.0 : 1000.0, zeroTangentialStiffness ? 0.0 : 500.0, zeroTangentialStiffness ? 0.0 : 200.0, zeroTangentialStiffness ? 0.0 : 100.0);
        pair.setDampingCoefficients(zeroDamping ? 0.0 : 13.0, zeroDamping ? 0.0 : 7.0, zeroDamping ? 0.0 : 3.0, zeroDamping ? 0.0 : 2.0);
        pair.setLinearStiffnessEnabled(true);
        pair.setDampingEnabled(true);
        if (mode == pairMode::linearDisabled || mode == pairMode::explicitDamping || mode == pairMode::zeroDamping)
            pair.setLinearStiffnessEnabled(false);
        if (mode == pairMode::linearDisabled || mode == pairMode::dampingDisabled)
            pair.setDampingEnabled(false);
    }
    particleContainer masters;
    hostAoSDeviceSoA<Slave, typename Slave::DeviceLayout> slaves;
    masters.host().push_back(makeParticle<particle>(materials, 0, {0.1, 0.0, 0.0}));
    slaves.host().push_back(makeParticle<Slave>(materials, 1, {-0.08, 0.0, 0.0}, infiniteMass));
    masters.host()[0].setVelocity({separating ? 1.0e6 : 0.2, 0.3, -0.1});
    masters.host()[0].setAngularVelocity({0.4, 0.5, 0.6});
    slaves.host()[0].setVelocity({-0.1, 0.05, 0.07});
    slaves.host()[0].setAngularVelocity({-0.2, -0.1, 0.15});
    const Real radius = std::is_same_v<Slave, particle> ? 0.05 : 0.1;
    const auto initial = makeContact(masters, 0, slaves, 0, radius);
    const auto results = evaluate(masters, slaves, materials, pairs, {initial});
    checkOrdinaryResponse(results[0], initial, masters.host()[0], slaves.host()[0], pair);
    if (separating)
        FUNDEM_TEST_REQUIRE(math::dot(results[0].force(), initial.normal()) < 0.0 && results[0].normalForceMagnitude() > 0.0);
}

void testNodalContact(bool infiniteMass, bool linearEnabled, bool dampingEnabled, bool separating = false)
{
    materialContainer materials;
    fillMaterials(materials);
    materialPairContainer pairs(2);
    auto& pair = pairs.at(0, 1);
    pair.setLevelSetStiffness(2400.0, 1200.0);
    pair.setLinearStiffness(9000.0, 8000.0, 7000.0, 6000.0); // LS nodes must ignore direct sphere stiffnesses.
    pair.setDampingCoefficients(12.0, 6.0, 3.0, 2.0);
    pair.setLinearStiffnessEnabled(linearEnabled);
    pair.setDampingEnabled(dampingEnabled);
    pair.setFrictionCoefficients(10.0, 10.0, 10.0);
    LSParticleContainer particles;
    particles.host().push_back(makeParticle<LSParticle>(materials, 1, {0.1, 0.0, 0.0}));
    particles.host().push_back(makeParticle<LSParticle>(materials, 0, {-0.08, 0.0, 0.0}, infiniteMass));
    particles.host()[0].setVelocity({separating ? 1.0e6 : 0.4, 0.2, 0.0});
    particles.host()[1].setVelocity({-0.1, -0.05, 0.0});
    const auto seed = makeContact(particles, 0, particles, 1, 0.0);
    std::vector<contact> initial;
    for (Real fraction : {0.25, 0.75})
    {
        auto value = seed;
        value.setArea(fraction);
        value.setEffectiveMass(seed.effectiveMass() * fraction);
        value.setSlidingSpringDeformation({0.0, 1.0e-4, 0.0});
        initial.push_back(value);
    }
    const auto results = evaluate(particles, particles, materials, pairs, initial);
    const Real beta = dissipationFactor(particles.host()[0], particles.host()[1]);
    const Real cn = dampingEnabled ? 12.0 : 2.0 * beta * std::sqrt(seed.effectiveMass() * 2400.0);
    const Real cs = dampingEnabled ? 6.0 : 2.0 * beta * std::sqrt(seed.effectiveMass() * 1200.0);
    const Real normalSpeed = particles.host()[0].velocity().x - particles.host()[1].velocity().x;
    for (int i = 0; i < 2; ++i)
    {
        const Real fraction = i == 0 ? 0.25 : 0.75;
        const Real spring = 1.0e-4 + timeStep * 0.25;
        requireNear(results[i].force(),
                    fraction * Vec3{2400.0 * overlap - cn * normalSpeed, -1200.0 * spring - cs * 0.25, 0.0},
                    "LS damping must scale with the nodal mass fraction and ignore the ordinary stiffness switch.");
        requireNear(results[i].normalForceMagnitude(), fraction * std::abs(2400.0 * overlap - cn * normalSpeed), "LS normal output must retain the absolute signed-response magnitude.");
        if (separating)
            FUNDEM_TEST_REQUIRE(results[i].force().x < 0.0 && results[i].normalForceMagnitude() > 0.0);
        requireNear(results[i].normalElasticEnergy(), 0.5 * 2400.0 * fraction * overlap * overlap, "LS normal stiffness must use node area.");
        requireNear(results[i].torque(), Vec3::zero(), "LS nodes must not apply direct rotational resistance.");
        requireNear(results[i].rollingSpringDeformation(), Vec3::zero(), "LS rolling history must be cleared.");
        requireNear(results[i].torsionalSpringDeformation(), Vec3::zero(), "LS torsional history must be cleared.");
        requireNear(results[i].rollingElasticEnergy(), 0.0, "LS rolling energy must be zero.");
        requireNear(results[i].torsionalElasticEnergy(), 0.0, "LS torsional energy must be zero.");
    }
}
} // namespace

int main()
{
    try
    {
        for (const auto mode : {pairMode::hertz,
                                pairMode::explicitLinear,
                                pairMode::activeZero,
                                pairMode::linearDisabled,
                                pairMode::dampingDisabled,
                                pairMode::explicitDamping,
                                pairMode::zeroDamping,
                                pairMode::damperOnly})
        {
            testSpherePair<particle>(mode);
            testSpherePair<LSParticle>(mode);
        }
        testSpherePair<LSParticle>(pairMode::hertz, true);
        testSpherePair<LSParticle>(pairMode::explicitLinear, true);
        for (const auto mode : {pairMode::hertz, pairMode::explicitLinear})
        {
            testSpherePair<particle>(mode, false, true);
            testSpherePair<LSParticle>(mode, false, true);
        }
        for (bool infiniteMass : {false, true})
            for (bool linearEnabled : {false, true})
                for (bool dampingEnabled : {false, true})
                    testNodalContact(infiniteMass, linearEnabled, dampingEnabled);
        testNodalContact(false, false, true, true);
        std::cout << "PASS CPU/CUDA two-switch pair laws, Hertz-Mindlin defaults, LS nodal weights, histories, energies and streamed scattering\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
