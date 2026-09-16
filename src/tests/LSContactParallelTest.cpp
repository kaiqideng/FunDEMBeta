#include "data/myLSObject.h"
#include "execution/contactDetection.h"
#include "execution/cpu/particleFunctions.h"
#include "interaction/contactSearch.h"
#include "material/LSMaterial.h"
#include "solver/LSDEM.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace
{

using namespace fundem;
using math::Real;
using math::Vec3;
using Key = std::tuple<int, int, int>;

void require(bool condition, const std::string& message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

bool close(Real first, Real second, Real tolerance = 2.0e-10)
{
    return std::isfinite(first) && std::isfinite(second) && std::abs(first - second) <= tolerance * (1.0 + std::max(std::abs(first), std::abs(second)));
}

bool close(const Vec3& first, const Vec3& second, Real tolerance = 2.0e-10)
{
    return close(first.x, second.x, tolerance) && close(first.y, second.y, tolerance) && close(first.z, second.z, tolerance);
}

void setThreads(int count)
{
#if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(count);
#else
    (void)count;
#endif
}

Key key(const contact& value) { return {value.masterParticleIndex(), value.slaveParticleIndex(), value.particleSurfaceNodeIndex()}; }

struct History {
    Vec3 sliding;
    Vec3 rolling;
    Vec3 torsional;
};

using Histories = std::map<Key, History>;

Histories readHistories(const contactContainer& contacts)
{
    Histories result;
    for (const contact& value : contacts.host())
    {
        require(result.emplace(key(value), History{value.slidingSpringDeformation(), value.rollingSpringDeformation(), value.torsionalSpringDeformation()}).second,
                "Contact keys must be unique.");
    }
    return result;
}

void seedHistories(contactContainer& contacts)
{
    for (contact& value : contacts.host())
    {
        const Real seed = 1.0e-7 * (1 + value.particleSurfaceNodeIndex() + 31 * value.masterParticleIndex() + 173 * value.slaveParticleIndex());
        value.setSlidingSpringDeformation({seed, -0.7 * seed, 0.3 * seed});
        value.setRollingSpringDeformation({-0.2 * seed, seed, 0.4 * seed});
        value.setTorsionalSpringDeformation({0.6 * seed, 0.1 * seed, -seed});
    }
}

struct Fixture {
    LSDEM definitions;
    LSParticleContainer particles;

    Fixture() { definitions.addMaterial(LSMaterial{2.0e5, 7.0e4, 0.6, 0.45, 1200.0}); }

    int sphere(int subdivision, Real radius)
    {
        levelset::Sphere geometry(radius);
        geometry.buildSurfaceNode(subdivision);
        geometry.buildLSGrid(20);
        return definitions.addGeometry(geometry);
    }

    int ellipsoid(int subdivision)
    {
        levelset::Superellipsoid geometry(0.46, 0.38, 0.52, 1.0, 1.0);
        geometry.buildSurfaceNode(subdivision);
        geometry.buildLSGrid(20);
        return definitions.addGeometry(geometry);
    }

    void add(int geometry, const Vec3& position, Real angle)
    {
        LSParticle value;
        value.setMaterial(definitions.materials(), 0);
        value.setGeometry(definitions.geometries(), geometry);
        value.setPosition(position);
        value.setOrientation(math::Quaternion::fromUnitAxisAngle(math::normalizedOrZero(Vec3{1.0, 2.0, -1.0}), angle));
        value.setVelocity({0.03 + 0.01 * position.x, -0.02, 0.01});
        value.setAngularVelocity({-0.02, 0.04, 0.01});
        particles.host().push_back(std::move(value));
    }

    void wall()
    {
        levelset::PlaneWall geometry(Vec3::unitZ(), 5.0);
        geometry.buildLSGrid(0.15, 3);
        const int index = definitions.addGeometry(geometry, true);
        add(index, {0.13, -0.17, 0.025}, 0.07);
        require(particles.host().back().inverseMass() == 0.0, "Wall fixture must have infinite mass.");
    }

    void move()
    {
        for (std::size_t index = 0; index < particles.hostSize(); ++index)
        {
            auto& value = particles.host()[index];
            const Real scale = static_cast<Real>(index + 1);
            value.setPosition(value.position() + Vec3{0.008 * scale, -0.003 * scale, 0.002 * scale});
            value.setOrientation(math::Quaternion::fromUnitAxisAngle(math::normalizedOrZero(Vec3{2.0, -1.0, 3.0}), 0.12 + 0.09 * scale));
        }
    }
};

// Deliberately retain the world-coordinate, per-node public query path as an
// independent oracle for cached relative transforms and flattened node tasks.
contactContainer referenceContacts(const LSParticleContainer& particles, const std::vector<particlePair>& pairs, const Histories& history)
{
    contactContainer result;
    auto& output = result.host();
    for (const particlePair& pair : pairs)
    {
        const LSParticle& master = particles.host()[pair.masterParticleIndex_];
        const LSParticle& slave = particles.host()[pair.slaveParticleIndex_];
        const std::size_t begin = output.size();
        for (std::size_t nodeIndex = 0; nodeIndex < master.surfaceNodes().size(); ++nodeIndex)
        {
            const auto& node = master.surfaceNodes()[nodeIndex];
            const Vec3 position = master.position() + math::rotateUnit(master.orientation(), node.position_);
            if (math::distanceSquared(position, slave.position()) > slave.boundingRadius() * slave.boundingRadius())
            {
                continue;
            }
            Real overlap = 0.0;
            Vec3 normal;
            if (!execution::detectLevelSetContact(overlap, normal, slave.gridNodes().data(), position, slave.position(), slave.orientation(), slave.gridNodeOrigin(),
                                                  slave.gridNodeInverseSpacing(), slave.gridNodeSize(), 0))
            {
                continue;
            }
            contact value;
            require(value.setMasterSlaveParticle(particles, pair.masterParticleIndex_, pair.slaveParticleIndex_), "Reference particle snapshots must be valid.");
            value.setParticleSurfaceNodeIndex(static_cast<int>(nodeIndex));
            value.setPoint(execution::surfaceNodeLevelSetContactPoint(position, normal, overlap));
            value.setNormal(normal);
            value.setOverlap(overlap);
            value.setArea(node.area_);
            value.setEffectiveRadius(0.0);
            const auto previous = history.find(key(value));
            if (previous != history.end())
            {
                value.setSlidingSpringDeformation(previous->second.sliding);
                value.setRollingSpringDeformation(previous->second.rolling);
                value.setTorsionalSpringDeformation(previous->second.torsional);
            }
            output.push_back(std::move(value));
        }
        Real area = 0.0;
        for (std::size_t index = begin; index < output.size(); ++index)
        {
            area += output[index].area();
        }
        const Real massPerArea = area > 0.0 ? execution::effectiveMass(master.inverseMass(), slave.inverseMass()) / area : 0.0;
        for (std::size_t index = begin; index < output.size(); ++index)
        {
            output[index].setEffectiveMass(output[index].area() * massPerArea);
        }
    }
    return result;
}

void compareContacts(const contactContainer& actual, const contactContainer& expected, const std::string& context, bool forces = false)
{
    require(actual.hostSize() == expected.hostSize(), context + ": contact count differs.");
    for (std::size_t index = 0; index < actual.hostSize(); ++index)
    {
        const auto& a = actual.host()[index];
        const auto& e = expected.host()[index];
        require(key(a) == key(e), context + ": deterministic pair/node order differs at " + std::to_string(index));
        require(close(a.point(), e.point()) && close(a.normal(), e.normal()) && close(a.overlap(), e.overlap()), context + ": contact geometry differs.");
        require(close(a.area(), e.area()) && close(a.effectiveMass(), e.effectiveMass()), context + ": area or nodal effective mass differs.");
        require(close(a.slidingSpringDeformation(), e.slidingSpringDeformation()) && close(a.rollingSpringDeformation(), e.rollingSpringDeformation()) &&
                    close(a.torsionalSpringDeformation(), e.torsionalSpringDeformation()), context + ": contact history differs.");
        if (forces)
        {
            require(close(a.force(), e.force(), 2.0e-8) && close(a.torque(), e.torque(), 2.0e-8) && close(a.normalElasticEnergy(), e.normalElasticEnergy(), 2.0e-8) &&
                        close(a.slidingElasticEnergy(), e.slidingElasticEnergy(), 2.0e-8), context + ": contact response differs.");
        }
    }
}

void checkPairMasses(const contactContainer& contacts, const LSParticleContainer& particles)
{
    struct Sum { Real mass{0.0}; Real area{0.0}; int lastNode{-1}; };
    std::map<std::pair<int, int>, Sum> sums;
    bool wallContact = false;
    for (const auto& value : contacts.host())
    {
        const auto& master = particles.host()[value.masterParticleIndex()];
        const auto& slave = particles.host()[value.slaveParticleIndex()];
        require(master.inverseMass() > 0.0, "Infinite-mass wall must never be the master.");
        wallContact = wallContact || slave.inverseMass() == 0.0;
        require(value.area() > 0.0 && value.effectiveMass() > 0.0, "Active nodal area and effective mass must be positive.");
        auto& sum = sums[{value.masterParticleIndex(), value.slaveParticleIndex()}];
        require(value.particleSurfaceNodeIndex() > sum.lastNode, "Surface-node order must increase across chunk boundaries.");
        sum.lastNode = value.particleSurfaceNodeIndex();
        sum.mass += value.effectiveMass();
        sum.area += value.area();
    }
    require(wallContact, "Fixture must exercise finite/infinite-mass contact.");
    for (const auto& entry : sums)
    {
        const auto& master = particles.host()[entry.first.first];
        const auto& slave = particles.host()[entry.first.second];
        const Real mass = execution::effectiveMass(master.inverseMass(), slave.inverseMass());
        require(close(entry.second.mass, mass), "Nodal masses must sum to the complete pair effective mass.");
    }
    for (const auto& value : contacts.host())
    {
        const auto& sum = sums.at({value.masterParticleIndex(), value.slaveParticleIndex()});
        require(close(value.effectiveMass(), sum.mass * value.area() / sum.area), "Nodal mass must be weighted by the whole active pair area, not one task chunk.");
    }
}

struct Loads { Vec3 force{Vec3::zero()}; Vec3 torque{Vec3::zero()}; };

std::vector<Loads> referenceLoads(contactContainer& contacts, const LSParticleContainer& particles)
{
    std::vector<Loads> loads(particles.hostSize());
    for (auto& value : contacts.host())
    {
        require(value.calculateForce(1.0e-5), "Reference force must consume a valid freshly refreshed contact.");
        const int master = value.masterParticleIndex();
        const int slave = value.slaveParticleIndex();
        loads[master].force += value.force();
        loads[master].torque += math::cross(value.point() - particles.host()[master].position(), value.force()) + value.torque();
        loads[slave].force -= value.force();
        loads[slave].torque += math::cross(value.point() - particles.host()[slave].position(), -value.force()) - value.torque();
    }
    return loads;
}

void compareLoads(contactContainer& contacts, LSParticleContainer& particles, const std::vector<Loads>& expected, int threads)
{
    setThreads(threads);
    cpu::clearForceAndTorque(particles);
    cpu::addContactForceAndTorque(contacts, particles, particles, 1.0e-5);
    Vec3 totalForce = Vec3::zero();
    for (std::size_t index = 0; index < particles.hostSize(); ++index)
    {
        const auto& value = particles.host()[index];
        require(close(value.force(), expected[index].force, 2.0e-8) && close(value.torque(), expected[index].torque, 2.0e-8), "Assembled particle/wall loads differ from direct reference scatter.");
        totalForce += value.force();
    }
    require(close(totalForce, Vec3::zero(), 1.0e-7), "Including wall reactions, total contact force must balance.");
}

void checkCacheLifecycle(Fixture& fixture, contactSearch& serialSearch, contactSearch& parallelSearch, contactContainer& serial, contactContainer& parallel)
{
    std::vector<Vec3> positions;
    for (std::size_t index = 0; index < fixture.particles.hostSize(); ++index)
    {
        auto& value = fixture.particles.host()[index];
        positions.push_back(value.position());
        value.setPosition({100.0 + 100.0 * static_cast<Real>(index), 100.0, 100.0});
    }
    setThreads(1);
    serialSearch.findContacts(serial, fixture.particles);
    setThreads(4);
    parallelSearch.findContacts(parallel, fixture.particles);
    require(serialSearch.candidatePairCount() == 0 && parallelSearch.candidatePairCount() == 0 && serial.hostSize() == 0 && parallel.hostSize() == 0,
            "A zero-work search must discard every previously cached contact.");

    for (std::size_t index = 0; index < fixture.particles.hostSize(); ++index)
    {
        fixture.particles.host()[index].setPosition(positions[index]);
    }
    setThreads(1);
    serialSearch.findContacts(serial, fixture.particles);
    setThreads(4);
    parallelSearch.findContacts(parallel, fixture.particles);
    auto reference = referenceContacts(fixture.particles, serialSearch.candidatePairs(), {});
    require(reference.hostSize() > 0, "Restoring particle positions must recreate contacts.");
    compareContacts(serial, reference, "serial cache regrowth with fresh history");
    compareContacts(parallel, reference, "parallel cache regrowth with fresh history");

    seedHistories(serial);
    seedHistories(parallel);
    const auto history = readHistories(serial);
    serialSearch.clear();
    parallelSearch.clear();
    require(serialSearch.candidatePairCount() == 0 && parallelSearch.candidatePairCount() == 0, "clear() must clear candidate work.");
    setThreads(1);
    serialSearch.findContacts(serial, fixture.particles);
    setThreads(4);
    parallelSearch.findContacts(parallel, fixture.particles);
    reference = referenceContacts(fixture.particles, serialSearch.candidatePairs(), history);
    compareContacts(serial, reference, "serial cleared search retains container-owned history");
    compareContacts(parallel, reference, "parallel cleared search retains container-owned history");
    checkPairMasses(parallel, fixture.particles);
}

void exercise(Fixture& fixture, const std::string& name)
{
    contactSearch serialSearch;
    contactSearch parallelSearch;
    contactContainer serial;
    contactContainer parallel;
    Histories history;
    for (int stage = 0; stage < 2; ++stage)
    {
        setThreads(1);
        serialSearch.findContacts(serial, fixture.particles);
        setThreads(4);
        parallelSearch.findContacts(parallel, fixture.particles);
        require(serialSearch.candidatePairCount() > 0 && serialSearch.candidatePairCount() < 256, "Fixture must exercise expensive nodes below the old pair threshold.");
        require(serial.hostSize() >= 256, "Fixture must also exercise staged contact force assembly.");
        auto reference = referenceContacts(fixture.particles, serialSearch.candidatePairs(), history);
        compareContacts(parallel, serial, name + " thread comparison");
        compareContacts(serial, reference, name + " world-coordinate reference");
        checkPairMasses(parallel, fixture.particles);
        if (stage == 0)
        {
            seedHistories(serial);
            seedHistories(parallel);
            history = readHistories(serial);
            fixture.move();
        }
        else
        {
            std::size_t retained = 0;
            for (const auto& value : serial.host())
            {
                retained += history.count(key(value));
            }
            require(retained > 0 && retained < serial.hostSize() && retained < history.size(), "Moving fixture must retain, create, and remove contact-history keys.");
            const auto loads = referenceLoads(reference, fixture.particles);
            compareLoads(serial, fixture.particles, loads, 1);
            compareLoads(parallel, fixture.particles, loads, 4);
            compareContacts(serial, reference, name + " serial response", true);
            compareContacts(parallel, reference, name + " parallel response", true);
        }
    }
    checkCacheLifecycle(fixture, serialSearch, parallelSearch, serial, parallel);
    std::cout << name << ": " << serialSearch.candidatePairCount() << " candidate pairs, " << serial.hostSize() << " contacts passed.\n";
}

void configureFewPairs(Fixture& fixture)
{
    fixture.wall(); // Deliberately insert the wall first to test pair reordering.
    const int geometry = fixture.sphere(5, 0.46);
    fixture.add(geometry, {-0.35, 0.04, 0.34}, 0.31);
    fixture.add(geometry, {0.36, -0.03, 0.37}, -0.27);
    require(fixture.particles.host()[1].surfaceNodes().size() == 10242, "Few-pair fixture must retain level-5 surface resolution.");
}

void benchmark(Fixture& fixture)
{
    constexpr int repeats = 30;
    for (int threads : {1, 4})
    {
        setThreads(threads);
        contactSearch search;
        contactContainer contacts;
        search.findContacts(contacts, fixture.particles);
        const auto begin = std::chrono::steady_clock::now();
        for (int iteration = 0; iteration < repeats; ++iteration)
        {
            search.findContacts(contacts, fixture.particles);
        }
        const Real milliseconds = std::chrono::duration<Real, std::milli>(std::chrono::steady_clock::now() - begin).count() / repeats;
        std::cout << "benchmark threads=" << threads << " pairs=" << search.candidatePairCount() << " contacts=" << contacts.hostSize() << " search_ms=" << milliseconds << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        setThreads(1);
        Fixture fewPairs;
        configureFewPairs(fewPairs);
        exercise(fewPairs, "level-5 few pairs");
        Fixture mixed;
        mixed.wall();
        const int low = mixed.sphere(2, 0.43);
        const int medium = mixed.ellipsoid(3);
        const int high = mixed.sphere(5, 0.46);
        const int geometries[] = {high, low, medium};
        for (int index = 0; index < 8; ++index)
        {
            mixed.add(geometries[index % 3], {-1.03 + 0.69 * (index % 4), -0.34 + 0.67 * (index / 4), 0.34 + 0.012 * index}, 0.13 * (index + 1));
        }
        exercise(mixed, "mixed surface resolutions");
        if (argc > 1)
        {
            require(argc == 2 && std::string(argv[1]) == "--benchmark", "Usage: fundem_ls_contact_parallel_test [--benchmark]");
            benchmark(fewPairs);
        }
#if !defined(_OPENMP)
        std::cout << "OpenMP unavailable: both thread configurations used the serial fallback.\n";
#endif
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "LS contact parallel regression failed: " << error.what() << '\n';
        return 1;
    }
}
