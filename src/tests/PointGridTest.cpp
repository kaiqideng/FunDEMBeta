#include "data/myLSObject.h"
#include "execution/cpu/pointGrid.h"
#include "interaction/contactSearch.h"
#include "solver/SphereDEM.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace
{
using fundem::math::Real;
using fundem::math::Vec3;
using Pair = std::pair<int, int>;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

template <class Operation> void invalid(Operation operation, const std::string& message)
{
    bool rejected = false;
    try { operation(); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, message);
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

Real exactDistance(const Vec3& first, const Vec3& second)
{
    const auto delta = first - second;
    return std::hypot(delta.x, delta.y, delta.z);
}

template <class Pairs> std::vector<Pair> indices(const Pairs& pairs)
{
    std::vector<Pair> result;
    for (const auto& pair : pairs)
        result.emplace_back(pair.firstIndex_, pair.secondIndex_);
    require(std::is_sorted(result.begin(), result.end()), "Point pairs lost deterministic original-index order.");
    require(std::adjacent_find(result.begin(), result.end()) == result.end(), "A fixed-radius query returned a duplicate pair.");
    return result;
}

std::vector<Pair> oracle(const std::vector<Vec3>& first, const std::vector<Vec3>& second, Real radius, bool same)
{
    std::vector<Pair> result;
    for (std::size_t a = 0; a < first.size(); ++a)
        for (std::size_t b = same ? a + 1 : 0; b < second.size(); ++b)
            if (exactDistance(first[a], second[b]) <= radius)
                result.emplace_back(static_cast<int>(a), static_cast<int>(b));
    return result;
}

void compare(const std::vector<Vec3>& first, const std::vector<Vec3>& second, Real radius)
{
    const auto sameExpected = oracle(first, first, radius, true), crossExpected = oracle(first, second, radius, false);
    for (int threads : {1, 2, 4})
    {
        setThreads(threads);
        require(indices(fundem::cpu::findPointPairs(first, radius)) == sameExpected, "Same-set fixed-radius pairs differ from exact inclusive brute force.");
        require(indices(fundem::cpu::findPointPairs(first, second, radius)) == crossExpected, "Cross-set fixed-radius pairs differ from exact inclusive brute force.");
    }
}

void checkInclusiveAndZeroRadius()
{
    const auto inside = std::nextafter(1.0, 0.0), outside = std::nextafter(1.0, std::numeric_limits<Real>::infinity());
    const std::vector<Vec3> points{{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0},
                                   {0, 0, 1}, {0, 0, -1}, {inside, 0, 0}, {outside, 0, 0},
                                   {0.6, 0.8, 0}, {0, 0, 0}, {-2, -2, -2}, {-1, -2, -2}};
    const std::vector<Vec3> other{{0, 0, 0}, {0, 0, 0}, {-1, 0, 0}, {2, 0, 0}, {-2, -2, -2}};
    compare(points, other, 1.0);
    compare(points, other, 0.0);
    compare(std::vector<Vec3>{{0, 0, 0}, {std::nextafter(0.0, 1.0), 0, 0}, {0, 0, 0}},
            std::vector<Vec3>{{0, 0, 0}, {std::nextafter(0.0, 1.0), 0, 0}}, 0.0);
    require(indices(fundem::cpu::findPointPairs(points, 0.0)) == std::vector<Pair>{{0, 10}},
            "Zero-radius same-set queries must retain only distinct coincident indices.");
    require(indices(fundem::cpu::findPointPairs(std::vector<Vec3>{{0, 0, 0}, {1, 0, 0}}, 1.0)) == std::vector<Pair>{{0, 1}},
            "A pair exactly at the maximum distance was excluded.");
    require(fundem::cpu::findPointPairs(std::vector<Vec3>{{0, 0, 0}, {outside, 0, 0}}, 1.0).empty(),
            "A pair one representable value beyond the threshold was included.");
}

void checkRandomAndCandidateCompleteness()
{
    std::mt19937_64 generator(0xC0111DE);
    std::uniform_real_distribution<Real> coordinate(-5.0, 5.0);
    std::vector<Vec3> first, second;
    for (int index = 0; index < 130; ++index)
        first.push_back({coordinate(generator), coordinate(generator), coordinate(generator)});
    for (int index = 0; index < 95; ++index)
        second.push_back({coordinate(generator), coordinate(generator), coordinate(generator)});
    first.insert(first.end(), {{-1, -1, -1}, {-1, -1, -1}, {0, 0, 0}, {1, 1, 1}});
    compare(first, second, 1.25);
    compare(first, first, 0.5); // Cross-set semantics retain both directions and equal local indices.
    for (Real cellSize : {0.25, 1.0, 3.0})
    {
        fundem::cpu::pointGrid grid;
        grid.build(first, cellSize);
        for (const Vec3 query : {Vec3{-1, -1, -1}, Vec3{}, Vec3{1, 1, 1}, Vec3{-5.01, 0, 0}})
            for (Real radius : {0.0, 0.5, 2.0})
            {
                std::vector<int> actual, candidates;
                grid.forEachCandidate(query, radius, [&](int index) {
                    require(index >= 0 && index < static_cast<int>(first.size()), "Grid returned an index outside its original input.");
                    candidates.push_back(index);
                    if (exactDistance(query, first[index]) <= radius)
                        actual.push_back(index);
                });
                std::sort(candidates.begin(), candidates.end());
                require(std::adjacent_find(candidates.begin(), candidates.end()) == candidates.end(), "A point was visited more than once in one query.");
                std::sort(actual.begin(), actual.end());
                std::vector<int> expected;
                for (int index = 0; index < static_cast<int>(first.size()); ++index)
                    if (exactDistance(query, first[index]) <= radius)
                        expected.push_back(index);
                require(actual == expected, "Cell boundary/radius expansion omitted an exact neighbor.");
            }
    }
}

void checkEmptyInvalidAndExtreme()
{
    const std::vector<Vec3> empty, one{{0, 0, 0}};
    require(fundem::cpu::findPointPairs(empty, 1.0).empty() && fundem::cpu::findPointPairs(one, empty, 1.0).empty() &&
                fundem::cpu::findPointPairs(empty, one, 1.0).empty(), "An empty query population produced pairs.");
    fundem::cpu::pointGrid grid;
    grid.build(empty, 1.0);
    require(grid.occupiedCellCount() == 0, "An empty grid retained occupied cells.");
    grid.forEachCandidate({}, 0.0, [](int) { throw std::runtime_error("An empty grid visited an index."); });
    const Real infinity = std::numeric_limits<Real>::infinity(), nan = std::numeric_limits<Real>::quiet_NaN();
    for (Real radius : {-1.0, infinity, nan})
    {
        invalid([&] { fundem::cpu::findPointPairs(one, radius); }, "Same-set query accepted an invalid radius.");
        invalid([&] { fundem::cpu::findPointPairs(one, empty, radius); }, "An empty cross-set query bypassed invalid-radius validation.");
        invalid([&] { grid.forEachCandidate({}, radius, [](int) {}); }, "Grid accepted an invalid query radius.");
    }
    for (Real cellSize : {0.0, -1.0, infinity, nan})
        invalid([&] { grid.build(one, cellSize); }, "Grid accepted a nonpositive/nonfinite cell size.");
    for (const Vec3 point : {Vec3{nan, 0, 0}, Vec3{0, infinity, 0}})
    {
        invalid([&] { grid.build({point}, 1.0); }, "Grid accepted a nonfinite point.");
        invalid([&] { fundem::cpu::findPointPairs({point}, one, 1.0); }, "Point-pair validation missed a nonfinite first-set point.");
        invalid([&] { fundem::cpu::findPointPairs(one, {point}, 1.0); }, "Point-pair validation missed a nonfinite second-set point.");
        invalid([&] { grid.forEachCandidate(point, 1.0, [](int) {}); }, "Grid accepted a nonfinite query position.");
    }
    // Finite input outside int64 cell coordinates must fall back conservatively,
    // never wrap indices, lose pairs, or enumerate an astronomical empty range.
    const std::vector<Vec3> extremes{{1e100, 0, 0}, {1e100, 0, 0}, {-1e100, 0, 0}, {0, 0, 0}};
    compare(extremes, std::vector<Vec3>{{1e100, 0, 0}, {0, 0, 0}}, 1.0);
    compare(std::vector<Vec3>{{-1e18, 0, 0}, {1, 2, 3}, {1.1, 2, 3}, {1.2, 2, 3}},
            std::vector<Vec3>{{1.05, 2, 3}, {1.15, 2, 3}}, 0.15);
    compare(std::vector<Vec3>{{-1e10, 0, 0}, {1e10, 0, 0}, {10000000000.000006, 0, 0}},
            std::vector<Vec3>{{1e10, 0, 0}, {10000000000.000006, 0, 0}}, 5.722045898437501e-6);
    compare(std::vector<Vec3>{{-1e10, 0, 0}, {1e-7, 0, 0}}, std::vector<Vec3>{{1e-7, 0, 0}}, 1e10);
    grid.build(one, 1.0);
    int visited = 0;
    grid.forEachCandidate({}, 1e100, [&](int) { ++visited; });
    require(visited == 1, "A huge-radius query lost or repeated its only occupied point.");
    grid.build({Vec3{2, 3, 4}, Vec3{2, 3, 4}}, std::numeric_limits<Real>::denorm_min());
    visited = 0;
    grid.forEachCandidate({2, 3, 4}, 0.0, [&](int) { ++visited; });
    require(visited == 2, "Very small positive cell size overflow lost coincident candidates.");
    // A grid retains its indices independently from temporary input storage.
    {
        const std::vector<Vec3> temporary{{4, 5, 6}, {-1, -2, -3}};
        grid.build(temporary, 1.0);
    }
    visited = 0;
    grid.forEachCandidate({4, 5, 6}, 0.0, [&](int index) { require(index == 0, "Retained cell index changed after source destruction."); ++visited; });
    require(visited == 1, "Grid borrowed a destroyed position array.");
}

void checkSparseCapacity()
{
    std::vector<Vec3> points;
    constexpr int sites = 10000;
    for (int index = 0; index < sites; ++index)
    {
        const Vec3 origin{10.0 * (index % 50), 10.0 * ((index / 50) % 20), 10.0 * (index / 1000)};
        points.push_back(origin);
        points.push_back(origin + Vec3{0.1, 0, 0});
    }
    fundem::cpu::pointGrid grid;
    grid.build(points, 0.2);
    std::size_t candidates = 0;
    for (const auto& point : points)
        grid.forEachCandidate(point, 0.2, [&](int) { ++candidates; });
    require(candidates <= 4 * points.size(), "Sparse fixed-radius queries degenerated into a full all-pairs scan.");
    std::vector<Pair> baseline;
    for (int threads : {1, 2, 4})
    {
        setThreads(threads);
        const auto result = indices(fundem::cpu::findPointPairs(points, 0.2));
        require(result.size() == sites, "A large sparse scene lost a close pair or included a distant site.");
        for (int index = 0; index < sites; ++index)
            require(result[index] == Pair{2 * index, 2 * index + 1}, "Large-scene original indices or ordering changed.");
        if (baseline.empty()) baseline = result;
        else require(result == baseline, "Thread count changed the fixed-radius pair order.");
    }
    std::cout << "Sparse point-grid candidates=" << candidates << ", all-pair upper bound=" << points.size() * points.size() << '\n';
}

using ContactKey = std::tuple<int, int, int>;
ContactKey contactKey(const fundem::contact& value)
{
    return {value.masterParticleIndex(), value.slaveParticleIndex(), value.particleSurfaceNodeIndex()};
}

void sameContacts(const fundem::contactContainer& first, const fundem::contactContainer& second)
{
    require(first.hostSize() == second.hostSize(), "Reused point-grid search changed the native contact count.");
    std::map<ContactKey, const fundem::contact*> expected;
    for (const auto& value : first.host())
        require(expected.emplace(contactKey(value), &value).second, "The direct contact oracle repeated an endpoint/node identity.");
    for (const auto& b : second.host())
    {
        const auto match = expected.find(contactKey(b));
        require(match != expected.end(), "The grid search introduced an endpoint/node identity absent from the direct oracle.");
        const auto& a = *match->second;
        require(exactDistance(a.point(), b.point()) < 1e-12 && exactDistance(a.normal(), b.normal()) < 1e-12 &&
                    std::abs(a.overlap() - b.overlap()) < 1e-12 && std::abs(a.effectiveMass() - b.effectiveMass()) < 1e-12 &&
                    exactDistance(a.slidingSpringDeformation(), b.slidingSpringDeformation()) < 1e-12,
                "Reused point-grid search changed contact geometry, mass or history.");
        expected.erase(match);
    }
}

void checkContactSearchRegression()
{
    fundem::SphereDEM solver;
    solver.addMaterial(fundem::material(1e8, 0.25, 0.9, 1000));
    for (int index = 0; index < 10; ++index)
    {
        fundem::particle value;
        value.setPosition({0.35 * (index % 5) - 0.7, 0.35 * (index / 5) - 0.2, 0});
        value.setRadius(0.2 + 0.02 * (index % 3));
        value.setMaterial(solver.materials(), 0);
        if (index == 0 || index == 9) value.setInfiniteMass();
        solver.addSphere(value);
    }
    fundem::levelset::Sphere geometry(0.5);
    geometry.buildSurfaceNode(1);
    geometry.buildLSGrid(0.1, 4, false);
    const int geometryIndex = solver.addGeometry(geometry);
    for (int index = 0; index < 4; ++index)
    {
        fundem::LSParticle value;
        value.setPosition({-0.8 + 0.6 * index, 0, 0});
        value.setGeometry(solver.geometries(), geometryIndex);
        value.setMaterial(solver.materials(), 0);
        if (index == 0) value.setInfiniteMass();
        solver.addLSParticle(value);
    }
    fundem::contactSearch direct, grid(0.3);
    direct.setDirectSearch(true);
    fundem::contactContainer expected, actual;
    const auto compareSearch = [&](auto search) {
        search(direct, expected); search(grid, actual);
        require(expected.hostSize() > 0, "Contact regression fixture has no actual contacts.");
        sameContacts(expected, actual);
        std::vector<ContactKey> order;
        for (const auto& value : actual.host()) order.push_back(contactKey(value));
        for (auto* contacts : {&expected, &actual})
            for (auto& contact : contacts->host())
                contact.setSlidingSpringDeformation({1e-4 * (contact.masterParticleIndex() + 1), -2e-4, 3e-4});
        search(direct, expected); search(grid, actual);
        sameContacts(expected, actual);
        std::vector<ContactKey> repeated;
        for (const auto& value : actual.host()) repeated.push_back(contactKey(value));
        require(order == repeated, "Repeated contact-grid queries changed stable pair/node traversal order.");
    };
    compareSearch([&](auto& search, auto& contacts) { search.findContacts(contacts, solver.spheres()); });
    compareSearch([&](auto& search, auto& contacts) { search.findContacts(contacts, solver.spheres(), solver.LSParticles()); });
    compareSearch([&](auto& search, auto& contacts) { search.findContacts(contacts, solver.LSParticles()); });

    fundem::SphereDEM separated;
    separated.addMaterial(fundem::material(1e8, 0.25, 0.9, 1000));
    for (Real x : {0.0, 1.5})
    {
        fundem::particle particle;
        particle.setPosition({x, 0, 0}); particle.setRadius(0.25); particle.setMaterial(separated.materials(), 0);
        separated.addSphere(particle);
    }
    require(fundem::cpu::findPointPairs(std::vector<Vec3>{{0, 0, 0}, {1.5, 0, 0}}, 1.6).size() == 1 &&
                grid.findContacts(actual, separated.spheres()) == 0, "Fixed-radius near neighbors were incorrectly restricted to physical contacts.");
}
} // namespace

int main()
{
    try
    {
        checkInclusiveAndZeroRadius();
        checkRandomAndCandidateCompleteness();
        checkEmptyInvalidAndExtreme();
        checkSparseCapacity();
        checkContactSearchRegression();
        std::cout << "POINT_GRID_CHECKS_PASSED: exact same/cross pairs, stable indices, guarded extremes, sparse candidates and unchanged contacts.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
