/** @file LSPrimitiveTest.cpp @brief Topology, distance and prepared-frame checks for finite LS primitives. */
#include "data/myLSObject.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{
using namespace fundem;
using math::Real;
using math::Vec3;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void near(Real actual, Real expected, Real tolerance, const char* message)
{
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}

template<class Function> void rejects(Function&& function, const char* message)
{
    bool rejected = false;
    try { function(); }
    catch (const std::exception&) { rejected = true; }
    require(rejected, message);
}

/** Every welded edge must have exactly two opposite uses; nodes must lie on the native analytic boundary. */
template<class NormalFunction> void checkSurface(const levelset::LSInfo& shape, Real spacing, int eulerCharacteristic, NormalFunction&& outward)
{
    const auto& nodes = shape.surfaceNodePosition();
    const auto& faces = shape.surfaceNodeConnectivity();
    require(!nodes.empty() && !faces.empty(), "Primitive generated an empty surface.");
    std::map<std::pair<int, int>, std::pair<int, int>> edges;
    std::vector<bool> used(nodes.size(), false);
    for (const Vec3& point : nodes)
        near(shape.signedDistance(point), 0.0, 2.0e-12, "Primitive vertex is off the exact surface.");
    for (const int3& face : faces)
    {
        const int indices[]{face.x, face.y, face.z};
        for (int index : indices)
        {
            require(index >= 0 && index < static_cast<int>(nodes.size()), "Primitive connectivity is out of range.");
            used[index] = true;
        }
        const Vec3 normal = math::cross(nodes[face.y] - nodes[face.x], nodes[face.z] - nodes[face.x]);
        require(math::norm(normal) > 0.0 && math::isFinite(normal), "Primitive has a degenerate face.");
        const Vec3 center = (nodes[face.x] + nodes[face.y] + nodes[face.z]) / 3.0;
        require(math::dot(normal, outward(center)) > 0.0, "Primitive face winding points inward.");
        for (int i = 0; i < 3; ++i)
        {
            const int first = indices[i];
            const int second = indices[(i + 1) % 3];
            auto& edge = edges[std::minmax(first, second)];
            ++edge.first;
            edge.second += first < second ? 1 : -1;
            require(math::norm(nodes[first] - nodes[second]) <= spacing * (1.0 + 1.0e-12), "Primitive mesh exceeds target edge spacing.");
        }
    }
    for (const auto& entry : edges)
        require(entry.second.first == 2 && entry.second.second == 0, "Primitive surface has an open, non-manifold or inconsistently oriented edge.");
    require(std::all_of(used.begin(), used.end(), [](bool value) { return value; }), "Primitive has an unused surface vertex.");
    require(static_cast<int>(nodes.size()) - static_cast<int>(edges.size()) + static_cast<int>(faces.size()) == eulerCharacteristic,
            "Primitive surface has incorrect topology.");
}

/** Grid samples and analytic queries must remain in one centroid-corrected local frame after repeated builds. */
void checkGrid(levelset::LSInfo& shape, Real expectedVolume)
{
    const Vec3 nativeMinimum = shape.gridSamplingBoundsMin();
    const Vec3 nativeMaximum = shape.gridSamplingBoundsMax();
    const Vec3 nativeFirstNode = shape.surfaceNodePosition().front();
    shape.buildLSGrid(0.05, 2);
    require(shape.volume() > 0.0 && std::isfinite(shape.boundingRadius()), "Primitive did not compute valid mass metadata.");
    near(shape.volume(), expectedVolume, expectedVolume * 0.15, "Primitive integrated volume disagrees with analytic volume.");
    const auto inertia = shape.unitDensityInertiaTensor();
    require(math::isFinite(inertia) && inertia(0, 0) > 0.0 && inertia(1, 1) > 0.0 && inertia(2, 2) > 0.0,
            "Primitive integrated inertia is invalid.");
    const int3 size = shape.gridNodeSize3D();
    const auto origin = shape.gridNodeOrigin();
    const auto values = shape.gridNodeSFD();
    const auto surface = shape.surfaceNodePosition();
    const Vec3 correction = nativeFirstNode - surface.front();
    near(math::norm(shape.gridSamplingBoundsMin() - (nativeMinimum - correction)), 0.0, 1.0e-12,
         "Primitive lower sampling extent lost its centroid-corrected frame.");
    near(math::norm(shape.gridSamplingBoundsMax() - (nativeMaximum - correction)), 0.0, 1.0e-12,
         "Primitive upper sampling extent lost its centroid-corrected frame.");
    for (std::size_t i = 0; i < values.size(); i += 71)
    {
        const Vec3 point = origin + shape.gridNodeSpacing() * Vec3{Real(i % size.x), Real((i / size.x) % size.y), Real(i / (size.x * size.y))};
        near(shape.signedDistance(point), values[i], 1.0e-12, "Primitive grid and analytic query use different local frames.");
    }
    shape.buildLSGrid(0.05, 2);
    require(values == shape.gridNodeSFD(), "Primitive rebuild changed native grid samples.");
    near(math::norm(shape.gridNodeOrigin() - origin), 0.0, 1.0e-12, "Primitive repeated centroid correction moved the grid twice.");
    require(surface.size() == shape.surfaceNodePosition().size(), "Primitive grid rebuild changed mesh resolution.");
    for (std::size_t i = 0; i < surface.size(); ++i)
        near(math::norm(surface[i] - shape.surfaceNodePosition()[i]), 0.0, 1.0e-12, "Primitive repeated centroid correction moved its surface twice.");
}

void testPolyhedra()
{
    constexpr Real spacing = 0.35;
    const auto radialNormal = [](const Vec3& point) { return point; };
    levelset::HexagonalPrism hexagon(1.0, 1.5);
    levelset::TriangularPrism triangle(1.0, 1.5);
    levelset::Octahedron octahedron({1.0, 0.75, 0.5});
    levelset::Pyramid pyramid(2.0, 1.5, 2.0);
    levelset::ConvexPolyhedron* polyhedra[]{&hexagon, &triangle, &octahedron, &pyramid};
    for (levelset::ConvexPolyhedron* shape : polyhedra)
    {
        shape->buildSurfaceNode(spacing);
        checkSurface(*shape, spacing, 2, radialNormal);
        require(shape->signedDistance(Vec3::zero()) < 0.0, "Convex primitive origin must be interior.");
        require(shape->signedDistance({10.0, 0.0, 0.0}) > 0.0, "Convex primitive exterior has wrong distance sign.");
        rejects([&] { shape->buildSurfaceNode(0.0); }, "Convex primitive accepted zero surface spacing.");
        rejects([&] { shape->buildSurfaceNode(1.0e-12); }, "Convex primitive did not enforce its mesh budget.");
    }
    near(hexagon.signedDistance({0.0, 0.0, 1.0}), 0.25, 1.0e-12, "Hexagonal prism cap distance is wrong.");
    near(hexagon.signedDistance({2.0, 0.0, 0.75}), 1.0, 1.0e-12, "Hexagonal prism edge distance is wrong.");
    near(triangle.signedDistance({-1.0, 0.0, 0.0}), 0.5, 1.0e-12, "Triangular prism planar side distance is wrong.");
    near(octahedron.signedDistance({2.0, 0.0, 0.0}), 1.0, 1.0e-12, "Octahedron vertex distance is wrong.");
    near(octahedron.signedDistance(Vec3::zero()), -1.0 / std::sqrt(1.0 + 1.0 / (0.75 * 0.75) + 4.0), 1.0e-12,
         "Anisotropic octahedron distance is not Euclidean.");
    near(pyramid.signedDistance({0.0, 0.0, -1.25}), 0.25, 1.0e-12, "Pyramid base distance is wrong.");
    checkGrid(hexagon, 1.5 * 3.0 * std::sqrt(3.0) / 2.0);
    checkGrid(triangle, 1.5 * 3.0 * std::sqrt(3.0) / 4.0);
    checkGrid(octahedron, 4.0 * 1.0 * 0.75 * 0.5 / 3.0);
    checkGrid(pyramid, 2.0 * 1.5 * 2.0 / 3.0);
    Real maximumZ = -std::numeric_limits<Real>::infinity();
    for (const Vec3& node : pyramid.surfaceNodePosition()) maximumZ = std::max(maximumZ, node.z);
    near(maximumZ, 1.5, 0.05, "Pyramid centroid is not at one-quarter height above the base.");
    pyramid.buildSurfaceNode(spacing);
    require(pyramid.gridNodeSFD().empty() && pyramid.volume() == 0.0, "Regenerating a primitive surface retained stale integrated data.");
    near(pyramid.gridSamplingBoundsMin().z, -1.0, 1.0e-12, "Pyramid lower sampling bound did not return to its native frame.");
    near(pyramid.gridSamplingBoundsMax().z, 1.0, 1.0e-12, "Pyramid upper sampling bound did not return to its native frame.");
    pyramid.buildLSGrid(0.05, 2, true);
    require(pyramid.volume() == 0.0 && pyramid.unitDensityInertiaTensor() == math::Mat3::zero(), "Fixed primitive did not skip mass integration.");
    near(pyramid.signedDistance({0.0, 0.0, 1.0}), 0.0, 1.0e-12, "Fixed pyramid lost its native apex position.");
}

void testRegularSolidsAndPrisms()
{
    constexpr Real spacing = 0.3;
    const auto outward = [](const Vec3& point) { return point; };
    levelset::Tetrahedron tetrahedron(1.0);
    levelset::Dodecahedron dodecahedron(1.0);
    levelset::Icosahedron icosahedron(1.0);
    require(tetrahedron.surfaceNodePosition().size() == 4 && tetrahedron.surfaceNodeConnectivity().size() == 4,
            "Regular tetrahedron base topology is wrong.");
    require(dodecahedron.surfaceNodePosition().size() == 20 && dodecahedron.surfaceNodeConnectivity().size() == 36,
            "Dodecahedron pentagonal facets were not uniquely triangulated.");
    require(icosahedron.surfaceNodePosition().size() == 12 && icosahedron.surfaceNodeConnectivity().size() == 20,
            "Regular icosahedron base topology is wrong.");
    levelset::ConvexPolyhedron* solids[]{&tetrahedron, &dodecahedron, &icosahedron};
    for (auto* shape : solids)
    {
        for (const auto& point : shape->surfaceNodePosition())
            near(math::norm(point), 1.0, 1.0e-12, "Regular solid circumradius does not match its vertices.");
        shape->buildSurfaceNode(spacing);
        checkSurface(*shape, spacing, 2, outward);
        require(shape->signedDistance(Vec3::zero()) < 0.0, "Regular solid origin is not inside.");
        near(shape->signedDistance(1.5 * shape->surfaceNodePosition().front()), 0.5, 1.0e-12,
             "Regular solid has an incorrect exact distance beyond a circumsphere vertex.");
    }
    const Real dodecahedronEdge = 4.0 / (std::sqrt(3.0) * (1.0 + std::sqrt(5.0)));
    const Real icosahedronEdge = 4.0 / std::sqrt(10.0 + 2.0 * std::sqrt(5.0));
    checkGrid(tetrahedron, 8.0 / (9.0 * std::sqrt(3.0)));
    checkGrid(dodecahedron, (15.0 + 7.0 * std::sqrt(5.0)) * std::pow(dodecahedronEdge, 3) / 4.0);
    checkGrid(icosahedron, 5.0 * (3.0 + std::sqrt(5.0)) * std::pow(icosahedronEdge, 3) / 12.0);
    for (int sides : {3, 5, 6, 8})
    {
        levelset::RegularPolygonPrism prism(1.0, 1.5, sides);
        prism.buildSurfaceNode(spacing);
        checkSurface(prism, spacing, 2, outward);
        near(prism.signedDistance({0.0, 0.0, 1.0}), 0.25, 1.0e-12, "General regular prism has the wrong cap distance.");
        checkGrid(prism, 0.5 * sides * std::sin(2.0 * math::pi / sides) * 1.5);
    }
    levelset::HexagonalPrism legacyHexagon(1.0, 1.5);
    levelset::TriangularPrism legacyTriangle(1.0, 1.5);
    for (int sides : {3, 6})
    {
        levelset::RegularPolygonPrism prism(1.0, 1.5, sides);
        const levelset::LSInfo& legacy = sides == 3 ? static_cast<const levelset::LSInfo&>(legacyTriangle)
                                                   : static_cast<const levelset::LSInfo&>(legacyHexagon);
        for (const Vec3& point : {Vec3::zero(), Vec3{1.5, 0.2, 0.3}, Vec3{-0.5, 0.7, 1.0}})
            near(prism.signedDistance(point), legacy.signedDistance(point), 1.0e-12, "General prism changed the legacy triangle/hexagon shape.");
    }
    rejects([] { levelset::RegularPolygonPrism shape(1.0, 1.0, 2); }, "Regular prism accepted fewer than three sides.");
    rejects([] { levelset::RegularPolygonPrism shape(1.0, 1.0, std::numeric_limits<int>::max()); }, "Prism base budget was not checked before allocation.");
    rejects([] { levelset::Tetrahedron shape(0.0); }, "Tetrahedron accepted zero radius.");
    rejects([] { levelset::Dodecahedron shape(-1.0); }, "Dodecahedron accepted negative radius.");
    rejects([] { levelset::Icosahedron shape(std::numeric_limits<Real>::quiet_NaN()); }, "Icosahedron accepted NaN radius.");
}

void testCutSphereAndLink()
{
    constexpr Real spacing = 0.3;
    for (Real cutHeight : {-0.6, 0.0, 0.6})
    {
        levelset::CutSphere sphere(1.0, cutHeight);
        sphere.buildSurfaceNode(spacing);
        checkSurface(sphere, spacing, 2, [cutHeight](const Vec3& point)
        { return std::abs(point.z - cutHeight) < 1.0e-12 ? Vec3::unitZ() : point; });
        near(sphere.signedDistance({0.0, 0.0, cutHeight + 0.2}), 0.2, 1.0e-12, "Cut sphere planar-cap distance is wrong.");
        require(sphere.signedDistance({0.0, 0.0, -0.8}) < 0.0, "Cut sphere retained the wrong half-space.");
        const Real capHeight = 1.0 + cutHeight;
        checkGrid(sphere, math::pi * capHeight * capHeight * (1.0 - capHeight / 3.0));
        if (cutHeight == 0.0)
            near(sphere.gridSamplingBoundsMax().z, 0.375, 0.025, "Hemisphere centroid correction is wrong.");
        sphere.buildSurfaceNode(spacing);
        near(sphere.gridSamplingBoundsMax().z, cutHeight, 1.0e-12, "Cut sphere regeneration retained an old centroid shift.");
        rejects([&] { sphere.buildSurfaceNode(1.0e-12); }, "Cut sphere surface budget was not enforced.");
    }
    levelset::CutSphere hemisphere(1.0, 0.0);
    near(hemisphere.signedDistance({2.0, 0.0, 1.0}), std::sqrt(2.0), 1.0e-12,
         "Cut sphere used a CSG max bound instead of exact distance at the rim.");
    near(hemisphere.signedDistance({0.0, 0.0, -0.25}), -0.25, 1.0e-12, "Hemisphere interior distance is wrong.");
    rejects([] { levelset::CutSphere shape(1.0, 1.0); }, "Cut sphere accepted a degenerate upper cut.");
    rejects([] { levelset::CutSphere shape(1.0, -1.0); }, "Cut sphere accepted an empty lower cut.");
    rejects([] { levelset::CutSphere shape(0.0, 0.0); }, "Cut sphere accepted zero radius.");
    for (Real straightLength : {0.0, 1.5})
    {
        levelset::Link link(1.0, 0.2, straightLength);
        link.buildSurfaceNode(spacing);
        checkSurface(link, spacing, 0, [straightLength](const Vec3& point)
        {
            const Real y = std::max(std::abs(point.y) - straightLength * 0.5, Real{0.0});
            const Real radial = std::hypot(point.x, y);
            return Vec3{point.x * (1.0 - 1.0 / radial), std::copysign(y, point.y) * (1.0 - 1.0 / radial), point.z};
        });
        near(link.signedDistance(Vec3::zero()), 0.8, 1.0e-12, "Link cavity was incorrectly filled.");
        near(link.signedDistance({1.0, 0.0, 0.0}), -0.2, 1.0e-12, "Link straight-side tube distance is wrong.");
        near(link.signedDistance({0.0, 0.5 * straightLength + 1.5, 0.0}), 0.3, 1.0e-12, "Link curved-end distance is wrong.");
        if (straightLength == 0.0)
        {
            levelset::Torus torus(1.0, 0.2);
            for (int i = 0; i < 30; ++i)
            {
                const Vec3 point{0.1 * i - 1.2, std::sin(Real(i)), std::cos(Real(i))};
                near(link.signedDistance(point), torus.signedDistance(point), 1.0e-12, "Zero-length link does not reduce to a torus.");
            }
        }
        checkGrid(link, math::pi * 0.2 * 0.2 * (2.0 * straightLength + 2.0 * math::pi));
        link.buildSurfaceNode(spacing);
        require(link.volume() == 0.0 && link.gridNodeSFD().empty(), "Link regeneration retained old grid metadata.");
        rejects([&] { link.buildSurfaceNode(1.0e-12); }, "Link surface budget was not enforced.");
    }
    rejects([] { levelset::Link shape(1.0, 1.0, 1.0); }, "Link accepted a tube as large as its bend radius.");
    rejects([] { levelset::Link shape(1.0, 0.1, -1.0); }, "Link accepted a negative straight length.");
}

void testTorusAndCapsule()
{
    constexpr Real spacing = 0.3;
    levelset::Torus torus(1.0, 0.25);
    torus.buildSurfaceNode(spacing);
    checkSurface(torus, spacing, 0, [](const Vec3& point)
    {
        const Real radial = std::hypot(point.x, point.y);
        return Vec3{point.x * (1.0 - 1.0 / radial), point.y * (1.0 - 1.0 / radial), point.z};
    });
    near(torus.signedDistance(Vec3::zero()), 0.75, 1.0e-12, "Torus hole was incorrectly filled.");
    near(torus.signedDistance({1.0, 0.0, 0.0}), -0.25, 1.0e-12, "Torus tube interior distance is wrong.");
    near(torus.signedDistance({1.5, 0.0, 0.0}), 0.25, 1.0e-12, "Torus outer distance is wrong.");
    checkGrid(torus, 2.0 * math::pi * math::pi * 0.25 * 0.25);
    for (Real length : {1.5, 0.0})
    {
        levelset::Capsule capsule(0.5, length);
        capsule.buildSurfaceNode(spacing);
        checkSurface(capsule, spacing, 2, [length](const Vec3& point)
        { return point - Vec3{0.0, 0.0, std::clamp(point.z, -0.5 * length, 0.5 * length)}; });
        near(capsule.signedDistance(Vec3::zero()), -0.5, 1.0e-12, "Capsule interior distance is wrong.");
        near(capsule.signedDistance({0.0, 0.0, 0.5 * length + 0.75}), 0.25, 1.0e-12, "Capsule end distance is wrong.");
        checkGrid(capsule, math::pi * 0.25 * length + 4.0 * math::pi * 0.125 / 3.0);
        capsule.buildSurfaceNode(spacing);
        require(capsule.volume() == 0.0 && capsule.gridNodeSFD().empty(), "Capsule surface rebuild retained grid state.");
        rejects([&] { capsule.buildSurfaceNode(1.0e-12); }, "Capsule mesh budget was not enforced.");
    }
    rejects([&] { torus.buildSurfaceNode(1.0e-12); }, "Torus mesh budget was not enforced.");
    rejects([&] { torus.buildSurfaceNode(-0.1); }, "Torus accepted negative surface spacing.");
    rejects([] { levelset::Capsule shape; shape.buildSurfaceNode(0.1); }, "Unconfigured capsule generated a surface.");
    rejects([] { levelset::HexagonalPrism shape(0.0, 1.0); }, "Prism accepted zero circumradius.");
    rejects([] { levelset::TriangularPrism shape(1.0, -1.0); }, "Prism accepted negative height.");
    rejects([] { levelset::Octahedron shape({1.0, 0.0, 1.0}); }, "Octahedron accepted zero half-length.");
    rejects([] { levelset::Pyramid shape(1.0, 1.0, 0.0); }, "Pyramid accepted zero height.");
    rejects([] { levelset::Torus shape(1.0, 1.0); }, "Horn torus was accepted.");
    rejects([] { levelset::Torus shape(0.5, 1.0); }, "Spindle torus was accepted.");
    rejects([] { levelset::Capsule shape(1.0, -0.1); }, "Capsule accepted negative straight length.");
    rejects([] { levelset::Capsule shape(std::numeric_limits<Real>::infinity(), 0.0); }, "Capsule accepted an infinite radius.");
}
} // namespace

int main()
{
    try
    {
        testPolyhedra();
        testTorusAndCapsule();
        testRegularSolidsAndPrisms();
        testCutSphereAndLink();
        std::cout << "Finite LS primitive geometry checks passed.\n";
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
