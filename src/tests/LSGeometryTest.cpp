#include "data/myLSObject.h"
#include "geometry/LSGeometry.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
using namespace fundem;
using math::Mat3;
using math::Real;
using math::Vec3;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void requireNear(Real actual, Real expected, const char* message, Real tolerance = 1.0e-10)
{
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}

void requireNear(const Vec3& actual, const Vec3& expected, const char* message, Real tolerance = 1.0e-10)
{
    requireNear(actual.x, expected.x, message, tolerance);
    requireNear(actual.y, expected.y, message, tolerance);
    requireNear(actual.z, expected.z, message, tolerance);
}

void requireNear(const Mat3& actual, const Mat3& expected, const char* message, Real tolerance = 1.0e-10)
{
    for (int i = 0; i < 9; ++i) requireNear(actual[i], expected[i], message, tolerance);
}

void requireSameSurface(const std::vector<Vec3>& actual, const std::vector<Vec3>& expected)
{
    require(actual.size() == expected.size(), "Surface node count changed unexpectedly.");
    for (std::size_t i = 0; i < actual.size(); ++i)
        requireNear(actual[i], expected[i], "Surface frame changed unexpectedly.");
}

struct GeometryState
{
    Vec3 origin;
    int3 size;
    Real spacing;
    Real radius;
    Real volume;
    Mat3 inertia;
    std::vector<Real> distances;
    std::vector<Vec3> surface;

    explicit GeometryState(const levelset::LSInfo& geometry)
        : origin(geometry.gridNodeOrigin()), size(geometry.gridNodeSize3D()),
          spacing(geometry.gridNodeSpacing()), radius(geometry.boundingRadius()),
          volume(geometry.volume()), inertia(geometry.unitDensityInertiaTensor()),
          distances(geometry.gridNodeSFD()), surface(geometry.surfaceNodePosition())
    {
    }

    void requireMatches(const levelset::LSInfo& geometry) const
    {
        requireNear(geometry.gridNodeOrigin(), origin, "Rebuilding changed the grid origin.");
        const int3 actualSize = geometry.gridNodeSize3D();
        require(actualSize.x == size.x && actualSize.y == size.y && actualSize.z == size.z,
                "Rebuilding changed the grid dimensions.");
        requireNear(geometry.gridNodeSpacing(), spacing, "Rebuilding changed grid spacing.");
        requireNear(geometry.boundingRadius(), radius, "Rebuilding changed bounding radius.");
        requireNear(geometry.volume(), volume, "Rebuilding changed integrated volume.");
        requireNear(geometry.unitDensityInertiaTensor(), inertia, "Rebuilding changed integrated inertia.");
        require(geometry.gridNodeSFD().size() == distances.size(), "Rebuilding changed sample count.");
        for (std::size_t i = 0; i < distances.size(); ++i)
            requireNear(geometry.gridNodeSFD()[i], distances[i], "Rebuilding changed signed distances.");
        requireSameSurface(geometry.surfaceNodePosition(), surface);
    }
};

levelset::TriangleMesh translatedBox(const Vec3& offset)
{
    std::vector<Vec3> nodes{{-1.0, -0.5, -0.5}, {1.0, -0.5, -0.5}, {1.0, 0.5, -0.5}, {-1.0, 0.5, -0.5},
                            {-1.0, -0.5, 0.5}, {1.0, -0.5, 0.5}, {1.0, 0.5, 0.5}, {-1.0, 0.5, 0.5}};
    for (Vec3& node : nodes) node += offset;
    return levelset::TriangleMesh(nodes, {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
                                         {0, 1, 5}, {0, 5, 4}, {3, 7, 6}, {3, 6, 2},
                                         {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}});
}

void requireQueryConsistency(const levelset::LSInfo& geometry)
{
    const int3 size = geometry.gridNodeSize3D();
    const Real spacing = geometry.gridNodeSpacing();
    const auto& distances = geometry.gridNodeSFD();
    for (std::size_t i = 0; i < distances.size(); i += 17)
    {
        const int x = static_cast<int>(i % size.x);
        const int y = static_cast<int>((i / size.x) % size.y);
        const int z = static_cast<int>(i / (size.x * size.y));
        const Vec3 point = geometry.gridNodeOrigin() + spacing * Vec3{Real(x), Real(y), Real(z)};
        requireNear(geometry.signedDistance(point), distances[i], "SDF query and corrected grid disagree.");
    }
}

void requireCleared(const levelset::LSInfo& geometry)
{
    require(geometry.gridNodeSFD().empty(), "Surface regeneration retained a stale grid.");
    const int3 size = geometry.gridNodeSize3D();
    require(size.x == 0 && size.y == 0 && size.z == 0, "Invalidated grid retained dimensions.");
    requireNear(geometry.gridNodeOrigin(), Vec3::zero(), "Invalidated grid retained its origin.");
    requireNear(geometry.gridNodeSpacing(), 0.0, "Invalidated grid retained spacing.");
    requireNear(geometry.boundingRadius(), 0.0, "Invalidated grid retained bounding radius.");
    requireNear(geometry.volume(), 0.0, "Invalidated grid retained volume.");
    requireNear(geometry.unitDensityInertiaTensor(), Mat3::zero(), "Invalidated grid retained inertia.");
}

void requireCopied(const LSGeometryContainer& container, int index, const levelset::LSInfo& geometry)
{
    const auto& descriptor = container.hostDescriptors().at(index);
    requireNear(descriptor.gridNodeOrigin_, geometry.gridNodeOrigin(), "Container shifted the grid again.", 0.0);
    requireNear(descriptor.gridNodeInverseSpacing_, 1.0 / geometry.gridNodeSpacing(), "Container changed spacing.", 0.0);
    const int3 size = geometry.gridNodeSize3D();
    require(descriptor.gridNodeSize_.x == size.x && descriptor.gridNodeSize_.y == size.y && descriptor.gridNodeSize_.z == size.z,
            "Container changed grid dimensions.");
    requireNear(descriptor.boundingRadius_, geometry.boundingRadius(), "Container recomputed bounding radius.", 0.0);
    requireNear(descriptor.volume_, geometry.volume(), "Container recomputed volume.", 0.0);
    requireNear(descriptor.unitDensityInertiaTensor_, geometry.unitDensityInertiaTensor(), "Container recomputed inertia.", 0.0);
    for (std::size_t i = 0; i < geometry.gridNodeSFD().size(); ++i)
        requireNear(container.hostGridNodes().at(descriptor.signedDistanceOffset_ + i).signedDistance_,
                    geometry.gridNodeSFD()[i], "Container changed a grid sample.", 0.0);
    for (std::size_t i = 0; i < geometry.surfaceNodePosition().size(); ++i)
    {
        const auto& node = container.hostSurfaceNodes().at(descriptor.surfaceNodeOffset_ + i);
        requireNear(node.position_, geometry.surfaceNodePosition()[i], "Container shifted the surface again.", 0.0);
        require(std::isfinite(node.area_) && node.area_ > 0.0, "Container did not compute valid nodal areas.");
    }
    for (std::size_t i = 0; i < geometry.surfaceNodeConnectivity().size(); ++i)
    {
        const int3 actual = container.hostSurfaceTriangles().at(descriptor.surfaceTriangleOffset_ + i).connectivity_;
        const int3 expected = geometry.surfaceNodeConnectivity()[i];
        require(actual.x == expected.x && actual.y == expected.y && actual.z == expected.z,
                "Container changed local triangle connectivity.");
    }
}

void testCentroidAndMetadata()
{
    const Vec3 offset{4.0, -3.0, 2.0};
    auto centered = translatedBox(Vec3::zero());
    auto translated = translatedBox(offset);
    requireNear(translated.signedDistance(offset), -0.5, "Unbuilt mesh lost its authored frame.");
    centered.buildLSGrid(0.125, 2);
    translated.buildLSGrid(0.125, 2);
    GeometryState(centered).requireMatches(translated);
    requireNear(translated.boundingRadius(), std::sqrt(1.5), "Incorrect centroidal bounding radius.");
    requireNear(translated.volume(), 2.0, "Box volume is not consistent with unit density.", 0.3);
    requireNear(translated.unitDensityInertiaTensor(), Mat3::diagonal({1.0 / 3.0, 5.0 / 6.0, 5.0 / 6.0}),
                "Box centroidal inertia is not consistent with unit density.", 0.15);
    requireNear(translated.signedDistance(Vec3::zero()), -0.5, "Mesh query did not follow the corrected frame.");
    requireNear(translated.signedDistance({1.25, 0.0, 0.0}), 0.25, "Corrected exterior mesh query is wrong.");
    requireQueryConsistency(translated);
    const GeometryState firstBuild(translated);
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        translated.buildLSGrid(0.125, 2);
        firstBuild.requireMatches(translated);
    }

    // Sign reversal changes stored samples only. Insertion must copy cached properties,
    // not integrate those changed samples or apply another centroid correction.
    translated.reverseSDFSign();
    LSGeometryContainer container;
    require(container.add(translated) == 0 && container.add(translated) == 1, "Geometry append indices are unstable.");
    requireCopied(container, 0, translated);
    requireCopied(container, 1, translated);
    translated.reverseSDFSign();
    firstBuild.requireMatches(translated);

    translated.fineMesh();
    requireCleared(translated);
    auto refinedReference = translatedBox(offset);
    refinedReference.fineMesh();
    require(translated.surfaceNodeConnectivity().size() == 48, "Mesh refinement did not subdivide all triangles.");
    requireSameSurface(translated.surfaceNodePosition(), refinedReference.surfaceNodePosition());
    translated.buildLSGrid(0.125, 2);
    refinedReference.buildLSGrid(0.125, 2);
    GeometryState(refinedReference).requireMatches(translated);
    requireQueryConsistency(translated);
}

void requireFixed(const levelset::LSInfo& geometry, const std::vector<Vec3>& authoredSurface)
{
    require(!geometry.gridNodeSFD().empty(), "Fixed geometry did not build a grid.");
    requireSameSurface(geometry.surfaceNodePosition(), authoredSurface);
    requireNear(geometry.volume(), 0.0, "Fixed geometry was volume-integrated.", 0.0);
    requireNear(geometry.unitDensityInertiaTensor(), Mat3::zero(), "Fixed geometry was inertia-integrated.", 0.0);
    Real radius = 0.0;
    for (const Vec3& point : authoredSurface) radius = std::max(radius, math::norm(point));
    requireNear(geometry.boundingRadius(), radius, "Fixed geometry radius is not in its authored frame.");
    requireQueryConsistency(geometry);
}

void testFixedGeometryAndShapeDefaults()
{
    const Vec3 offset{4.0, -3.0, 2.0};
    auto mesh = translatedBox(offset);
    const auto authoredSurface = mesh.surfaceNodePosition();
    mesh.buildLSGrid(0.125, 2, true);
    requireFixed(mesh, authoredSurface);
    requireNear(mesh.gridNodeOrigin(), offset + Vec3{-1.25, -0.75, -0.75}, "Fixed grid origin was recentered.");
    LSGeometryContainer container;
    requireCopied(container, container.add(mesh), mesh);
    mesh.buildLSGrid(16, false);
    require(mesh.volume() > 0.0, "Explicit movable mesh failed to integrate.");
    mesh.buildLSGrid(16, true);
    requireFixed(mesh, authoredSurface);

    levelset::PlaneWall plane(Vec3::unitZ(), 1.0);
    levelset::BoxWall box({1.0, 1.0, 1.0});
    levelset::CylinderWall cylinder({3.0, -2.0, -0.5}, {3.0, -2.0, 0.5}, 0.5);
    levelset::ConeWall cone({3.0, -2.0, -0.5}, {3.0, -2.0, 0.5}, 0.5, 0.25);
    cylinder.setCircumferentialSegments(12);
    cone.setCircumferentialSegments(12);
    for (levelset::LSInfo* wall : std::vector<levelset::LSInfo*>{&plane, &box, &cylinder, &cone})
    {
        const auto originalSurface = wall->surfaceNodePosition();
        wall->buildLSGrid(0.25, 2);
        requireFixed(*wall, originalSurface);
        wall->buildLSGrid(8);
        requireFixed(*wall, originalSurface);
        wall->buildLSGrid(0.25, 2, false);
        require(wall->volume() > 0.0, "Explicit movable wall did not override its fixed default.");
        wall->buildLSGrid(0.25, 2);
        requireFixed(*wall, originalSurface);
    }

    levelset::BoxParticle particle({1.0, 1.0, 1.0});
    particle.buildSurfaceNode(0.5);
    particle.buildLSGrid(0.25, 2);
    require(particle.volume() > 0.0, "BoxParticle inherited the wall's fixed default.");
    levelset::BoxWall& asWall = particle;
    asWall.buildLSGrid(8);
    require(particle.volume() > 0.0, "BoxParticle default changed through a BoxWall reference.");
    levelset::LSInfo& asInfo = particle;
    asInfo.buildLSGrid(0.25, 2);
    require(particle.volume() > 0.0, "BoxParticle default changed through an LSInfo reference.");
}

void testRegenerationAndRandomShape()
{
    auto random = levelset::makeRandomShape(1.0, -0.15, 0.2, 1, 2026);
    random->buildLSGrid(12);
    const GeometryState firstRandomBuild(*random);
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        random->buildLSGrid(12);
        firstRandomBuild.requireMatches(*random);
        requireQueryConsistency(*random);
    }

    levelset::Sphere sphere(1.0);
    sphere.buildSurfaceNode(0);
    sphere.buildLSGrid(0.25, 2);
    sphere.buildSurfaceNode(1);
    requireCleared(sphere);
    sphere.buildLSGrid(0.25, 2);
    requireQueryConsistency(sphere);
    sphere.setParameter(0.75);
    requireCleared(sphere);

    levelset::Superellipsoid ellipsoid(1.0, 0.75, 0.5, 1.0, 1.0);
    ellipsoid.buildSurfaceNode(0);
    ellipsoid.buildLSGrid(0.25, 2);
    ellipsoid.buildSurfaceNode(1);
    requireCleared(ellipsoid);
    ellipsoid.buildLSGrid(0.25, 2);
    requireQueryConsistency(ellipsoid);

    levelset::BoxParticle box({1.0, 1.0, 1.0});
    box.buildSurfaceNode(0.5);
    box.buildLSGrid(0.25, 2);
    box.buildSurfaceNode(0.25);
    requireCleared(box);
    box.buildLSGrid(0.25, 2);
    requireQueryConsistency(box);

    const Vec3 bottom{3.0, -2.0, -0.5};
    const Vec3 top{3.0, -2.0, 0.5};
    levelset::CylinderWall cylinder(bottom, top, 0.5);
    cylinder.setCircumferentialSegments(8);
    cylinder.buildLSGrid(0.25, 2, false);
    cylinder.setCircumferentialSegments(12);
    requireCleared(cylinder);
    levelset::CylinderWall freshCylinder(bottom, top, 0.5);
    freshCylinder.setCircumferentialSegments(12);
    requireSameSurface(cylinder.surfaceNodePosition(), freshCylinder.surfaceNodePosition());
    cylinder.buildLSGrid(0.25, 2);
    freshCylinder.buildLSGrid(0.25, 2);
    GeometryState(freshCylinder).requireMatches(cylinder);

    levelset::ConeWall cone(bottom, top, 0.5, 0.25);
    cone.setCircumferentialSegments(8);
    cone.buildLSGrid(0.25, 2, false);
    cone.setCircumferentialSegments(12);
    requireCleared(cone);
    levelset::ConeWall freshCone(bottom, top, 0.5, 0.25);
    freshCone.setCircumferentialSegments(12);
    requireSameSurface(cone.surfaceNodePosition(), freshCone.surfaceNodePosition());
    cone.buildLSGrid(0.25, 2);
    freshCone.buildLSGrid(0.25, 2);
    GeometryState(freshCone).requireMatches(cone);
}
} // namespace

int main()
{
    try
    {
        testCentroidAndMetadata();
        testFixedGeometryAndShapeDefaults();
        testRegenerationAndRandomShape();
        std::cout << "LS geometry ownership and lifecycle tests passed.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "LS geometry test failed: " << error.what() << '\n';
        return 1;
    }
}
