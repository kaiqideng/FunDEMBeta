#include "data/myLSObject.h"
#include "geometry/LSGeometry.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
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

template <class Operation> void requireLogicError(Operation operation, const char* message)
{
    bool rejected = false;
    try { operation(); }
    catch (const std::logic_error&) { rejected = true; }
    require(rejected, message);
}

template <class Operation> void requireRejected(Operation operation, const char* message)
{
    bool rejected = false;
    try { operation(); }
    catch (const std::exception&) { rejected = true; }
    require(rejected, message);
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

void requireSameSize(const int3& actual, const int3& expected, const char* message)
{
    require(actual.x == expected.x && actual.y == expected.y && actual.z == expected.z, message);
}

void requireSurfaceBounds(const levelset::LSInfo& geometry)
{
    const auto& surface = geometry.surfaceNodePosition();
    require(!surface.empty(), "Surface bounds fixture has no surface nodes.");
    Vec3 minimum = surface.front();
    Vec3 maximum = surface.front();
    for (const Vec3& point : surface)
        for (int axis = 0; axis < 3; ++axis)
        {
            minimum[axis] = std::min(minimum[axis], point[axis]);
            maximum[axis] = std::max(maximum[axis], point[axis]);
        }
    requireNear(geometry.surfaceAABBMin(), minimum, "Public lower bound is not the current surface minimum.", 0.0);
    requireNear(geometry.surfaceAABBMax(), maximum, "Public upper bound is not the current surface maximum.", 0.0);
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

void requireSurfaceInsideGrid(const levelset::LSInfo& geometry)
{
    requireSurfaceBounds(geometry);
    const Vec3 gridMinimum = geometry.gridNodeOrigin();
    const int3 size = geometry.gridNodeSize3D();
    const Vec3 gridMaximum = gridMinimum + geometry.gridNodeSpacing() * Vec3{Real(size.x - 1), Real(size.y - 1), Real(size.z - 1)};
    for (const Vec3& point : geometry.surfaceNodePosition())
        for (int axis = 0; axis < 3; ++axis)
        {
            require(point[axis] >= gridMinimum[axis] - 1.0e-12 && point[axis] <= gridMaximum[axis] + 1.0e-12,
                    "A plane corner lies outside its padded SDF grid.");
        }
}

void testPlaneSamplingBounds()
{
    for (const Vec3& normal : std::vector<Vec3>{Vec3::unitX(), -Vec3::unitX(), Vec3::unitY(), -Vec3::unitY(), Vec3::unitZ(), -Vec3::unitZ()})
    {
        levelset::PlaneWall axisPlane(normal, 1.0);
        Vec3 halfExtent{0.5};
        for (int axis = 0; axis < 3; ++axis)
            if (normal[axis] != 0.0) halfExtent[axis] = 0.0;
        requireNear(axisPlane.surfaceAABBMin(), -halfExtent, "Axis-plane surface bounds include artificial sampling thickness.", 0.0);
        requireNear(axisPlane.surfaceAABBMax(), halfExtent, "Axis-plane surface bounds include artificial sampling thickness.", 0.0);
        const int3 estimate = axisPlane.estimateGridNodeSize(0.125, 2);
        requireSameSize(estimate, {13, 13, 13}, "Axis-plane grid estimate lost its original sampling cube.");
        axisPlane.buildLSGrid(0.125, 2);
        requireSameSize(axisPlane.gridNodeSize3D(), estimate, "Axis-plane estimate differs from the actual grid.");
        requireSurfaceInsideGrid(axisPlane);
    }

    levelset::PlaneWall plane({1.0, 1.0, 1.0}, 1.0);
    const Vec3 nativeMinimum = plane.surfaceAABBMin();
    const Vec3 nativeMaximum = plane.surfaceAABBMax();
    const auto nativeSurface = plane.surfaceNodePosition();
    const int3 estimate = plane.estimateGridNodeSize(0.1, 2);
    plane.buildLSGrid(100);
    requireFixed(plane, nativeSurface);
    requireSurfaceInsideGrid(plane);
    requireNear(plane.surfaceAABBMin(), nativeMinimum, "Fixed plane shifted its lower surface bound.");
    requireNear(plane.surfaceAABBMax(), nativeMaximum, "Fixed plane shifted its upper surface bound.");

    plane.buildLSGrid(0.1, 2, false);
    requireSameSize(plane.gridNodeSize3D(), estimate, "Movable-plane estimate differs from its actual grid.");
    const Vec3 correction = nativeSurface.front() - plane.surfaceNodePosition().front();
    require(math::norm(correction) > 0.01, "Movable oblique-plane fixture did not exercise centroid correction.");
    requireNear(plane.surfaceAABBMin(), nativeMinimum - correction, "Movable-plane lower bound has an inconsistent local frame.");
    requireNear(plane.surfaceAABBMax(), nativeMaximum - correction, "Movable-plane upper bound has an inconsistent local frame.");
    requireSurfaceInsideGrid(plane);
    const GeometryState firstBuild(plane);
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        requireSameSize(plane.estimateGridNodeSize(0.1, 2), estimate, "Centroid correction changed the native grid dimensions.");
        firstBuild.requireMatches(plane);
        plane.buildLSGrid(0.1, 2, false);
        firstBuild.requireMatches(plane);
        requireNear(plane.surfaceAABBMin(), nativeMinimum - correction, "Rebuilding applied the plane's centroid offset twice.");
        requireNear(plane.surfaceAABBMax(), nativeMaximum - correction, "Rebuilding applied the plane's centroid offset twice.");
        requireSurfaceInsideGrid(plane);
        requireQueryConsistency(plane);
    }
    plane.buildLSGrid(0.1, 2);
    requireFixed(plane, nativeSurface);
    requireNear(plane.surfaceAABBMin(), nativeMinimum, "Returning to fixed policy retained a shifted lower bound.");
    requireNear(plane.surfaceAABBMax(), nativeMaximum, "Returning to fixed policy retained a shifted upper bound.");
}

void requireUnavailableMesh(levelset::TriangleMesh& geometry)
{
    require(!geometry.valid(), "An empty or moved-from mesh remained configured.");
    requireLogicError([&] { (void)geometry.surfaceAABBMin(); }, "An unavailable mesh returned a lower surface bound.");
    requireLogicError([&] { (void)geometry.surfaceAABBMax(); }, "An unavailable mesh returned an upper surface bound.");
    requireLogicError([&] { (void)geometry.estimateGridNodeSize(0.125, 2); }, "An unavailable mesh accepted a grid estimate.");
    requireLogicError([&] { geometry.buildLSGrid(16); }, "An unavailable mesh accepted a resolution-based grid build.");
    requireLogicError([&] { geometry.buildLSGrid(0.125, 2); }, "An unavailable mesh accepted a spacing-based grid build.");
    requireLogicError([&] { (void)geometry.signedDistance(Vec3::zero()); }, "An unavailable mesh accepted a signed-distance query.");
}

void testUnavailableMeshBounds()
{
    levelset::TriangleMesh empty;
    requireUnavailableMesh(empty);
    auto original = translatedBox({4.0, -3.0, 2.0});
    original.buildLSGrid(0.125, 2);
    const GeometryState accepted(original);
    const Vec3 minimum = original.surfaceAABBMin();
    const Vec3 maximum = original.surfaceAABBMax();

    auto moved = std::move(original);
    requireUnavailableMesh(original);
    require(moved.valid(), "Move construction invalidated the destination mesh.");
    accepted.requireMatches(moved);
    requireNear(moved.surfaceAABBMin(), minimum, "Move construction changed the prepared lower frame.");
    requireNear(moved.surfaceAABBMax(), maximum, "Move construction changed the prepared upper frame.");
    requireQueryConsistency(moved);

    levelset::TriangleMesh assigned;
    assigned = std::move(moved);
    requireUnavailableMesh(moved);
    require(assigned.valid(), "Move assignment invalidated the destination mesh.");
    accepted.requireMatches(assigned);
    requireNear(assigned.surfaceAABBMin(), minimum, "Move assignment changed the prepared lower frame.");
    requireNear(assigned.surfaceAABBMax(), maximum, "Move assignment changed the prepared upper frame.");
    requireQueryConsistency(assigned);
}

void requirePreparedBounds(const levelset::LSInfo& source)
{
    auto prepared = levelset::LSInfo::fromPreparedData(source.gridNodeOrigin(), source.gridNodeSpacing(), source.gridNodeSize3D(), source.gridNodeSFD(),
                                                     source.surfaceNodePosition(), source.surfaceNodeConnectivity(), source.boundingRadius(), source.volume(),
                                                     source.unitDensityInertiaTensor());
    const GeometryState accepted(source);
    const int3 size = source.gridNodeSize3D();
    const Vec3 lastNode = source.gridNodeOrigin() + source.gridNodeSpacing() * Vec3{Real(size.x - 1), Real(size.y - 1), Real(size.z - 1)};
    requireSurfaceBounds(*prepared);
    requireNear(prepared->surfaceAABBMin(), source.surfaceAABBMin(), "Prepared lower bound changed the accepted surface frame.", 0.0);
    requireNear(prepared->surfaceAABBMax(), source.surfaceAABBMax(), "Prepared upper bound changed the accepted surface frame.", 0.0);
    require(math::norm(prepared->surfaceAABBMin() - source.gridNodeOrigin()) > 0.0 &&
            math::norm(prepared->surfaceAABBMax() - lastNode) > 0.0,
            "Prepared surface bounds were replaced by padded grid bounds.");
    require(prepared->gridNodeSFD() == source.gridNodeSFD() && prepared->surfaceNodePosition() == source.surfaceNodePosition(),
            "Restoring prepared geometry resampled or recentered its data.");
    require(prepared->surfaceNodeConnectivity().size() == source.surfaceNodeConnectivity().size(), "Prepared triangle count changed.");
    for (std::size_t index = 0; index < source.surfaceNodeConnectivity().size(); ++index)
    {
        const int3 actual = prepared->surfaceNodeConnectivity()[index];
        const int3 expected = source.surfaceNodeConnectivity()[index];
        require(actual.x == expected.x && actual.y == expected.y && actual.z == expected.z, "Prepared triangle connectivity changed.");
    }
    requireNear(prepared->boundingRadius(), source.boundingRadius(), "Prepared radius changed.", 0.0);
    requireNear(prepared->volume(), source.volume(), "Prepared volume changed.", 0.0);
    requireNear(prepared->unitDensityInertiaTensor(), source.unitDensityInertiaTensor(), "Prepared inertia changed.", 0.0);
    accepted.requireMatches(*prepared);
    requireQueryConsistency(*prepared);
    requireNear(prepared->signedDistance(source.gridNodeOrigin()), source.gridNodeSFD().front(), "Prepared lower boundary is not queryable.");
    requireNear(prepared->signedDistance(lastNode), source.gridNodeSFD().back(), "Prepared upper boundary is not queryable.");
    requireLogicError([&] { prepared->buildLSGrid(16); }, "Prepared geometry accepted a source-shape rebuild.");
    requireLogicError([&] { (void)prepared->estimateGridNodeSize(source.gridNodeSpacing(), 2); },
                      "Prepared geometry accepted a source-shape estimate.");
    accepted.requireMatches(*prepared);
}

void testGridSizeEstimation()
{
    auto mesh = translatedBox({4.0, -3.0, 2.0});
    const GeometryState native(mesh);
    const int3 estimate = mesh.estimateGridNodeSize(0.125, 2);
    requireSameSize(estimate, {21, 13, 13}, "Grid estimate changed the existing ceil/padding scheme.");
    requireSameSize(mesh.estimateGridNodeSize(0.125, 0), {17, 9, 9}, "Zero-padding estimate has incorrect dimensions.");
    native.requireMatches(mesh);
    requireSurfaceBounds(mesh);
    mesh.buildLSGrid(0.125, 2);
    requireSameSize(mesh.gridNodeSize3D(), estimate, "Translated mesh estimate differs from its actual grid.");
    requireSurfaceBounds(mesh);
    const GeometryState accepted(mesh);
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        requireSameSize(mesh.estimateGridNodeSize(0.125, 2), estimate, "Recentered mesh estimate changed dimensions.");
        accepted.requireMatches(mesh);
        mesh.buildLSGrid(0.125, 2);
        accepted.requireMatches(mesh);
        requireSurfaceBounds(mesh);
    }
    for (Real spacing : {0.0, -0.125, std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN(),
                         std::numeric_limits<Real>::denorm_min(), 1.0e-20, 0.0005})
    {
        requireRejected([&] { (void)mesh.estimateGridNodeSize(spacing, 2); }, "Invalid or overflowing grid estimate was accepted.");
        accepted.requireMatches(mesh);
    }
    for (int padding : {-1, std::numeric_limits<int>::max()})
    {
        requireRejected([&] { (void)mesh.estimateGridNodeSize(0.125, padding); }, "Invalid or overflowing grid padding was accepted.");
        accepted.requireMatches(mesh);
    }

    // Analytic parameters suffice for sizing, but cannot substitute for a missing discrete surface AABB.
    levelset::Sphere sphere(1.0);
    require(sphere.valid() && sphere.surfaceNodePosition().empty(), "Analytic no-surface fixture is invalid.");
    requireLogicError([&] { (void)sphere.surfaceAABBMin(); }, "A shape without surface nodes invented a lower surface bound.");
    requireLogicError([&] { (void)sphere.surfaceAABBMax(); }, "A shape without surface nodes invented an upper surface bound.");
    requireSameSize(sphere.estimateGridNodeSize(0.25, 2), {13, 13, 13}, "Analytic no-surface estimate is incorrect.");
    require(sphere.surfaceNodePosition().empty() && sphere.gridNodeSFD().empty(), "Grid estimation mutated an unbuilt analytic shape.");

    // Public surface bounds include every retained entry; querying must not prune or remap supplied vertices.
    auto authored = translatedBox(Vec3::zero());
    auto points = authored.surfaceNodePosition();
    points.push_back({10.0, 20.0, 30.0});
    const auto faces = authored.surfaceNodeConnectivity();
    levelset::TriangleMesh retained(points, faces);
    requireSurfaceBounds(retained);
    requireNear(retained.surfaceAABBMax(), points.back(), "Surface bounds ignored a retained surface entry.", 0.0);
    require(retained.surfaceNodePosition() == points && retained.surfaceNodeConnectivity().size() == faces.size(),
            "Surface bounds query modified supplied mesh data.");
}

void testPreparedGeometryBounds()
{
    auto corrected = translatedBox({4.0, -3.0, 2.0});
    corrected.buildLSGrid(0.125, 2);
    requirePreparedBounds(corrected);
    levelset::PlaneWall plane({1.0, 1.0, 1.0}, 1.0);
    plane.buildLSGrid(0.1, 2);
    requirePreparedBounds(plane);
}
} // namespace

int main()
{
    try
    {
        testCentroidAndMetadata();
        testFixedGeometryAndShapeDefaults();
        testRegenerationAndRandomShape();
        testPlaneSamplingBounds();
        testUnavailableMeshBounds();
        testPreparedGeometryBounds();
        testGridSizeEstimation();
        std::cout << "LS geometry ownership and lifecycle tests passed.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "LS geometry test failed: " << error.what() << '\n';
        return 1;
    }
}
