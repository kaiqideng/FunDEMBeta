#include "data/myLSObject.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
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

void validateMesh(const levelset::TriangleMesh& mesh, math::Real minimumRadius, math::Real maximumRadius)
{
    require(mesh.valid(), "Random shape must be valid before grid construction.");
    const auto& vertices = mesh.surfaceNodePosition();
    std::map<std::pair<int, int>, std::pair<int, int>> edges;
    math::Real smallest = std::numeric_limits<math::Real>::max();
    math::Real largest = 0.0;
    for (const auto& vertex : vertices)
    {
        require(math::isFinite(vertex), "All random-shape vertices must be finite.");
        const math::Real radius = math::norm(vertex);
        smallest = std::min(smallest, radius);
        largest = std::max(largest, radius);
        require(radius >= minimumRadius - 1.0e-12 && radius <= maximumRadius + 1.0e-12, "Vertex radii must respect the supplied height bounds.");
    }
    require(std::abs(smallest - minimumRadius) < 1.0e-12 && std::abs(largest - maximumRadius) < 1.0e-12, "Sampled heights must span both requested endpoints.");
    for (const auto& face : mesh.surfaceNodeConnectivity())
    {
        const auto normal = math::cross(vertices.at(face.y) - vertices.at(face.x), vertices.at(face.z) - vertices.at(face.x));
        require(math::isFinite(normal) && math::dot(normal, vertices.at(face.x)) > 0.0, "Every triangle must be nondegenerate and oriented outward.");
        for (const auto& edge : {std::pair<int, int>{face.x, face.y}, {face.y, face.z}, {face.z, face.x}})
        {
            auto& entry = edges[std::minmax(edge.first, edge.second)];
            ++entry.first;
            entry.second += edge.first < edge.second ? 1 : -1;
        }
    }
    for (const auto& [edge, incidence] : edges)
    {
        (void)edge;
        require(incidence.first == 2 && incidence.second == 0, "Mesh edges must join two consistently oriented faces.");
    }
    require(mesh.signedDistance(math::Vec3::zero()) < 0.0, "The origin must be inside the star-shaped particle.");
    require(mesh.signedDistance({2.0 * maximumRadius, 0.0, 0.0}) > 0.0, "A point beyond the enclosing sphere must be outside.");
}

void testGeometryAndDistance()
{
    for (std::uint64_t seed = 0; seed < 4; ++seed)
    {
        const auto mesh = levelset::makeRandomShape(1.0, -0.35, 0.2, 2, seed);
        validateMesh(*mesh, 0.65, 1.2);
    }
    const auto first = levelset::makeRandomShape(1.0, -0.2, 0.3, 3, 1234);
    const auto repeat = levelset::makeRandomShape(1.0, -0.2, 0.3, 3, 1234);
    const auto different = levelset::makeRandomShape(1.0, -0.2, 0.3, 3, 1235);
    require(first->surfaceNodePosition() == repeat->surfaceNodePosition(), "A fixed seed must exactly reproduce the sampled vertices.");
    require(first->surfaceNodePosition() != different->surfaceNodePosition(), "Different seeds must produce different shapes.");
    require(first->gridNodeSFD().empty(), "The factory must not allocate an LS grid.");
    require(first->surfaceNodePosition().size() == 642 && first->surfaceNodeConnectivity().size() == 1280, "Subdivision level three must retain the expected icosphere resolution.");
    validateMesh(*first, 0.8, 1.3);

    const auto face = first->surfaceNodeConnectivity().front();
    const auto& vertices = first->surfaceNodePosition();
    const auto center = (vertices[face.x] + vertices[face.y] + vertices[face.z]) / 3.0;
    const auto normal = math::normalizedOrZero(math::cross(vertices[face.y] - vertices[face.x], vertices[face.z] - vertices[face.x]));
    constexpr math::Real offset = 1.0e-6;
    require(std::abs(first->signedDistance(center + offset * normal) - offset) < 1.0e-10, "Outside SDF must equal Euclidean distance to the mesh face.");
    require(std::abs(first->signedDistance(center - offset * normal) + offset) < 1.0e-10, "Inside SDF must equal negative Euclidean distance to the mesh face.");

    first->buildLSGrid(8);
    require(!first->gridNodeSFD().empty(), "The random mesh must support subsequent LS-grid construction.");
    require(std::all_of(first->gridNodeSFD().begin(), first->gridNodeSFD().end(), [](math::Real value) { return math::isFinite(value); }), "Random-mesh LS-grid samples must remain finite.");

    const auto sphere = levelset::makeRandomShape(1.0, 0.1, 0.1, 0);
    validateMesh(*sphere, 1.1, 1.1);
}

void expectInvalid(math::Real radius, math::Real minimum, math::Real maximum, int subdivision = 1)
{
    bool rejected = false;
    try
    {
        (void)levelset::makeRandomShape(radius, minimum, maximum, subdivision);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "Invalid random-shape parameters must fail before returning a mesh.");
}

void testValidation()
{
    expectInvalid(0.0, 0.0, 0.1);
    expectInvalid(1.0, 0.2, 0.1);
    expectInvalid(1.0, -1.0, 0.1);
    expectInvalid(1.0, 0.0, 0.1, -1);
    expectInvalid(1.0, 0.0, 0.1, std::numeric_limits<int>::max());
    expectInvalid(std::numeric_limits<math::Real>::infinity(), 0.0, 0.1);
    expectInvalid(1.0, std::numeric_limits<math::Real>::quiet_NaN(), 0.1);
    expectInvalid(1.0, 0.0, std::numeric_limits<math::Real>::infinity());
    expectInvalid(std::numeric_limits<math::Real>::max(), 0.0, std::numeric_limits<math::Real>::max());
}

void testSuperellipsoidPerturbation()
{
    const math::Vec3 semiAxes{1.4, 0.8, 1.0};
    constexpr math::Real equatorialExponent = 0.8;
    constexpr math::Real polarExponent = 1.2;
    constexpr math::Real minimumHeight = -0.12;
    constexpr math::Real maximumHeight = 0.18;
    levelset::Superellipsoid base(semiAxes.x, semiAxes.y, semiAxes.z, equatorialExponent, polarExponent);
    base.buildSurfaceNode(2);
    auto mesh = levelset::makeRandomShape(semiAxes, equatorialExponent, polarExponent, minimumHeight, maximumHeight, 2, 537);
    const auto repeat = levelset::makeRandomShape(semiAxes, equatorialExponent, polarExponent, minimumHeight, maximumHeight, 2, 537);
    const auto different = levelset::makeRandomShape(semiAxes, equatorialExponent, polarExponent, minimumHeight, maximumHeight, 2, 538);
    const auto unperturbed = levelset::makeRandomShape(semiAxes, equatorialExponent, polarExponent, 0.0, 0.0, 2, 537);
    require(mesh->surfaceNodePosition() == repeat->surfaceNodePosition(), "Superellipsoid perturbations must be reproducible from their seed.");
    require(mesh->surfaceNodePosition() != different->surfaceNodePosition(), "Superellipsoid perturbations must depend on their seed.");
    require(mesh->surfaceNodePosition().size() == base.surfaceNodePosition().size(), "Random geometry must retain the requested base subdivision.");
    math::Real lowestHeight = std::numeric_limits<math::Real>::max();
    math::Real highestHeight = std::numeric_limits<math::Real>::lowest();
    for (std::size_t index = 0; index < mesh->surfaceNodePosition().size(); ++index)
    {
        const auto& original = base.surfaceNodePosition()[index];
        const auto& deformed = mesh->surfaceNodePosition()[index];
        const math::Real height = math::norm(deformed) - math::norm(original);
        lowestHeight = std::min(lowestHeight, height);
        highestHeight = std::max(highestHeight, height);
        require(height >= minimumHeight - 1.0e-12 && height <= maximumHeight + 1.0e-12, "Perturbations must stay within radial height bounds relative to the superellipsoid.");
        require(math::norm(math::normalizedOrZero(original) - math::normalizedOrZero(deformed)) < 1.0e-12, "Superellipsoid perturbations must be radial.");
        require(math::norm(unperturbed->surfaceNodePosition()[index] - original) < 1.0e-12, "Zero height must preserve the superellipsoid base, including its exponents.");
    }
    require(std::abs(lowestHeight - minimumHeight) < 1.0e-12 && std::abs(highestHeight - maximumHeight) < 1.0e-12, "Superellipsoid offsets must span both requested endpoints.");
    for (const auto& face : mesh->surfaceNodeConnectivity())
    {
        const auto& vertices = mesh->surfaceNodePosition();
        const auto normal = math::cross(vertices[face.y] - vertices[face.x], vertices[face.z] - vertices[face.x]);
        require(math::dot(normal, vertices[face.x]) > 0.0, "Superellipsoid deformation must preserve outward winding.");
    }
    const auto nativePositions = mesh->surfaceNodePosition();
    mesh->buildLSGrid(0.16, 2);
    require(mesh->volume() > 0.0 && mesh->boundingRadius() > 0.0, "Irregular superellipsoids must have integrated volume and a bounding radius.");
    const auto correctedPositions = mesh->surfaceNodePosition();
    require(math::norm(nativePositions.front() - correctedPositions.front()) > 1.0e-5, "Asymmetric shapes must use LSInfo centroid correction.");
    mesh->buildLSGrid(0.16, 2);
    for (std::size_t index = 0; index < correctedPositions.size(); ++index)
        require(math::norm(correctedPositions[index] - mesh->surfaceNodePosition()[index]) < 1.0e-12, "Rebuilding must not accumulate centroid corrections.");

    for (const auto& invalid : {math::Vec3{0.0, 1.0, 1.0}, math::Vec3{1.0, -1.0, 1.0}})
    {
        bool rejected = false;
        try { (void)levelset::makeRandomShape(invalid, 1.0, 1.0, -0.1, 0.1, 1); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Invalid superellipsoid semi-axes must be rejected.");
    }
    bool rejected = false;
    try { (void)levelset::makeRandomShape(semiAxes, 1.0, 1.0, -2.0, -2.0, 1); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "Offsets that collapse the base must be rejected.");
}

} // namespace

int main()
{
    try
    {
        testGeometryAndDistance();
        testValidation();
        testSuperellipsoidPerturbation();
        std::cout << "Random-shape geometry tests passed.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
