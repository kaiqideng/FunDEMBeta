#include "myLSObject.h"
#include "math/Indexing.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace fundem::levelset
{
namespace
{

using Real = math::Real;
using Vec3 = math::Vec3;

int3 triangle(int first, int second, int third) noexcept { return {first, second, third}; }

std::uint64_t edgeKey(int first, int second) noexcept
{
    const std::uint32_t low = static_cast<std::uint32_t>(std::min(first, second));
    const std::uint32_t high = static_cast<std::uint32_t>(std::max(first, second));
    return (static_cast<std::uint64_t>(low) << 32U) | high;
}

void orthogonalBasis(const Vec3& normal, Vec3& tangent1, Vec3& tangent2) noexcept
{
    tangent1 = std::abs(normal.x) > std::abs(normal.z) ? Vec3{-normal.y, normal.x, 0.0} : Vec3{0.0, -normal.z, normal.y};
    tangent1 = math::normalizedOrZero(tangent1);
    tangent2 = math::cross(normal, tangent1);
}

bool finitePositive(Real value) noexcept { return math::isFinite(value) && value > 0.0; }

} // namespace

void LSInfo::clearGrid() noexcept
{
    if (centroidOffset_ != Vec3::zero())
        for (Vec3& point : surfaceNodePositions_)
            point += centroidOffset_;
    centroidOffset_ = Vec3::zero();
    boundingRadius_ = 0.0;
    volume_ = 0.0;
    unitDensityInertiaTensor_ = Mat3::zero();
    gridNodeOrigin_ = Vec3::zero();
    gridNodeSize_ = {0, 0, 0};
    gridNodeSpacing_ = 0.0;
    gridNodeSignedDistance_.clear();
}

void LSInfo::buildLSGrid(int resolutionPerDiameter) { buildLSGrid(resolutionPerDiameter, defaultFixedGeometry()); }

void LSInfo::buildLSGrid(Real spacing, int paddingSize) { buildLSGrid(spacing, paddingSize, defaultFixedGeometry()); }

void LSInfo::buildLSGrid(int resolutionPerDiameter, bool isFixed)
{
    if (!isValid())
    {
        throw std::logic_error("Configure the level-set geometry before building its grid.");
    }
    if (resolutionPerDiameter <= 0)
    {
        throw std::invalid_argument("Level-set resolution must be positive.");
    }

    const Vec3 extent = boundingBoxMax() - boundingBoxMin();
    const Real referenceDiameter = std::max({extent.x, extent.y, extent.z});
    if (!finitePositive(referenceDiameter))
    {
        throw std::domain_error("Level-set bounding box must have a positive finite extent.");
    }
    buildLSGrid(referenceDiameter / resolutionPerDiameter, 2, isFixed);
}

void LSInfo::buildLSGrid(Real spacing, int paddingSize, bool isFixed)
{
    if (!isValid())
    {
        throw std::logic_error("Configure the level-set geometry before building its grid.");
    }
    if (!finitePositive(spacing) || paddingSize < 0)
    {
        throw std::invalid_argument("Level-set spacing must be positive and padding must be non-negative.");
    }
    clearGrid();
    gridNodeSpacing_ = spacing;
    try
    {
        buildLSGridKernel(paddingSize);
        updateGeometryProperties(isFixed);
    }
    catch (...)
    {
        clearGrid();
        throw;
    }
}

void LSInfo::buildLSGridKernel(int paddingSize)
{
    const Vec3 boxMin = boundingBoxMin();
    const Vec3 boxMax = boundingBoxMax();
    if (!math::isFinite(boxMin) || !math::isFinite(boxMax) || boxMax.x <= boxMin.x || boxMax.y <= boxMin.y || boxMax.z <= boxMin.z)
    {
        throw std::domain_error("Level-set bounding box is invalid.");
    }

    const Real padding = paddingSize * gridNodeSpacing_;
    const Vec3 center = 0.5 * (boxMin + boxMax);
    const Vec3 halfExtent = 0.5 * (boxMax - boxMin) + Vec3{padding};
    const int halfSizeX = static_cast<int>(std::ceil(halfExtent.x / gridNodeSpacing_));
    const int halfSizeY = static_cast<int>(std::ceil(halfExtent.y / gridNodeSpacing_));
    const int halfSizeZ = static_cast<int>(std::ceil(halfExtent.z / gridNodeSpacing_));

    gridNodeSize_ = {2 * halfSizeX + 1, 2 * halfSizeY + 1, 2 * halfSizeZ + 1};
    gridNodeOrigin_ = center - gridNodeSpacing_ * Vec3{Real(halfSizeX), Real(halfSizeY), Real(halfSizeZ)};
    const int gridNodeCount = gridNodeSize_.x * gridNodeSize_.y * gridNodeSize_.z;
    gridNodeSignedDistance_.resize(gridNodeCount);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (gridNodeCount >= 4096)
#endif
    for (int index = 0; index < gridNodeCount; ++index)
    {
        const int planeIndex = index / gridNodeSize_.x;
        const int x = index % gridNodeSize_.x;
        const int y = planeIndex % gridNodeSize_.y;
        const int z = planeIndex / gridNodeSize_.y;
        const Vec3 point = gridNodeOrigin_ + gridNodeSpacing_ * Vec3{Real(x), Real(y), Real(z)};
        gridNodeSignedDistance_[index] = evaluateSFD(point);
    }

    if (std::any_of(gridNodeSignedDistance_.begin(), gridNodeSignedDistance_.end(), [](Real value) { return !math::isFinite(value); }))
    {
        clearGrid();
        throw std::domain_error("Level-set evaluation produced a non-finite signed distance.");
    }
}

void LSInfo::updateGeometryProperties(bool isFixed)
{
    const Real inverseSpacing = 1.0 / gridNodeSpacing_;
    const Real nodeVolume = gridNodeSpacing_ * gridNodeSpacing_ * gridNodeSpacing_;
    Real volume = 0.0;
    Vec3 centroid = Vec3::zero();
    Mat3 inertia = Mat3::zero();
    if (!isFixed)
    {
        const auto interiorFraction = [](Real distance) noexcept
        {
            constexpr Real smoothingWidth = 1.5;
            if (distance < -smoothingWidth)
                return Real{1.0};
            if (distance > smoothingWidth)
                return Real{0.0};
            const Real normalizedDistance = -distance / smoothingWidth;
            return 0.5 * (1.0 + normalizedDistance + std::sin(math::pi * normalizedDistance) / math::pi);
        };
        const int sizeX = gridNodeSize_.x;
        const int sizeY = gridNodeSize_.y;
        const int sizeZ = gridNodeSize_.z;
        Real occupancySum = 0.0;
        Vec3 firstMoment = Vec3::zero();
        for (int z = 0; z < sizeZ; ++z)
            for (int y = 0; y < sizeY; ++y)
                for (int x = 0; x < sizeX; ++x)
                {
                    const int index = math::linearIndex(x, y, z, sizeX, sizeY);
                    const Real occupancy = interiorFraction(gridNodeSignedDistance_[index] * inverseSpacing);
                    const Vec3 nodePosition = gridNodeOrigin_ + gridNodeSpacing_ * Vec3{Real(x), Real(y), Real(z)};
                    occupancySum += occupancy;
                    firstMoment += occupancy * nodePosition;
                }
        volume = occupancySum * nodeVolume;
        if (!math::isFinite(occupancySum) || occupancySum <= math::defaultTolerance || !math::isFinite(volume) || volume <= math::defaultTolerance)
            throw std::domain_error("Level-set grid has no finite interior volume.");
        centroid = firstMoment / occupancySum;
        if (!math::isFinite(centroid))
            throw std::overflow_error("Level-set centroid integration is not finite.");
        for (int z = 0; z < sizeZ; ++z)
            for (int y = 0; y < sizeY; ++y)
                for (int x = 0; x < sizeX; ++x)
                {
                    const int index = math::linearIndex(x, y, z, sizeX, sizeY);
                    const Real pointMass = interiorFraction(gridNodeSignedDistance_[index] * inverseSpacing) * nodeVolume;
                    const Vec3 position = gridNodeOrigin_ + gridNodeSpacing_ * Vec3{Real(x), Real(y), Real(z)} - centroid;
                    inertia(0, 0) += pointMass * (position.y * position.y + position.z * position.z);
                    inertia(1, 1) += pointMass * (position.x * position.x + position.z * position.z);
                    inertia(2, 2) += pointMass * (position.x * position.x + position.y * position.y);
                    inertia(0, 1) -= pointMass * position.x * position.y;
                    inertia(0, 2) -= pointMass * position.x * position.z;
                    inertia(1, 2) -= pointMass * position.y * position.z;
                }
        inertia(1, 0) = inertia(0, 1);
        inertia(2, 0) = inertia(0, 2);
        inertia(2, 1) = inertia(1, 2);
        if (!math::isFinite(inertia))
            throw std::overflow_error("Level-set inertia integration is not finite.");
    }
    Real radius = 0.0;
    const auto includeRadius = [&](const Vec3& point)
    {
        const Real candidate = math::norm(point - centroid);
        if (!math::isFinite(candidate))
            throw std::overflow_error("Level-set bounding radius is not finite.");
        radius = std::max(radius, candidate);
    };
    for (const Vec3& point : surfaceNodePositions_)
        includeRadius(point);
    if (surfaceNodePositions_.empty())
    {
        // Grid-only analytic queries remain supported before surface sampling.
        const Vec3 minimum = boundingBoxMin();
        const Vec3 maximum = boundingBoxMax();
        for (int corner = 0; corner < 8; ++corner)
            includeRadius({corner & 1 ? maximum.x : minimum.x, corner & 2 ? maximum.y : minimum.y, corner & 4 ? maximum.z : minimum.z});
    }
    if (radius <= 0.0)
        throw std::domain_error("Level-set bounding radius must be positive.");
    for (Vec3& point : surfaceNodePositions_)
        point -= centroid;
    gridNodeOrigin_ -= centroid;
    centroidOffset_ = centroid;
    boundingRadius_ = radius;
    volume_ = volume;
    unitDensityInertiaTensor_ = inertia;
}

void LSInfo::reverseSDFSign() noexcept
{
    for (Real& value : gridNodeSignedDistance_)
    {
        value = -value;
    }
}

Real LSInfo::signedDistance(const Vec3& point) const
{
    if (!isValid())
    {
        throw std::logic_error("Configure the level-set geometry before evaluating its signed distance.");
    }
    if (!math::isFinite(point))
    {
        throw std::invalid_argument("Cannot evaluate a non-finite point.");
    }
    const Real value = evaluateSFD(point + centroidOffset_);
    if (!math::isFinite(value))
    {
        throw std::domain_error("Level-set evaluation produced a non-finite signed distance.");
    }
    return value;
}

Real LSInfo::radius() const noexcept
{
    Real result = 0.0;
    for (const Vec3& point : surfaceNodePositions_)
    {
        result = std::max(result, math::norm(point));
    }
    return result;
}

Real LSInfo::meanSurfaceNodeArea() const noexcept
{
    if (surfaceNodePositions_.empty())
    {
        return 0.0;
    }

    Real area = 0.0;
    for (const int3& face : surfaceTriangles_)
    {
        const Vec3& first = surfaceNodePositions_[face.x];
        const Vec3& second = surfaceNodePositions_[face.y];
        const Vec3& third = surfaceNodePositions_[face.z];
        area += 0.5 * math::norm(math::cross(second - first, third - first));
    }
    return area > 0.0 ? area / static_cast<int>(surfaceNodePositions_.size()) : gridNodeSpacing_ * gridNodeSpacing_;
}

void LSInfo::outputGridVTI(const std::string& fileName) const
{
    if (gridNodeSignedDistance_.empty())
    {
        throw std::logic_error("Build the level-set grid before writing it.");
    }

    std::filesystem::path path(fileName);
    if (path.extension() != ".vti")
    {
        path += ".vti";
    }
    if (!path.parent_path().empty())
    {
        std::error_code error;
        std::filesystem::create_directories(path.parent_path(), error);
        if (error)
        {
            throw std::runtime_error("Cannot create VTI output directory: " + error.message());
        }
    }

    std::ofstream output(path, std::ios::trunc);
    if (!output)
    {
        throw std::runtime_error("Cannot open VTI file for writing: " + path.string());
    }
    output << std::setprecision(17) << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           << "<ImageData Origin=\"" << gridNodeOrigin_.x << ' ' << gridNodeOrigin_.y << ' ' << gridNodeOrigin_.z << "\" Spacing=\"" << gridNodeSpacing_ << ' ' << gridNodeSpacing_ << ' '
           << gridNodeSpacing_ << "\" WholeExtent=\"0 " << gridNodeSize_.x - 1 << " 0 " << gridNodeSize_.y - 1 << " 0 " << gridNodeSize_.z - 1 << "\">\n"
           << "<Piece Extent=\"0 " << gridNodeSize_.x - 1 << " 0 " << gridNodeSize_.y - 1 << " 0 " << gridNodeSize_.z - 1 << "\">\n"
           << "<PointData Scalars=\"signedDistance\">\n<DataArray type=\"Float64\" Name=\"signedDistance\" format=\"ascii\">\n";
    for (Real value : gridNodeSignedDistance_)
    {
        output << value << '\n';
    }
    output << "</DataArray>\n</PointData>\n</Piece>\n</ImageData>\n</VTKFile>\n";
    if (!output)
    {
        throw std::runtime_error("Failed while writing VTI file: " + path.string());
    }
}

void LSInfo::buildImplicitSurfaceNode(int subdivisionLevel)
{
    if (!isValid())
    {
        throw std::logic_error("Configure the level-set geometry before building its surface.");
    }
    if (subdivisionLevel < 0)
    {
        throw std::invalid_argument("Subdivision level must be non-negative.");
    }
    clearGrid();
    const Vec3 boxMin = boundingBoxMin();
    const Vec3 boxMax = boundingBoxMax();
    if (!math::isFinite(boxMin) || !math::isFinite(boxMax) || boxMax.x <= boxMin.x || boxMax.y <= boxMin.y || boxMax.z <= boxMin.z)
    {
        throw std::domain_error("Cannot build a surface from an invalid bounding box.");
    }
    buildIcosahedron();
    for (Vec3& point : surfaceNodePositions_)
    {
        point = projectToSurface(point);
    }
    for (int level = 0; level < subdivisionLevel; ++level)
    {
        subdivideSurface();
    }
    orientSurfaceOutward();
}

void LSInfo::buildIcosahedron()
{
    const Real goldenRatio = 0.5 * (1.0 + std::sqrt(5.0));
    surfaceNodePositions_ = {{-1.0, goldenRatio, 0.0},
                             {1.0, goldenRatio, 0.0},
                             {-1.0, -goldenRatio, 0.0},
                             {1.0, -goldenRatio, 0.0},
                             {0.0, -1.0, goldenRatio},
                             {0.0, 1.0, goldenRatio},
                             {0.0, -1.0, -goldenRatio},
                             {0.0, 1.0, -goldenRatio},
                             {goldenRatio, 0.0, -1.0},
                             {goldenRatio, 0.0, 1.0},
                             {-goldenRatio, 0.0, -1.0},
                             {-goldenRatio, 0.0, 1.0}};
    surfaceTriangles_ = {triangle(0, 11, 5),  triangle(0, 5, 1),  triangle(0, 1, 7),  triangle(0, 7, 10), triangle(0, 10, 11), triangle(1, 5, 9), triangle(5, 11, 4),
                         triangle(11, 10, 2), triangle(10, 7, 6), triangle(7, 1, 8),  triangle(3, 9, 4),  triangle(3, 4, 2),   triangle(3, 2, 6), triangle(3, 6, 8),
                         triangle(3, 8, 9),   triangle(4, 9, 5),  triangle(2, 4, 11), triangle(6, 2, 10), triangle(8, 6, 7),   triangle(9, 8, 1)};
}

void LSInfo::subdivideSurface()
{
    clearGrid();
    std::unordered_map<std::uint64_t, int> midpointIndices;
    std::vector<int3> refinedTriangles;
    refinedTriangles.reserve(4 * surfaceTriangles_.size());

    const auto midpointIndex = [&](int first, int second)
    {
        const std::uint64_t key = edgeKey(first, second);
        const auto iterator = midpointIndices.find(key);
        if (iterator != midpointIndices.end())
        {
            return iterator->second;
        }
        const int index = static_cast<int>(surfaceNodePositions_.size());
        surfaceNodePositions_.push_back(projectToSurface(0.5 * (surfaceNodePositions_[first] + surfaceNodePositions_[second])));
        midpointIndices.emplace(key, index);
        return index;
    };

    for (const int3& face : surfaceTriangles_)
    {
        const int firstSecond = midpointIndex(face.x, face.y);
        const int secondThird = midpointIndex(face.y, face.z);
        const int thirdFirst = midpointIndex(face.z, face.x);
        refinedTriangles.push_back(triangle(face.x, firstSecond, thirdFirst));
        refinedTriangles.push_back(triangle(face.y, secondThird, firstSecond));
        refinedTriangles.push_back(triangle(face.z, thirdFirst, secondThird));
        refinedTriangles.push_back(triangle(firstSecond, secondThird, thirdFirst));
    }
    surfaceTriangles_ = std::move(refinedTriangles);
}

void LSInfo::orientSurfaceOutward()
{
    for (int3& face : surfaceTriangles_)
    {
        const Vec3& first = surfaceNodePositions_[face.x];
        const Vec3& second = surfaceNodePositions_[face.y];
        const Vec3& third = surfaceNodePositions_[face.z];
        if (math::dot(math::cross(second - first, third - first), first + second + third) < 0.0)
        {
            std::swap(face.y, face.z);
        }
    }
}

struct TriangleMesh::accelerationData {
    struct triangleReference {
        int triangleIndex_{0};
        Vec3 centroid_{Vec3::zero()};
    };

    struct node {
        Vec3 minimum_{Vec3::zero()};
        Vec3 maximum_{Vec3::zero()};
        Vec3 centroid_{Vec3::zero()};
        Vec3 areaNormal_{Vec3::zero()};
        int left_{-1};
        int right_{-1};
        int begin_{0};
        int count_{0};
        bool leaf() const noexcept { return count_ > 0; }
    };

    std::vector<Vec3> vertices_; ///< Native-frame vertex copy, independent of centroid-corrected output surface positions.
    const int3* triangles_{nullptr};
    std::vector<triangleReference> references_;
    std::vector<node> nodes_;

    void build(const std::vector<Vec3>& vertices, const std::vector<int3>& triangles)
    {
        vertices_ = vertices;
        triangles_ = triangles.data();
        references_.resize(triangles.size());
        for (int index = 0; index < static_cast<int>(triangles.size()); ++index)
        {
            const int3 face = triangles[index];
            references_[index] = {index, (vertices[face.x] + vertices[face.y] + vertices[face.z]) / 3.0};
        }
        nodes_.clear();
        nodes_.reserve(2 * triangles.size());
        buildNode(0, static_cast<int>(references_.size()));
    }

    int buildNode(int begin, int end)
    {
        node value;
        computeBounds(value.minimum_, value.maximum_, begin, end);
        const int count = end - begin;
        if (count <= 8)
        {
            value.begin_ = begin;
            value.count_ = count;
            nodes_.push_back(value);
            return static_cast<int>(nodes_.size()) - 1;
        }

        const Vec3 extent = value.maximum_ - value.minimum_;
        const int axis = extent.x > extent.y && extent.x > extent.z ? 0 : (extent.y > extent.z ? 1 : 2);
        const int middle = (begin + end) / 2;
        std::nth_element(references_.begin() + begin,
                         references_.begin() + middle,
                         references_.begin() + end,
                         [axis](const triangleReference& first, const triangleReference& second) { return first.centroid_[axis] < second.centroid_[axis]; });

        const int index = static_cast<int>(nodes_.size());
        nodes_.push_back(value);
        const int left = buildNode(begin, middle);
        const int right = buildNode(middle, end);
        nodes_[index].left_ = left;
        nodes_[index].right_ = right;
        computeCluster(nodes_[index].centroid_, nodes_[index].areaNormal_, begin, end);
        return index;
    }

    void computeBounds(Vec3& minimum, Vec3& maximum, int begin, int end) const noexcept
    {
        const int3 firstFace = triangles_[references_[begin].triangleIndex_];
        minimum = vertices_[firstFace.x];
        maximum = minimum;
        for (int index = begin; index < end; ++index)
        {
            const int3 face = triangles_[references_[index].triangleIndex_];
            minimum = math::componentMin(minimum, math::componentMin(vertices_[face.x], math::componentMin(vertices_[face.y], vertices_[face.z])));
            maximum = math::componentMax(maximum, math::componentMax(vertices_[face.x], math::componentMax(vertices_[face.y], vertices_[face.z])));
        }
    }

    void computeCluster(Vec3& centroid, Vec3& areaNormal, int begin, int end) const noexcept
    {
        centroid = Vec3::zero();
        areaNormal = Vec3::zero();
        Real areaSum = 0.0;
        for (int index = begin; index < end; ++index)
        {
            const int3 face = triangles_[references_[index].triangleIndex_];
            const Vec3& first = vertices_[face.x];
            const Vec3& second = vertices_[face.y];
            const Vec3& third = vertices_[face.z];
            const Vec3 twiceAreaNormal = math::cross(second - first, third - first);
            const Real area = 0.5 * math::norm(twiceAreaNormal);
            centroid += area * (first + second + third) / 3.0;
            areaNormal += 0.5 * twiceAreaNormal; // Sum A * n for the far-field solid-angle approximation.
            areaSum += area;
        }
        if (areaSum > 0.0)
        {
            centroid /= areaSum;
        }
    }

    static Real distanceToBox(const Vec3& point, const Vec3& minimum, const Vec3& maximum) noexcept
    {
        const Real dx = std::max({minimum.x - point.x, 0.0, point.x - maximum.x});
        const Real dy = std::max({minimum.y - point.y, 0.0, point.y - maximum.y});
        const Real dz = std::max({minimum.z - point.z, 0.0, point.z - maximum.z});
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    static Real pointTriangleDistance(const Vec3& point, const Vec3& first, const Vec3& second, const Vec3& third) noexcept
    {
        const Vec3 firstSecond = second - first;
        const Vec3 firstThird = third - first;
        const Vec3 firstPoint = point - first;
        const Real d1 = math::dot(firstSecond, firstPoint);
        const Real d2 = math::dot(firstThird, firstPoint);
        if (d1 <= 0.0 && d2 <= 0.0)
        {
            return math::norm(firstPoint);
        }

        const Vec3 secondPoint = point - second;
        const Real d3 = math::dot(firstSecond, secondPoint);
        const Real d4 = math::dot(firstThird, secondPoint);
        if (d3 >= 0.0 && d4 <= d3)
        {
            return math::norm(secondPoint);
        }

        const Real vc = d1 * d4 - d3 * d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            return math::norm(firstPoint - firstSecond * (d1 / (d1 - d3)));
        }

        const Vec3 thirdPoint = point - third;
        const Real d5 = math::dot(firstSecond, thirdPoint);
        const Real d6 = math::dot(firstThird, thirdPoint);
        if (d6 >= 0.0 && d5 <= d6)
        {
            return math::norm(thirdPoint);
        }

        const Real vb = d5 * d2 - d1 * d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            return math::norm(firstPoint - firstThird * (d2 / (d2 - d6)));
        }

        const Real va = d3 * d6 - d5 * d4;
        if (va <= 0.0 && d4 >= d3 && d5 >= d6)
        {
            const Vec3 secondThird = third - second;
            const Real weight = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return math::norm(secondPoint - weight * secondThird);
        }

        const Vec3 normal = math::cross(firstSecond, firstThird);
        return std::abs(math::dot(firstPoint, normal)) / math::norm(normal);
    }

    Real queryDistance(const Vec3& point, int nodeIndex, Real best) const noexcept
    {
        const node& value = nodes_[nodeIndex];
        if (distanceToBox(point, value.minimum_, value.maximum_) > best)
        {
            return best;
        }
        if (value.leaf())
        {
            for (int offset = 0; offset < value.count_; ++offset)
            {
                const int3 face = triangles_[references_[value.begin_ + offset].triangleIndex_];
                best = std::min(best, pointTriangleDistance(point, vertices_[face.x], vertices_[face.y], vertices_[face.z]));
            }
            return best;
        }
        best = queryDistance(point, value.left_, best);
        return queryDistance(point, value.right_, best);
    }

    static Real triangleSolidAngle(const Vec3& point, const Vec3& first, const Vec3& second, const Vec3& third) noexcept
    {
        const Vec3 firstVector = first - point;
        const Vec3 secondVector = second - point;
        const Vec3 thirdVector = third - point;
        const Real firstLength = math::norm(firstVector);
        const Real secondLength = math::norm(secondVector);
        const Real thirdLength = math::norm(thirdVector);
        const Real numerator = math::dot(firstVector, math::cross(secondVector, thirdVector));
        const Real denominator = firstLength * secondLength * thirdLength + math::dot(firstVector, secondVector) * thirdLength + math::dot(secondVector, thirdVector) * firstLength +
                                 math::dot(thirdVector, firstVector) * secondLength;
        return 2.0 * std::atan2(numerator, denominator);
    }

    Real solidAngle(const Vec3& point, int nodeIndex) const noexcept
    {
        const node& value = nodes_[nodeIndex];
        if (!value.leaf())
        {
            const Vec3 displacement = value.centroid_ - point;
            const Real distance = math::norm(displacement);
            const Real extent = math::norm(value.maximum_ - value.minimum_);
            if (distance > 4.0 * extent)
            {
                return math::dot(value.areaNormal_, displacement) / (distance * distance * distance);
            }
            return solidAngle(point, value.left_) + solidAngle(point, value.right_);
        }

        Real result = 0.0;
        for (int offset = 0; offset < value.count_; ++offset)
        {
            const int3 face = triangles_[references_[value.begin_ + offset].triangleIndex_];
            result += triangleSolidAngle(point, vertices_[face.x], vertices_[face.y], vertices_[face.z]);
        }
        return result;
    }
};

TriangleMesh::TriangleMesh() : acceleration_(std::make_unique<accelerationData>()) {}

TriangleMesh::TriangleMesh(const std::vector<Vec3>& vertexPositions, const std::vector<int3>& triangles) : TriangleMesh()
{
    surfaceNodePositions_ = vertexPositions;
    surfaceTriangles_ = triangles;
    initializeMesh();
}

TriangleMesh::~TriangleMesh() = default;
TriangleMesh::TriangleMesh(TriangleMesh&&) noexcept = default;
TriangleMesh& TriangleMesh::operator=(TriangleMesh&&) noexcept = default;

bool TriangleMesh::isValid() const noexcept { return configured_; }

Vec3 TriangleMesh::boundingBoxMin() const noexcept { return acceleration_->nodes_.front().minimum_; }

Vec3 TriangleMesh::boundingBoxMax() const noexcept { return acceleration_->nodes_.front().maximum_; }

Real TriangleMesh::evaluateSFD(const Vec3& point) const noexcept
{
    const Real distance = acceleration_->queryDistance(point, 0, std::numeric_limits<Real>::max());
    const Real windingNumber = acceleration_->solidAngle(point, 0) / (4.0 * math::pi);
    return windingNumber > 0.5 ? -distance : distance;
}

void TriangleMesh::initializeMesh()
{
    configured_ = false;
    clearGrid();
    if (surfaceNodePositions_.empty() || surfaceTriangles_.empty())
    {
        throw std::invalid_argument("Triangle mesh requires vertices and triangles.");
    }
    for (const Vec3& point : surfaceNodePositions_)
    {
        if (!math::isFinite(point))
        {
            throw std::invalid_argument("Triangle mesh contains a non-finite vertex.");
        }
    }

    struct adjacency {
        int neighbor_{0};
        int first_{0};
        int second_{0};
    };
    std::unordered_map<std::uint64_t, std::vector<int>> edgeTriangles;
    for (int triangleIndex = 0; triangleIndex < static_cast<int>(surfaceTriangles_.size()); ++triangleIndex)
    {
        const int3 face = surfaceTriangles_[triangleIndex];
        const int vertexCount = static_cast<int>(surfaceNodePositions_.size());
        if (face.x < 0 || face.x >= vertexCount || face.y < 0 || face.y >= vertexCount || face.z < 0 || face.z >= vertexCount || face.x == face.y || face.y == face.z || face.z == face.x)
        {
            throw std::invalid_argument("Triangle mesh contains invalid connectivity.");
        }
        if (math::norm(math::cross(surfaceNodePositions_[face.y] - surfaceNodePositions_[face.x], surfaceNodePositions_[face.z] - surfaceNodePositions_[face.x])) <= math::defaultTolerance)
        {
            throw std::invalid_argument("Triangle mesh contains a degenerate triangle.");
        }
        edgeTriangles[edgeKey(face.x, face.y)].push_back(triangleIndex);
        edgeTriangles[edgeKey(face.y, face.z)].push_back(triangleIndex);
        edgeTriangles[edgeKey(face.z, face.x)].push_back(triangleIndex);
    }

    std::vector<std::vector<adjacency>> neighbors(surfaceTriangles_.size());
    for (const auto& [key, incidentTriangles] : edgeTriangles)
    {
        if (incidentTriangles.size() != 2)
        {
            throw std::invalid_argument("Triangle mesh must be a closed two-manifold surface.");
        }
        const int first = static_cast<int>(key >> 32U);
        const int second = static_cast<int>(key & 0xFFFFFFFFU);
        neighbors[incidentTriangles[0]].push_back({incidentTriangles[1], first, second});
        neighbors[incidentTriangles[1]].push_back({incidentTriangles[0], first, second});
    }

    const auto edgeDirection = [](const int3& face, int first, int second) noexcept
    {
        if ((face.x == first && face.y == second) || (face.y == first && face.z == second) || (face.z == first && face.x == second))
        {
            return 1;
        }
        return -1;
    };

    std::vector<unsigned char> visited(surfaceTriangles_.size(), 0);
    std::vector<std::vector<int>> components;
    for (int start = 0; start < static_cast<int>(surfaceTriangles_.size()); ++start)
    {
        if (visited[start] != 0)
        {
            continue;
        }
        components.emplace_back();
        std::queue<int> pending;
        pending.push(start);
        visited[start] = 1;
        while (!pending.empty())
        {
            const int current = pending.front();
            pending.pop();
            components.back().push_back(current);
            for (const adjacency& adjacent : neighbors[current])
            {
                const bool sameDirection =
                    edgeDirection(surfaceTriangles_[current], adjacent.first_, adjacent.second_) == edgeDirection(surfaceTriangles_[adjacent.neighbor_], adjacent.first_, adjacent.second_);
                if (visited[adjacent.neighbor_] == 0)
                {
                    if (sameDirection)
                    {
                        std::swap(surfaceTriangles_[adjacent.neighbor_].y, surfaceTriangles_[adjacent.neighbor_].z);
                    }
                    visited[adjacent.neighbor_] = 1;
                    pending.push(adjacent.neighbor_);
                }
                else if (sameDirection)
                {
                    throw std::invalid_argument("Triangle mesh orientation is inconsistent.");
                }
            }
        }
    }

    for (const std::vector<int>& component : components)
    {
        Real signedVolumeTimesSix = 0.0;
        for (int triangleIndex : component)
        {
            const int3 face = surfaceTriangles_[triangleIndex];
            signedVolumeTimesSix += math::dot(surfaceNodePositions_[face.x], math::cross(surfaceNodePositions_[face.y], surfaceNodePositions_[face.z]));
        }
        if (std::abs(signedVolumeTimesSix) <= math::defaultTolerance)
        {
            throw std::invalid_argument("Triangle mesh encloses no measurable volume.");
        }
        if (signedVolumeTimesSix < 0.0)
        {
            for (int triangleIndex : component)
            {
                std::swap(surfaceTriangles_[triangleIndex].y, surfaceTriangles_[triangleIndex].z);
            }
        }
    }
    acceleration_->build(surfaceNodePositions_, surfaceTriangles_);
    configured_ = true;
}

void TriangleMesh::loadOBJ(const std::string& fileName)
{
    std::ifstream input(fileName);
    if (!input)
    {
        throw std::runtime_error("Cannot open OBJ file: " + fileName);
    }

    std::vector<Vec3> vertices;
    std::vector<int3> triangles;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line))
    {
        ++lineNumber;
        std::istringstream stream(line);
        std::string tag;
        stream >> tag;
        if (tag.empty() || tag[0] == '#')
        {
            continue;
        }
        if (tag == "v")
        {
            Vec3 point;
            if (!(stream >> point.x >> point.y >> point.z) || !math::isFinite(point))
            {
                throw std::runtime_error("Invalid OBJ vertex at line " + std::to_string(lineNumber) + '.');
            }
            vertices.push_back(point);
        }
        else if (tag == "f")
        {
            std::vector<int> face;
            std::string token;
            while (stream >> token)
            {
                const std::string indexText = token.substr(0, token.find('/'));
                if (indexText.empty())
                {
                    throw std::runtime_error("Invalid OBJ face at line " + std::to_string(lineNumber) + '.');
                }
                int index = 0;
                const char* indexEnd = indexText.data() + indexText.size();
                const auto [parsedEnd, error] = std::from_chars(indexText.data(), indexEnd, index);
                if (error != std::errc{} || parsedEnd != indexEnd)
                {
                    throw std::runtime_error("Invalid OBJ index at line " + std::to_string(lineNumber) + '.');
                }
                index = index > 0 ? index - 1 : static_cast<int>(vertices.size()) + index;
                if (index < 0 || index >= static_cast<int>(vertices.size()))
                {
                    throw std::runtime_error("OBJ face index is out of range at line " + std::to_string(lineNumber) + '.');
                }
                face.push_back(index);
            }
            if (face.size() < 3)
            {
                throw std::runtime_error("OBJ face has fewer than three vertices at line " + std::to_string(lineNumber) + '.');
            }
            for (int index = 1; index + 1 < static_cast<int>(face.size()); ++index)
            {
                triangles.push_back(triangle(face[0], face[index], face[index + 1]));
            }
        }
    }

    TriangleMesh loaded(vertices, triangles);
    *this = std::move(loaded);
}

void TriangleMesh::fineMesh()
{
    if (surfaceNodePositions_.empty() || surfaceTriangles_.empty() || !acceleration_ || acceleration_->nodes_.empty())
    {
        throw std::logic_error("Load a valid triangle mesh before refining it.");
    }
    configured_ = false;
    subdivideSurface();
    initializeMesh();
}

std::unique_ptr<TriangleMesh> makeRandomShape(Real radius, Real minimumSurfaceHeight, Real maximumSurfaceHeight, int subdivisionLevel, std::uint64_t seed)
{
    const Real minimumRadius = radius + minimumSurfaceHeight;
    const Real maximumRadius = radius + maximumSurfaceHeight;
    if (!finitePositive(radius) || !math::isFinite(minimumSurfaceHeight) || !math::isFinite(maximumSurfaceHeight) || minimumSurfaceHeight > maximumSurfaceHeight ||
        !finitePositive(minimumRadius) || !finitePositive(maximumRadius))
    {
        throw std::invalid_argument("Random shape requires a positive finite radius, ordered finite height bounds, and positive finite deformed radii.");
    }
    if (subdivisionLevel < 0)
    {
        throw std::invalid_argument("Subdivision level must be non-negative.");
    }
    // The mesh acceleration structure uses int triangle indices. Check growth
    // before sampling, without imposing an unrelated rendering/detail limit.
    int faceCount = 20;
    for (int level = 0; level < subdivisionLevel; ++level)
    {
        if (faceCount > std::numeric_limits<int>::max() / 4)
        {
            throw std::invalid_argument("Random shape subdivision exceeds the mesh index capacity.");
        }
        faceCount *= 4;
    }

    Sphere sphere(1.0);
    sphere.buildSurfaceNode(subdivisionLevel);
    std::vector<Vec3> vertices = sphere.surfaceNodePosition();

    // A few low-frequency waves form a continuous correlated field on the unit
    // sphere. Unlike independent vertex noise, refinement preserves smooth lobes.
    struct surfaceWave {
        Vec3 waveVector;
        Real phase;
        Real amplitude;
    };
    std::mt19937_64 random(seed);
    const auto uniform = [&random]() { return static_cast<Real>(random() >> 11U) * 0x1.0p-53; };
    std::array<surfaceWave, 12> waves;
    for (std::size_t index = 0; index < waves.size(); ++index)
    {
        const Real z = 2.0 * uniform() - 1.0;
        const Real azimuth = math::twoPi * uniform();
        const Real equatorialRadius = std::sqrt(std::max(0.0, 1.0 - z * z));
        const Real band = 1.0 + static_cast<Real>(index / 4);
        waves[index] = {(1.0 + band) * Vec3{equatorialRadius * std::cos(azimuth), equatorialRadius * std::sin(azimuth), z},
                        math::twoPi * uniform(),
                        (0.5 + 0.5 * uniform()) / (band * band)};
    }
    std::vector<Real> heights(vertices.size(), 0.0);
    for (std::size_t index = 0; index < vertices.size(); ++index)
    {
        for (const surfaceWave& wave : waves)
        {
            heights[index] += wave.amplitude * std::sin(math::dot(wave.waveVector, vertices[index]) + wave.phase);
        }
    }
    const auto [minimum, maximum] = std::minmax_element(heights.begin(), heights.end());
    const Real heightRange = *maximum - *minimum;
    for (std::size_t index = 0; index < vertices.size(); ++index)
    {
        const Real fraction = heightRange > 0.0 ? std::clamp((heights[index] - *minimum) / heightRange, 0.0, 1.0) : 0.5;
        vertices[index] *= minimumRadius + fraction * (maximumRadius - minimumRadius);
    }

    // Positive radial scaling preserves the icosphere's topology and winding.
    // Reject arithmetic overflow before it can enter mesh acceleration queries;
    // TriangleMesh also checks closure and its existing area/volume tolerances.
    Real signedVolumeTimesSix = 0.0;
    for (const int3& face : sphere.surfaceNodeConnectivity())
    {
        const Vec3& first = vertices[face.x];
        const Vec3& second = vertices[face.y];
        const Vec3& third = vertices[face.z];
        const Vec3 normal = math::cross(second - first, third - first);
        const Real volumeTimesSix = math::dot(first, math::cross(second, third));
        signedVolumeTimesSix += volumeTimesSix;
        if (!math::isFinite(normal) || !finitePositive(volumeTimesSix) || !math::isFinite(signedVolumeTimesSix))
        {
            throw std::invalid_argument("Random shape dimensions exceed the numerical range of a non-degenerate triangle mesh.");
        }
    }
    return std::make_unique<TriangleMesh>(vertices, sphere.surfaceNodeConnectivity());
}

Sphere::Sphere(Real radius) { setParameter(radius); }

void Sphere::setParameter(Real radius)
{
    if (!finitePositive(radius))
    {
        throw std::invalid_argument("Sphere radius must be positive and finite.");
    }
    radius_ = radius;
    configured_ = true;
    surfaceNodePositions_.clear();
    surfaceTriangles_.clear();
    clearGrid();
}

void Sphere::buildSurfaceNode(int subdivisionLevel) { buildImplicitSurfaceNode(subdivisionLevel); }

bool Sphere::isValid() const noexcept { return configured_; }
Vec3 Sphere::boundingBoxMin() const noexcept { return Vec3{-radius_}; }
Vec3 Sphere::boundingBoxMax() const noexcept { return Vec3{radius_}; }
Real Sphere::evaluateSFD(const Vec3& point) const noexcept { return math::norm(point) - radius_; }

Vec3 Sphere::projectToSurface(const Vec3& point) const noexcept
{
    const Vec3 direction = math::normalizedOrZero(point);
    return direction == Vec3::zero() ? Vec3::unitX() * radius_ : direction * radius_;
}

Superellipsoid::Superellipsoid(Real radiusX, Real radiusY, Real radiusZ, Real equatorialExponent, Real polarExponent) { setParameter(radiusX, radiusY, radiusZ, equatorialExponent, polarExponent); }

void Superellipsoid::setParameter(Real radiusX, Real radiusY, Real radiusZ, Real equatorialExponent, Real polarExponent)
{
    if (!finitePositive(radiusX) || !finitePositive(radiusY) || !finitePositive(radiusZ) || !finitePositive(equatorialExponent) || !finitePositive(polarExponent))
    {
        throw std::invalid_argument("Superellipsoid radii and exponents must be positive and finite.");
    }
    radiusX_ = radiusX;
    radiusY_ = radiusY;
    radiusZ_ = radiusZ;
    equatorialExponent_ = equatorialExponent;
    polarExponent_ = polarExponent;
    configured_ = true;
    surfaceNodePositions_.clear();
    surfaceTriangles_.clear();
    clearGrid();
}

void Superellipsoid::buildSurfaceNode(int subdivisionLevel) { buildImplicitSurfaceNode(subdivisionLevel); }

Real Superellipsoid::radius() const noexcept { return std::sqrt(radiusX_ * radiusX_ + radiusY_ * radiusY_ + radiusZ_ * radiusZ_); }

bool Superellipsoid::isValid() const noexcept { return configured_; }

Vec3 Superellipsoid::boundingBoxMin() const noexcept { return {-radiusX_, -radiusY_, -radiusZ_}; }
Vec3 Superellipsoid::boundingBoxMax() const noexcept { return {radiusX_, radiusY_, radiusZ_}; }

Real Superellipsoid::evaluateSFD(const Vec3& point) const noexcept
{
    const Real pointRadius = math::norm(point);
    if (pointRadius <= math::defaultTolerance)
    {
        return -std::min({radiusX_, radiusY_, radiusZ_});
    }

    const Vec3 surfacePoint = projectToSurface(point);
    const Real surfaceRadius = math::norm(surfacePoint);
    const Vec3 gradient = evaluateImplicitGradient(surfacePoint);
    const Vec3 normal = math::normalizedOrZero(gradient);
    if (normal == Vec3::zero())
    {
        return pointRadius - surfaceRadius;
    }
    const Vec3 radialDirection = point / pointRadius;
    return (pointRadius - surfaceRadius) * std::min(1.0, std::abs(math::dot(radialDirection, normal)));
}

Vec3 Superellipsoid::projectToSurface(const Vec3& point) const noexcept
{
    const Vec3 direction = math::normalizedOrZero(point);
    if (direction == Vec3::zero())
    {
        return Vec3::zero();
    }

    const Real negativeInfinity = -std::numeric_limits<Real>::infinity();
    const auto logPower = [negativeInfinity](Real value, Real exponent) { return value > 0.0 ? exponent * std::log(value) : negativeInfinity; };
    const auto logAddExp = [negativeInfinity](Real first, Real second)
    {
        if (first == negativeInfinity)
        {
            return second;
        }
        if (second == negativeInfinity)
        {
            return first;
        }
        const Real maximum = std::max(first, second);
        return maximum + std::log(std::exp(first - maximum) + std::exp(second - maximum));
    };

    const Real exponentXY = 2.0 / equatorialExponent_;
    const Real exponentZ = 2.0 / polarExponent_;
    const Real logXY = logAddExp(logPower(std::abs(direction.x / radiusX_), exponentXY), logPower(std::abs(direction.y / radiusY_), exponentXY));
    const Real logXYBlock = logXY == negativeInfinity ? negativeInfinity : (equatorialExponent_ / polarExponent_) * logXY;
    const Real logMeasure = logAddExp(logXYBlock, logPower(std::abs(direction.z / radiusZ_), exponentZ));
    return direction * std::exp(-0.5 * polarExponent_ * logMeasure);
}

Vec3 Superellipsoid::evaluateImplicitGradient(const Vec3& point) const noexcept
{
    const Real absoluteX = std::abs(point.x / radiusX_);
    const Real absoluteY = std::abs(point.y / radiusY_);
    const Real absoluteZ = std::abs(point.z / radiusZ_);
    const Real exponentXY = 2.0 / equatorialExponent_;
    const Real exponentZ = 2.0 / polarExponent_;
    const Real xy = std::pow(absoluteX, exponentXY) + std::pow(absoluteY, exponentXY);

    Real gradientX = 0.0;
    Real gradientY = 0.0;
    if (xy > 0.0)
    {
        const Real common = (2.0 / polarExponent_) * std::pow(xy, equatorialExponent_ / polarExponent_ - 1.0);
        if (absoluteX > 0.0)
        {
            gradientX = std::copysign(common * std::pow(absoluteX, exponentXY - 1.0) / radiusX_, point.x);
        }
        if (absoluteY > 0.0)
        {
            gradientY = std::copysign(common * std::pow(absoluteY, exponentXY - 1.0) / radiusY_, point.y);
        }
    }
    const Real gradientZ = absoluteZ > 0.0 ? std::copysign((2.0 / polarExponent_) * std::pow(absoluteZ, exponentZ - 1.0) / radiusZ_, point.z) : 0.0;
    return {gradientX, gradientY, gradientZ};
}

PlaneWall::PlaneWall(const Vec3& outwardNormal, Real size) { setParameter(outwardNormal, size); }

void PlaneWall::setParameter(const Vec3& outwardNormal, Real size)
{
    const Vec3 normalized = math::normalizedOrZero(outwardNormal);
    if (normalized == Vec3::zero() || !finitePositive(size))
    {
        throw std::invalid_argument("Plane wall normal and size are invalid.");
    }
    outwardNormal_ = normalized;
    size_ = size;
    configured_ = true;
    clearGrid();
    buildSurfaceNode();
}

void PlaneWall::buildSurfaceNode()
{
    clearGrid();
    Vec3 tangent1;
    Vec3 tangent2;
    orthogonalBasis(outwardNormal_, tangent1, tangent2);
    const Real halfSize = 0.5 * size_;
    surfaceNodePositions_ = {-halfSize * tangent1 - halfSize * tangent2,
                             halfSize * tangent1 - halfSize * tangent2,
                             halfSize * tangent1 + halfSize * tangent2,
                             -halfSize * tangent1 + halfSize * tangent2};
    surfaceTriangles_ = {triangle(0, 1, 2), triangle(0, 2, 3)};
}

bool PlaneWall::isValid() const noexcept { return configured_; }
Vec3 PlaneWall::boundingBoxMin() const noexcept { return Vec3{-0.5 * size_}; }
Vec3 PlaneWall::boundingBoxMax() const noexcept { return Vec3{0.5 * size_}; }
Real PlaneWall::evaluateSFD(const Vec3& point) const noexcept { return math::dot(point, outwardNormal_); }

BoxWall::BoxWall(const Vec3& size) { setParameter(size); }

void BoxWall::setParameter(const Vec3& size)
{
    if (!math::isFinite(size) || size.x <= 0.0 || size.y <= 0.0 || size.z <= 0.0)
    {
        throw std::invalid_argument("Box dimensions must be positive and finite.");
    }
    size_ = size;
    configured_ = true;
    clearGrid();
    buildSurfaceNode();
}

void BoxWall::buildSurfaceNode()
{
    clearGrid();
    const Vec3 half = 0.5 * size_;
    surfaceNodePositions_ = {{-half.x, -half.y, -half.z},
                             {half.x, -half.y, -half.z},
                             {half.x, half.y, -half.z},
                             {-half.x, half.y, -half.z},
                             {-half.x, -half.y, half.z},
                             {half.x, -half.y, half.z},
                             {half.x, half.y, half.z},
                             {-half.x, half.y, half.z}};
    surfaceTriangles_ = {triangle(0, 2, 1),
                         triangle(0, 3, 2),
                         triangle(4, 5, 6),
                         triangle(4, 6, 7),
                         triangle(0, 1, 5),
                         triangle(0, 5, 4),
                         triangle(3, 7, 6),
                         triangle(3, 6, 2),
                         triangle(0, 4, 7),
                         triangle(0, 7, 3),
                         triangle(1, 2, 6),
                         triangle(1, 6, 5)};
}

bool BoxWall::isValid() const noexcept { return configured_; }
Vec3 BoxWall::boundingBoxMin() const noexcept { return -0.5 * size_; }
Vec3 BoxWall::boundingBoxMax() const noexcept { return 0.5 * size_; }

Real BoxWall::evaluateSFD(const Vec3& point) const noexcept
{
    const Vec3 half = 0.5 * size_;
    const Vec3 distance{std::abs(point.x) - half.x, std::abs(point.y) - half.y, std::abs(point.z) - half.z};
    const Vec3 outside{std::max(distance.x, 0.0), std::max(distance.y, 0.0), std::max(distance.z, 0.0)};
    return math::norm(outside) + std::min(std::max({distance.x, distance.y, distance.z}), 0.0);
}

void BoxParticle::buildSurfaceNode(Real surfaceSpacing)
{
    const Vec3 size = boxSize3D();
    if (!finitePositive(surfaceSpacing))
    {
        throw std::invalid_argument("Box surface spacing must be positive and finite.");
    }
    clearGrid();

    const int countX = std::max(1, static_cast<int>(std::ceil(size.x / surfaceSpacing)));
    const int countY = std::max(1, static_cast<int>(std::ceil(size.y / surfaceSpacing)));
    const int countZ = std::max(1, static_cast<int>(std::ceil(size.z / surfaceSpacing)));
    const Vec3 spacing{size.x / countX, size.y / countY, size.z / countZ};
    const Vec3 half = 0.5 * size;
    const auto storageIndex = [countX, countY](int x, int y, int z) { return (z * (countY + 1) + y) * (countX + 1) + x; };
    std::vector<int> nodeIndices((countX + 1) * (countY + 1) * (countZ + 1), -1);
    surfaceNodePositions_.clear();
    surfaceTriangles_.clear();

    for (int z = 0; z <= countZ; ++z)
    {
        for (int y = 0; y <= countY; ++y)
        {
            for (int x = 0; x <= countX; ++x)
            {
                if (x != 0 && x != countX && y != 0 && y != countY && z != 0 && z != countZ)
                {
                    continue;
                }
                nodeIndices[storageIndex(x, y, z)] = static_cast<int>(surfaceNodePositions_.size());
                surfaceNodePositions_.push_back(-half + Vec3{x * spacing.x, y * spacing.y, z * spacing.z});
            }
        }
    }

    const auto addTriangle = [&](int x0, int y0, int z0, int x1, int y1, int z1, int x2, int y2, int z2)
    { surfaceTriangles_.push_back(triangle(nodeIndices[storageIndex(x0, y0, z0)], nodeIndices[storageIndex(x1, y1, z1)], nodeIndices[storageIndex(x2, y2, z2)])); };
    for (int x = 0; x < countX; ++x)
    {
        for (int y = 0; y < countY; ++y)
        {
            addTriangle(x, y, 0, x + 1, y + 1, 0, x + 1, y, 0);
            addTriangle(x, y, 0, x, y + 1, 0, x + 1, y + 1, 0);
            addTriangle(x, y, countZ, x + 1, y, countZ, x + 1, y + 1, countZ);
            addTriangle(x, y, countZ, x + 1, y + 1, countZ, x, y + 1, countZ);
        }
    }
    for (int x = 0; x < countX; ++x)
    {
        for (int z = 0; z < countZ; ++z)
        {
            addTriangle(x, 0, z, x + 1, 0, z, x + 1, 0, z + 1);
            addTriangle(x, 0, z, x + 1, 0, z + 1, x, 0, z + 1);
            addTriangle(x, countY, z, x + 1, countY, z + 1, x + 1, countY, z);
            addTriangle(x, countY, z, x, countY, z + 1, x + 1, countY, z + 1);
        }
    }
    for (int y = 0; y < countY; ++y)
    {
        for (int z = 0; z < countZ; ++z)
        {
            addTriangle(0, y, z, 0, y + 1, z + 1, 0, y + 1, z);
            addTriangle(0, y, z, 0, y, z + 1, 0, y + 1, z + 1);
            addTriangle(countX, y, z, countX, y + 1, z, countX, y + 1, z + 1);
            addTriangle(countX, y, z, countX, y + 1, z + 1, countX, y, z + 1);
        }
    }
}

Real BoxParticle::meanSurfaceNodeArea() const noexcept
{
    const Vec3 size = boxSize3D();
    return surfaceNodePositions_.empty() ? 0.0 : 2.0 * (size.x * size.y + size.x * size.z + size.y * size.z) / static_cast<int>(surfaceNodePositions_.size());
}

CylinderWall::CylinderWall(const Vec3& bottomCenter, const Vec3& topCenter, Real radius) { setParameter(bottomCenter, topCenter, radius); }

void CylinderWall::setParameter(const Vec3& bottomCenter, const Vec3& topCenter, Real radius)
{
    if (!math::isFinite(bottomCenter) || !math::isFinite(topCenter) || !finitePositive(radius) || math::norm(topCenter - bottomCenter) <= math::defaultTolerance)
    {
        throw std::invalid_argument("Cylinder endpoints and radius are invalid.");
    }
    bottomCenter_ = bottomCenter;
    topCenter_ = topCenter;
    radius_ = radius;
    configured_ = true;
    clearGrid();
    buildSurfaceNode(segmentCount_);
}

void CylinderWall::setCircumferentialSegments(int segmentCount)
{
    if (segmentCount < 3)
    {
        throw std::invalid_argument("Cylinder requires at least three circumferential segments.");
    }
    segmentCount_ = segmentCount;
    buildSurfaceNode(segmentCount_);
}

void CylinderWall::buildSurfaceNode(int segmentCount)
{
    clearGrid();
    const Vec3 axis = math::normalizedOrZero(topCenter_ - bottomCenter_);
    Vec3 tangent1;
    Vec3 tangent2;
    orthogonalBasis(axis, tangent1, tangent2);
    surfaceNodePositions_.clear();
    surfaceTriangles_.clear();
    surfaceNodePositions_.reserve(2 * segmentCount + 2);
    surfaceTriangles_.reserve(4 * segmentCount);

    for (int index = 0; index < segmentCount; ++index)
    {
        const Real angle = math::twoPi * index / segmentCount;
        const Vec3 direction = std::cos(angle) * tangent1 + std::sin(angle) * tangent2;
        surfaceNodePositions_.push_back(bottomCenter_ + radius_ * direction);
        surfaceNodePositions_.push_back(topCenter_ + radius_ * direction);
    }
    const int bottomCenterIndex = static_cast<int>(surfaceNodePositions_.size());
    surfaceNodePositions_.push_back(bottomCenter_);
    const int topCenterIndex = static_cast<int>(surfaceNodePositions_.size());
    surfaceNodePositions_.push_back(topCenter_);

    for (int index = 0; index < segmentCount; ++index)
    {
        const int next = (index + 1) % segmentCount;
        const int bottom = 2 * index;
        const int top = bottom + 1;
        const int nextBottom = 2 * next;
        const int nextTop = nextBottom + 1;
        surfaceTriangles_.push_back(triangle(bottom, nextBottom, nextTop));
        surfaceTriangles_.push_back(triangle(bottom, nextTop, top));
        surfaceTriangles_.push_back(triangle(bottomCenterIndex, nextBottom, bottom));
        surfaceTriangles_.push_back(triangle(topCenterIndex, top, nextTop));
    }
}

bool CylinderWall::isValid() const noexcept { return configured_; }

Vec3 CylinderWall::boundingBoxMin() const noexcept
{
    const Vec3 axis = math::normalizedOrZero(topCenter_ - bottomCenter_);
    const Vec3 radialExtent{radius_ * std::sqrt(std::max(0.0, 1.0 - axis.x * axis.x)),
                            radius_ * std::sqrt(std::max(0.0, 1.0 - axis.y * axis.y)),
                            radius_ * std::sqrt(std::max(0.0, 1.0 - axis.z * axis.z))};
    return math::componentMin(bottomCenter_, topCenter_) - radialExtent;
}

Vec3 CylinderWall::boundingBoxMax() const noexcept
{
    const Vec3 axis = math::normalizedOrZero(topCenter_ - bottomCenter_);
    const Vec3 radialExtent{radius_ * std::sqrt(std::max(0.0, 1.0 - axis.x * axis.x)),
                            radius_ * std::sqrt(std::max(0.0, 1.0 - axis.y * axis.y)),
                            radius_ * std::sqrt(std::max(0.0, 1.0 - axis.z * axis.z))};
    return math::componentMax(bottomCenter_, topCenter_) + radialExtent;
}

Real CylinderWall::evaluateSFD(const Vec3& point) const noexcept
{
    const Vec3 axis = topCenter_ - bottomCenter_;
    const Vec3 relative = point - bottomCenter_;
    const Real axisSquared = math::dot(axis, axis);
    const Real projection = math::dot(relative, axis);
    const Vec3 radialTerm = relative * axisSquared - axis * projection;
    const Real radialDistance = math::norm(radialTerm) - radius_ * axisSquared;
    const Real axialDistance = std::abs(projection - 0.5 * axisSquared) - 0.5 * axisSquared;
    const Real radialSquared = radialDistance * radialDistance;
    const Real axialSquared = axialDistance * axialDistance * axisSquared;
    const Real signedSquaredDistance = std::max(radialDistance, axialDistance) < 0.0
                                           ? -std::min(radialSquared, axialSquared)
                                           : std::max(radialDistance, 0.0) * std::max(radialDistance, 0.0) + std::max(axialDistance, 0.0) * std::max(axialDistance, 0.0) * axisSquared;
    return std::copysign(std::sqrt(std::abs(signedSquaredDistance)) / axisSquared, signedSquaredDistance);
}

ConeWall::ConeWall(const Vec3& bottomCenter, const Vec3& topCenter, Real bottomRadius, Real topRadius) { setParameter(bottomCenter, topCenter, bottomRadius, topRadius); }

void ConeWall::setParameter(const Vec3& bottomCenter, const Vec3& topCenter, Real bottomRadius, Real topRadius)
{
    if (!math::isFinite(bottomCenter) || !math::isFinite(topCenter) || !math::isFinite(bottomRadius) || !math::isFinite(topRadius) || bottomRadius < 0.0 || topRadius < 0.0 ||
        (bottomRadius <= math::defaultTolerance && topRadius <= math::defaultTolerance) || math::norm(topCenter - bottomCenter) <= math::defaultTolerance)
    {
        throw std::invalid_argument("Cone endpoints and radii are invalid.");
    }
    bottomCenter_ = bottomCenter;
    topCenter_ = topCenter;
    bottomRadius_ = bottomRadius;
    topRadius_ = topRadius;
    configured_ = true;
    clearGrid();
    buildSurfaceNode(segmentCount_);
}

void ConeWall::setCircumferentialSegments(int segmentCount)
{
    if (segmentCount < 3)
    {
        throw std::invalid_argument("Cone requires at least three circumferential segments.");
    }
    segmentCount_ = segmentCount;
    buildSurfaceNode(segmentCount_);
}

void ConeWall::buildSurfaceNode(int segmentCount)
{
    clearGrid();
    const Vec3 axis = math::normalizedOrZero(topCenter_ - bottomCenter_);
    Vec3 tangent1;
    Vec3 tangent2;
    orthogonalBasis(axis, tangent1, tangent2);
    surfaceNodePositions_.clear();
    surfaceTriangles_.clear();

    const bool hasBottomRing = bottomRadius_ > math::defaultTolerance;
    const bool hasTopRing = topRadius_ > math::defaultTolerance;
    const int bottomBegin = hasBottomRing ? 0 : -1;
    if (hasBottomRing)
    {
        for (int index = 0; index < segmentCount; ++index)
        {
            const Real angle = math::twoPi * index / segmentCount;
            surfaceNodePositions_.push_back(bottomCenter_ + bottomRadius_ * (std::cos(angle) * tangent1 + std::sin(angle) * tangent2));
        }
    }
    const int topBegin = hasTopRing ? static_cast<int>(surfaceNodePositions_.size()) : -1;
    if (hasTopRing)
    {
        for (int index = 0; index < segmentCount; ++index)
        {
            const Real angle = math::twoPi * index / segmentCount;
            surfaceNodePositions_.push_back(topCenter_ + topRadius_ * (std::cos(angle) * tangent1 + std::sin(angle) * tangent2));
        }
    }

    int bottomApex = -1;
    int topApex = -1;
    if (!hasBottomRing)
    {
        bottomApex = static_cast<int>(surfaceNodePositions_.size());
        surfaceNodePositions_.push_back(bottomCenter_);
    }
    if (!hasTopRing)
    {
        topApex = static_cast<int>(surfaceNodePositions_.size());
        surfaceNodePositions_.push_back(topCenter_);
    }

    for (int index = 0; index < segmentCount; ++index)
    {
        const int next = (index + 1) % segmentCount;
        if (hasBottomRing && hasTopRing)
        {
            surfaceTriangles_.push_back(triangle(bottomBegin + index, bottomBegin + next, topBegin + next));
            surfaceTriangles_.push_back(triangle(bottomBegin + index, topBegin + next, topBegin + index));
        }
        else if (hasBottomRing)
        {
            surfaceTriangles_.push_back(triangle(bottomBegin + index, bottomBegin + next, topApex));
        }
        else
        {
            surfaceTriangles_.push_back(triangle(bottomApex, topBegin + next, topBegin + index));
        }
    }

    if (hasBottomRing)
    {
        const int center = static_cast<int>(surfaceNodePositions_.size());
        surfaceNodePositions_.push_back(bottomCenter_);
        for (int index = 0; index < segmentCount; ++index)
        {
            surfaceTriangles_.push_back(triangle(center, bottomBegin + (index + 1) % segmentCount, bottomBegin + index));
        }
    }
    if (hasTopRing)
    {
        const int center = static_cast<int>(surfaceNodePositions_.size());
        surfaceNodePositions_.push_back(topCenter_);
        for (int index = 0; index < segmentCount; ++index)
        {
            surfaceTriangles_.push_back(triangle(center, topBegin + index, topBegin + (index + 1) % segmentCount));
        }
    }
}

bool ConeWall::isValid() const noexcept { return configured_; }

Vec3 ConeWall::boundingBoxMin() const noexcept
{
    const Vec3 axis = math::normalizedOrZero(topCenter_ - bottomCenter_);
    const Vec3 radialFactor{std::sqrt(std::max(0.0, 1.0 - axis.x * axis.x)), std::sqrt(std::max(0.0, 1.0 - axis.y * axis.y)), std::sqrt(std::max(0.0, 1.0 - axis.z * axis.z))};
    return math::componentMin(bottomCenter_ - bottomRadius_ * radialFactor, topCenter_ - topRadius_ * radialFactor);
}

Vec3 ConeWall::boundingBoxMax() const noexcept
{
    const Vec3 axis = math::normalizedOrZero(topCenter_ - bottomCenter_);
    const Vec3 radialFactor{std::sqrt(std::max(0.0, 1.0 - axis.x * axis.x)), std::sqrt(std::max(0.0, 1.0 - axis.y * axis.y)), std::sqrt(std::max(0.0, 1.0 - axis.z * axis.z))};
    return math::componentMax(bottomCenter_ + bottomRadius_ * radialFactor, topCenter_ + topRadius_ * radialFactor);
}

Real ConeWall::evaluateSFD(const Vec3& point) const noexcept
{
    const Vec3 endpointVector = topCenter_ - bottomCenter_;
    const Real height = math::norm(endpointVector);
    const Vec3 axis = endpointVector / height;
    const Vec3 relative = point - 0.5 * (bottomCenter_ + topCenter_);
    const Real axial = math::dot(relative, axis);
    const Real radial = math::norm(relative - axial * axis);
    const Real halfHeight = 0.5 * height;
    const Real capRadius = axial < 0.0 ? bottomRadius_ : topRadius_;
    const Real capX = radial - std::min(radial, capRadius);
    const Real capY = std::abs(axial) - halfHeight;
    const Real slopeX = topRadius_ - bottomRadius_;
    const Real slopeY = height;
    const Real fraction = std::clamp(((topRadius_ - radial) * slopeX + (halfHeight - axial) * slopeY) / (slopeX * slopeX + slopeY * slopeY), 0.0, 1.0);
    const Real sideX = radial - topRadius_ + slopeX * fraction;
    const Real sideY = axial - halfHeight + slopeY * fraction;
    const Real sign = sideX < 0.0 && capY < 0.0 ? -1.0 : 1.0;
    return sign * std::sqrt(std::min(capX * capX + capY * capY, sideX * sideX + sideY * sideY));
}

} // namespace fundem::levelset
