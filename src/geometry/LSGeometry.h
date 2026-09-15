/**
 * @file LSGeometry.h
 * @brief Defines reusable level-set geometry records and their combined container.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "math/Indexing.h"
#include "math/Matrix3.h"
#include "math/Vector3.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace fundem
{

/** One descriptor per reusable level-set geometry. */
class LSGeometryDescriptor
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;
    using Mat3 = math::Mat3;

    Vec3 gridNodeOrigin_{Vec3::zero()};           ///< First grid-node position in the centered local frame.
    Real gridNodeInverseSpacing_{0.0};            ///< Reciprocal uniform grid spacing.
    int3 gridNodeSize_{0, 0, 0};                  ///< Node counts along x, y, and z.
    int signedDistanceOffset_{0};                 ///< First signed-distance sample in the combined node array.
    int surfaceNodeOffset_{0};                    ///< First node in the combined surface-node array.
    int surfaceTriangleOffset_{0};                ///< First triangle in the combined connectivity array.
    Real boundingRadius_{0.0};                    ///< Maximum local surface-node radius.
    Real volume_{0.0};                            ///< Integrated solid volume; zero for fixed geometry.
    Mat3 unitDensityInertiaTensor_{Mat3::zero()}; ///< Centroidal inertia tensor at unit density.

    /** Read-only structure-of-arrays view passed to device algorithms. */
    struct const_device_type {
        const Vec3* gridNodeOrigin_{nullptr};           ///< First local grid-node positions.
        const Real* gridNodeInverseSpacing_{nullptr};   ///< Reciprocal grid spacings.
        const int3* gridNodeSize_{nullptr};             ///< Grid dimensions.
        const int* signedDistanceOffset_{nullptr};      ///< Signed-distance array offsets.
        const int* surfaceNodeOffset_{nullptr};         ///< Surface-node array offsets.
        const int* surfaceTriangleOffset_{nullptr};     ///< Surface-triangle array offsets.
        const Real* boundingRadius_{nullptr};           ///< Local bounding radii.
        const Real* volume_{nullptr};                   ///< Integrated geometry volumes.
        const Mat3* unitDensityInertiaTensor_{nullptr}; ///< Unit-density inertia tensors.
    };

    struct gridNodeOriginField {
        inline static constexpr auto member = &LSGeometryDescriptor::gridNodeOrigin_;
    };
    struct gridNodeInverseSpacingField {
        inline static constexpr auto member = &LSGeometryDescriptor::gridNodeInverseSpacing_;
    };
    struct gridNodeSizeField {
        inline static constexpr auto member = &LSGeometryDescriptor::gridNodeSize_;
    };
    struct signedDistanceOffsetField {
        inline static constexpr auto member = &LSGeometryDescriptor::signedDistanceOffset_;
    };
    struct surfaceNodeOffsetField {
        inline static constexpr auto member = &LSGeometryDescriptor::surfaceNodeOffset_;
    };
    struct surfaceTriangleOffsetField {
        inline static constexpr auto member = &LSGeometryDescriptor::surfaceTriangleOffset_;
    };
    struct boundingRadiusField {
        inline static constexpr auto member = &LSGeometryDescriptor::boundingRadius_;
    };
    struct volumeField {
        inline static constexpr auto member = &LSGeometryDescriptor::volume_;
    };
    struct unitDensityInertiaTensorField {
        inline static constexpr auto member = &LSGeometryDescriptor::unitDensityInertiaTensor_;
    };

    using DeviceLayout = deviceLayout<gridNodeOriginField,
                                      gridNodeInverseSpacingField,
                                      gridNodeSizeField,
                                      signedDistanceOffsetField,
                                      surfaceNodeOffsetField,
                                      surfaceTriangleOffsetField,
                                      boundingRadiusField,
                                      volumeField,
                                      unitDensityInertiaTensorField>;
};

/** One entry per signed-distance grid node. */
class LSGridNode
{
public:
    using Real = math::Real;

    Real signedDistance_{0.0}; ///< Signed distance at one Cartesian grid node.

    /** Read-only signed-distance array view. */
    struct const_device_type {
        const Real* signedDistance_{nullptr}; ///< Flattened signed-distance samples.
    };

    struct signedDistanceField {
        inline static constexpr auto member = &LSGridNode::signedDistance_;
    };

    using DeviceLayout = deviceLayout<signedDistanceField>;
};

/** One entry per level-set surface node. */
class LSSurfaceNode
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    Vec3 position_{Vec3::zero()}; ///< Geometry-local surface-node position.
    Real area_{0.0};              ///< Tributary surface area represented by this node.

    /** Read-only surface-node array view. */
    struct const_device_type {
        const Vec3* position_{nullptr}; ///< Geometry-local positions.
        const Real* area_{nullptr};     ///< Tributary nodal areas.
    };

    struct positionField {
        inline static constexpr auto member = &LSSurfaceNode::position_;
    };
    struct areaField {
        inline static constexpr auto member = &LSSurfaceNode::area_;
    };

    using DeviceLayout = deviceLayout<positionField, areaField>;
};

/** One entry per level-set surface triangle. */
class LSSurfaceTriangle
{
public:
    int3 connectivity_{-1, -1, -1}; ///< Three local surface-node indices.

    /** Read-only triangle-connectivity array view. */
    struct const_device_type {
        const int3* connectivity_{nullptr}; ///< Three local node indices per triangle.
    };

    struct connectivityField {
        inline static constexpr auto member = &LSSurfaceTriangle::connectivity_;
    };

    using DeviceLayout = deviceLayout<connectivityField>;
};

using LSGeometryDescriptorContainer = hostAoSDeviceSoA<LSGeometryDescriptor, LSGeometryDescriptor::DeviceLayout>;
using LSGridNodeContainer = hostAoSDeviceSoA<LSGridNode, LSGridNode::DeviceLayout>;
using LSSurfaceNodeContainer = hostAoSDeviceSoA<LSSurfaceNode, LSSurfaceNode::DeviceLayout>;
using LSSurfaceTriangleContainer = hostAoSDeviceSoA<LSSurfaceTriangle, LSSurfaceTriangle::DeviceLayout>;

/**
 * Owns the four combined arrays that describe all reusable level-set geometries.
 * Geometry indices and offsets remain stable because insertion is append-only.
 */
class LSGeometryContainer
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;
    using Mat3 = math::Mat3;

    LSGeometryContainer() = default;
    LSGeometryContainer(const LSGeometryContainer&) = delete;
    LSGeometryContainer& operator=(const LSGeometryContainer&) = delete;
    LSGeometryContainer(LSGeometryContainer&&) noexcept = default;
    LSGeometryContainer& operator=(LSGeometryContainer&&) noexcept = default;

    const std::vector<LSGeometryDescriptor>& hostDescriptors() const noexcept { return descriptors_.host(); }

    const std::vector<LSGridNode>& hostGridNodes() const noexcept { return gridNodes_.host(); }

    const std::vector<LSSurfaceNode>& hostSurfaceNodes() const noexcept { return surfaceNodes_.host(); }

    const std::vector<LSSurfaceTriangle>& hostSurfaceTriangles() const noexcept { return surfaceTriangles_.host(); }

    /**
     * Validates and appends one signed-distance geometry.
     * @param gridNodeOrigin Position of grid node (0,0,0) in the input frame.
     * @param gridNodeSpacing Positive uniform grid spacing.
     * @param gridNodeSize Number of nodes along each axis.
     * @param gridNodeSignedDistance Flattened signed-distance samples.
     * @param surfaceNodePositions Input-frame surface-node positions.
     * @param connectivities Triangle indices into @p surfaceNodePositions.
     * @param isFixed Skip centroid, volume, and inertia integration when true.
     * @return Stable descriptor index of the appended geometry.
     * @throws std::invalid_argument for malformed input arrays.
     * @throws std::domain_error when a movable geometry has no usable interior
     * volume or a surface node has no positive area.
     */
    int add(const Vec3& gridNodeOrigin,
            Real gridNodeSpacing,
            const int3& gridNodeSize,
            const std::vector<Real>& gridNodeSignedDistance,
            const std::vector<Vec3>& surfaceNodePositions,
            const std::vector<int3>& connectivities,
            bool isFixed = false)
    {
        if (gridNodeSize.x < 2 || gridNodeSize.y < 2 || gridNodeSize.z < 2)
        {
            throw std::invalid_argument("Invalid level-set grid size.");
        }

        const int gridNodeCount = gridNodeSize.x * gridNodeSize.y * gridNodeSize.z;
        if (!math::isFinite(gridNodeOrigin) || !math::isFinite(gridNodeSpacing) || gridNodeSpacing <= 0.0 || gridNodeCount != static_cast<int>(gridNodeSignedDistance.size()))
        {
            throw std::invalid_argument("Invalid level-set grid definition.");
        }

        for (Real signedDistance : gridNodeSignedDistance)
        {
            if (!math::isFinite(signedDistance))
            {
                throw std::invalid_argument("Level-set signed distances must be finite.");
            }
        }

        if (surfaceNodePositions.empty())
        {
            throw std::invalid_argument("Invalid level-set surface storage size.");
        }
        for (const Vec3& position : surfaceNodePositions)
        {
            if (!math::isFinite(position))
            {
                throw std::invalid_argument("Level-set surface node positions must be finite.");
            }
        }

        const int surfaceNodeCount = static_cast<int>(surfaceNodePositions.size());
        for (const int3& connectivity : connectivities)
        {
            if (connectivity.x < 0 || connectivity.x >= surfaceNodeCount || connectivity.y < 0 || connectivity.y >= surfaceNodeCount || connectivity.z < 0 || connectivity.z >= surfaceNodeCount)
            {
                throw std::invalid_argument("Level-set triangle connectivity is invalid.");
            }
        }

        const Real inverseSpacing = 1.0 / gridNodeSpacing;
        Real volume = 0.0;
        Vec3 centroid = Vec3::zero();
        Mat3 unitDensityInertiaTensor = Mat3::zero();

        if (!isFixed)
        {
            const auto interiorFraction = [](Real dimensionlessDistance) noexcept
            {
                constexpr Real smoothingWidth = 1.5;
                if (dimensionlessDistance < -smoothingWidth)
                {
                    return Real{1.0};
                }
                if (dimensionlessDistance > smoothingWidth)
                {
                    return Real{0.0};
                }

                const Real normalizedDistance = -dimensionlessDistance / smoothingWidth;
                return 0.5 * (1.0 + normalizedDistance + std::sin(math::pi * normalizedDistance) / math::pi);
            };

            const int sizeX = gridNodeSize.x;
            const int sizeY = gridNodeSize.y;
            const int sizeZ = gridNodeSize.z;
            const Real nodeVolume = gridNodeSpacing * gridNodeSpacing * gridNodeSpacing;

            Real occupancySum = 0.0;
            Vec3 firstMoment = Vec3::zero();
            for (int z = 0; z < sizeZ; ++z)
            {
                for (int y = 0; y < sizeY; ++y)
                {
                    for (int x = 0; x < sizeX; ++x)
                    {
                        const int index = math::linearIndex(x, y, z, sizeX, sizeY);
                        const Real occupancy = interiorFraction(gridNodeSignedDistance[index] * inverseSpacing);
                        const Vec3 nodePosition = gridNodeOrigin + gridNodeSpacing * Vec3{Real(x), Real(y), Real(z)};
                        occupancySum += occupancy;
                        firstMoment += occupancy * nodePosition;
                    }
                }
            }

            volume = occupancySum * nodeVolume;
            if (!math::isFinite(occupancySum) || occupancySum <= math::defaultTolerance || !math::isFinite(volume) || volume <= math::defaultTolerance)
            {
                throw std::domain_error("Level-set grid has no finite interior volume.");
            }

            centroid = firstMoment / occupancySum;
            if (!math::isFinite(centroid))
            {
                throw std::overflow_error("Level-set centroid integration is not finite.");
            }

            for (int z = 0; z < sizeZ; ++z)
            {
                for (int y = 0; y < sizeY; ++y)
                {
                    for (int x = 0; x < sizeX; ++x)
                    {
                        const int index = math::linearIndex(x, y, z, sizeX, sizeY);
                        const Real occupancy = interiorFraction(gridNodeSignedDistance[index] * inverseSpacing);
                        const Real pointMass = occupancy * nodeVolume;
                        const Vec3 pointPosition = gridNodeOrigin + gridNodeSpacing * Vec3{Real(x), Real(y), Real(z)} - centroid;

                        unitDensityInertiaTensor(0, 0) += pointMass * (pointPosition.y * pointPosition.y + pointPosition.z * pointPosition.z);
                        unitDensityInertiaTensor(1, 1) += pointMass * (pointPosition.x * pointPosition.x + pointPosition.z * pointPosition.z);
                        unitDensityInertiaTensor(2, 2) += pointMass * (pointPosition.x * pointPosition.x + pointPosition.y * pointPosition.y);
                        unitDensityInertiaTensor(0, 1) -= pointMass * pointPosition.x * pointPosition.y;
                        unitDensityInertiaTensor(0, 2) -= pointMass * pointPosition.x * pointPosition.z;
                        unitDensityInertiaTensor(1, 2) -= pointMass * pointPosition.y * pointPosition.z;
                    }
                }
            }
            unitDensityInertiaTensor(1, 0) = unitDensityInertiaTensor(0, 1);
            unitDensityInertiaTensor(2, 0) = unitDensityInertiaTensor(0, 2);
            unitDensityInertiaTensor(2, 1) = unitDensityInertiaTensor(1, 2);
            if (!math::isFinite(unitDensityInertiaTensor))
            {
                throw std::overflow_error("Level-set inertia integration is not finite.");
            }
        }

        std::vector<LSGridNode> centeredGridNodes;
        centeredGridNodes.reserve(gridNodeSignedDistance.size());
        for (Real signedDistance : gridNodeSignedDistance)
        {
            centeredGridNodes.push_back({signedDistance});
        }

        std::vector<LSSurfaceNode> centeredSurfaceNodes;
        centeredSurfaceNodes.reserve(surfaceNodePositions.size());
        for (const Vec3& position : surfaceNodePositions)
        {
            centeredSurfaceNodes.push_back({position, 0.0});
        }
        Real boundingRadius = 0.0;
        for (LSSurfaceNode& node : centeredSurfaceNodes)
        {
            if (!isFixed)
            {
                node.position_ -= centroid;
            }
            node.area_ = 0.0;
            const Real candidate = math::norm(node.position_);
            if (!math::isFinite(candidate))
            {
                throw std::overflow_error("Level-set bounding radius is not finite.");
            }
            boundingRadius = candidate > boundingRadius ? candidate : boundingRadius;
        }
        if (boundingRadius <= 0.0)
        {
            throw std::domain_error("Level-set bounding radius must be positive.");
        }

        if (connectivities.empty())
        {
            const Real fallbackArea = 4.0 * math::pi * boundingRadius * boundingRadius / surfaceNodeCount;
            for (LSSurfaceNode& node : centeredSurfaceNodes)
            {
                node.area_ = fallbackArea;
            }
        }
        else
        {
            for (const int3& connectivity : connectivities)
            {
                const Vec3& a = centeredSurfaceNodes[connectivity.x].position_;
                const Vec3& b = centeredSurfaceNodes[connectivity.y].position_;
                const Vec3& c = centeredSurfaceNodes[connectivity.z].position_;
                const Real nodalArea = math::norm(math::cross(b - a, c - a)) / 6.0;

                centeredSurfaceNodes[connectivity.x].area_ += nodalArea;
                centeredSurfaceNodes[connectivity.y].area_ += nodalArea;
                centeredSurfaceNodes[connectivity.z].area_ += nodalArea;
            }
        }

        for (const LSSurfaceNode& node : centeredSurfaceNodes)
        {
            if (!math::isFinite(node.area_) || node.area_ <= 0.0)
            {
                throw std::domain_error("Level-set surface node areas must be finite and positive.");
            }
        }

        auto& descriptorHost = descriptors_.host();
        auto& gridNodeHost = gridNodes_.host();
        auto& surfaceNodeHost = surfaceNodes_.host();
        auto& surfaceTriangleHost = surfaceTriangles_.host();

        const int geometryIndex = static_cast<int>(descriptorHost.size());
        const int signedDistanceOffset = static_cast<int>(gridNodeHost.size());
        const int surfaceNodeOffset = static_cast<int>(surfaceNodeHost.size());
        const int surfaceTriangleOffset = static_cast<int>(surfaceTriangleHost.size());

        descriptorHost.reserve(descriptorHost.size() + 1);
        gridNodeHost.reserve(gridNodeHost.size() + centeredGridNodes.size());
        surfaceNodeHost.reserve(surfaceNodeHost.size() + centeredSurfaceNodes.size());
        surfaceTriangleHost.reserve(surfaceTriangleHost.size() + connectivities.size());

        gridNodeHost.insert(gridNodeHost.end(), centeredGridNodes.begin(), centeredGridNodes.end());
        surfaceNodeHost.insert(surfaceNodeHost.end(), centeredSurfaceNodes.begin(), centeredSurfaceNodes.end());
        for (const int3& connectivity : connectivities)
        {
            surfaceTriangleHost.push_back({connectivity});
        }
        descriptorHost.push_back({isFixed ? gridNodeOrigin : gridNodeOrigin - centroid,
                                  inverseSpacing,
                                  gridNodeSize,
                                  signedDistanceOffset,
                                  surfaceNodeOffset,
                                  surfaceTriangleOffset,
                                  boundingRadius,
                                  volume,
                                  unitDensityInertiaTensor});

        return geometryIndex;
    }

    LSGeometryDescriptor::const_device_type descriptors() const noexcept
    {
        return {descriptors_.device<LSGeometryDescriptor::gridNodeOriginField>(),
                descriptors_.device<LSGeometryDescriptor::gridNodeInverseSpacingField>(),
                descriptors_.device<LSGeometryDescriptor::gridNodeSizeField>(),
                descriptors_.device<LSGeometryDescriptor::signedDistanceOffsetField>(),
                descriptors_.device<LSGeometryDescriptor::surfaceNodeOffsetField>(),
                descriptors_.device<LSGeometryDescriptor::surfaceTriangleOffsetField>(),
                descriptors_.device<LSGeometryDescriptor::boundingRadiusField>(),
                descriptors_.device<LSGeometryDescriptor::volumeField>(),
                descriptors_.device<LSGeometryDescriptor::unitDensityInertiaTensorField>()};
    }

    LSGridNode::const_device_type gridNodes() const noexcept { return {gridNodes_.device<LSGridNode::signedDistanceField>()}; }

    LSSurfaceNode::const_device_type surfaceNodes() const noexcept { return {surfaceNodes_.device<LSSurfaceNode::positionField>(), surfaceNodes_.device<LSSurfaceNode::areaField>()}; }

    LSSurfaceTriangle::const_device_type surfaceTriangles() const noexcept { return {surfaceTriangles_.device<LSSurfaceTriangle::connectivityField>()}; }

    /** Returns allocated device memory for all four arrays in GiB. */
    double deviceMemoryGB() const noexcept { return descriptors_.deviceMemoryGB() + gridNodes_.deviceMemoryGB() + surfaceNodes_.deviceMemoryGB() + surfaceTriangles_.deviceMemoryGB(); }

    /**
     * Enqueues all four host-to-device copies on @p stream.
     * @param stream CUDA stream receiving the copy operations.
     */
    void copyHostToDeviceAsync(cudaStream_t stream = nullptr)
    {
        descriptors_.copyHostToDeviceAsync(stream);
        gridNodes_.copyHostToDeviceAsync(stream);
        surfaceNodes_.copyHostToDeviceAsync(stream);
        surfaceTriangles_.copyHostToDeviceAsync(stream);
    }

    /** Copies all geometry arrays to the device and waits for @p stream. */
    void copyHostToDevice(cudaStream_t stream = nullptr)
    {
        copyHostToDeviceAsync(stream);
        host_device_detail::synchronize(stream);
    }

private:
    LSGeometryDescriptorContainer descriptors_;   ///< Per-geometry metadata.
    LSGridNodeContainer gridNodes_;               ///< Combined signed-distance samples.
    LSSurfaceNodeContainer surfaceNodes_;         ///< Combined local surface nodes and areas.
    LSSurfaceTriangleContainer surfaceTriangles_; ///< Combined local triangle connectivity.
};

} // namespace fundem
