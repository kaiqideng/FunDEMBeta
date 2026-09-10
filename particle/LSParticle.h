/**
 * @file LSParticle.h
 * @brief Defines rigid level-set particles with host and device geometry state.
 */
#pragma once

#include "geometry/LSGeometry.h"
#include "particle.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace fundem
{

/**
 * Rigid level-set particle that stores a host geometry copy and a stable index
 * into the shared device geometry container.
 */
class LSParticle : public particle
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;
    using Mat3 = math::Mat3;
    using Quaternion = math::Quaternion;

    LSParticle() = default;
    ~LSParticle() override = default;
    LSParticle(const LSParticle&) = default;
    LSParticle& operator=(const LSParticle&) = default;
    LSParticle(LSParticle&&) noexcept = default;
    LSParticle& operator=(LSParticle&&) noexcept = default;

    /** Creates LS rigid kinematics before material and geometry are bound. */
    LSParticle(const Vec3& position, const Quaternion& orientation, const Vec3& velocity, const Vec3& angularVelocity) : particle(position, orientation, velocity, angularVelocity) {}

    Real boundingRadius() const noexcept { return particle::radius(); }
    const Vec3& gridNodeOrigin() const noexcept { return gridNodeOrigin_; }
    Real gridNodeInverseSpacing() const noexcept { return gridNodeInverseSpacing_; }
    const int3& gridNodeSize() const noexcept { return gridNodeSize_; }
    const std::vector<LSGridNode>& gridNodes() const noexcept { return gridNodes_; }
    const std::vector<LSSurfaceNode>& surfaceNodes() const noexcept { return surfaceNodes_; }
    const std::vector<LSSurfaceTriangle>& surfaceTriangles() const noexcept { return surfaceTriangles_; }
    Real volume() const noexcept { return volume_; }
    const Mat3& unitDensityInertiaTensor() const noexcept { return unitDensityInertiaTensor_; }
    int geometryIndex() const noexcept { return geometryIndex_; }
    /** Tests whether the stored stable geometry index belongs to @p container. */
    bool geometryBelongsTo(const LSGeometryContainer& container) const noexcept { return geometryContainer_ == &container; }

    /** Binds only level-set material definitions; standard materials are rejected. */
    void setMaterial(const materialContainer& container, int index) noexcept
    {
        if (index < 0 || index >= static_cast<int>(container.hostSize()) || !container.host()[index].isLevelSet())
        {
            return;
        }
        particle::setMaterial(container, index);
    }

    /**
     * Copies the selected geometry into host-local state and records its stable
     * shared-container index for device execution.
     * @param container Source geometry container.
     * @param index Descriptor index in @p container.
     */
    void setGeometry(const LSGeometryContainer& container, int index)
    {
        const auto& descriptors = container.hostDescriptors();
        const auto& gridNodes = container.hostGridNodes();
        const auto& surfaceNodes = container.hostSurfaceNodes();
        const auto& surfaceTriangles = container.hostSurfaceTriangles();

        const int geometryCount = static_cast<int>(descriptors.size());
        if (index < 0 || index >= geometryCount)
        {
            return;
        }

        const LSGeometryDescriptor& descriptor = descriptors[index];
        const int gridNodeCount = descriptor.gridNodeSize_.x * descriptor.gridNodeSize_.y * descriptor.gridNodeSize_.z;
        const int gridNodeBegin = descriptor.signedDistanceOffset_;
        const int gridNodeEnd = gridNodeBegin + gridNodeCount;
        const int surfaceNodeBegin = descriptor.surfaceNodeOffset_;
        const int surfaceNodeEnd = index + 1 < geometryCount ? descriptors[index + 1].surfaceNodeOffset_ : static_cast<int>(surfaceNodes.size());
        const int surfaceTriangleBegin = descriptor.surfaceTriangleOffset_;
        const int surfaceTriangleEnd = index + 1 < geometryCount ? descriptors[index + 1].surfaceTriangleOffset_ : static_cast<int>(surfaceTriangles.size());

        std::vector<LSGridNode> newGridNodes(gridNodes.begin() + gridNodeBegin, gridNodes.begin() + gridNodeEnd);

        std::vector<LSSurfaceNode> newSurfaceNodes(surfaceNodes.begin() + surfaceNodeBegin, surfaceNodes.begin() + surfaceNodeEnd);
        std::vector<LSSurfaceTriangle> newSurfaceTriangles(surfaceTriangles.begin() + surfaceTriangleBegin, surfaceTriangles.begin() + surfaceTriangleEnd);

        if (!setRadiusValue(descriptor.boundingRadius_))
        {
            return;
        }

        gridNodeOrigin_ = descriptor.gridNodeOrigin_;
        gridNodeInverseSpacing_ = descriptor.gridNodeInverseSpacing_;
        gridNodeSize_ = descriptor.gridNodeSize_;
        gridNodes_ = std::move(newGridNodes);
        surfaceNodes_ = std::move(newSurfaceNodes);
        surfaceTriangles_ = std::move(newSurfaceTriangles);
        volume_ = descriptor.volume_;
        unitDensityInertiaTensor_ = descriptor.unitDensityInertiaTensor_;
        geometryContainer_ = &container;
        geometryIndex_ = index;
        geometrySet_ = true;
        setMassAndInertia();
    }

    /** Reports whether particle state, level-set material, and geometry are ready for insertion. */
    bool isValid() const noexcept { return particle::isValid() && geometrySet_; }

    struct geometryIndexField;

    /** Mutable device view extending sphere state with geometry indices. */
    struct device_type : particle::device_type {
        int* geometryIndex_{nullptr}; ///< Indices into `LSGeometryContainer`.
    };

private:
    /** Recomputes mass and inertia from geometry integrals and material density. */
    virtual void setMassAndInertia() noexcept override
    {
        const Real density = getMaterial().density();
        const Real mass = density * volume_;
        const Mat3 inertiaTensor = density * unitDensityInertiaTensor_;
        rigidBody::setMassAndInertia(mass, inertiaTensor);
    }

    using particle::isValid;
    using particle::radius;
    using particle::setRadius;

    Vec3 gridNodeOrigin_{Vec3::zero()};               ///< Local position of grid node (0,0,0).
    Real gridNodeInverseSpacing_{0.0};                ///< Reciprocal Cartesian grid spacing.
    int3 gridNodeSize_{0, 0, 0};                      ///< Grid-node counts along each axis.
    Real volume_{0.0};                                ///< Integrated geometry volume.
    Mat3 unitDensityInertiaTensor_{Mat3::zero()};     ///< Centroidal inertia at unit density.
    std::vector<LSGridNode> gridNodes_;               ///< Host-local signed-distance grid nodes.
    std::vector<LSSurfaceNode> surfaceNodes_;         ///< Host-local surface nodes and areas.
    std::vector<LSSurfaceTriangle> surfaceTriangles_; ///< Host-local triangle connectivity.

    bool geometrySet_{false};                               ///< Whether geometry binding succeeded.
    const LSGeometryContainer* geometryContainer_{nullptr}; ///< Owning shared geometry container.
    int geometryIndex_{-1};                                 ///< Stable shared geometry index.

public:
    struct geometryIndexField {
        inline static constexpr auto member = &LSParticle::geometryIndex_;
    };

    using DeviceLayout = concatDeviceLayoutsT<particle::DeviceLayout, deviceLayout<geometryIndexField>>;
};

using LSParticleContainer = hostAoSDeviceSoA<LSParticle, LSParticle::DeviceLayout>;

/** Builds the mutable device view of a level-set particle container. */
inline LSParticle::device_type deviceFields(LSParticleContainer& particles) noexcept
{
    LSParticle::device_type result;
    static_cast<particle::device_type&>(result) = particleDeviceFields(particles);
    result.geometryIndex_ = particles.device<LSParticle::geometryIndexField>();
    return result;
}

} // namespace fundem
