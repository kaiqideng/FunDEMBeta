/**
 * @file contactDetectionKernel.cuh
 * @brief Declares streamed CUDA contact-detection and prefix-sum launch interfaces.
 */
#pragma once

#include "geometry/LSGeometry.h"
#include "interaction/interactionContainer.h"
#include "particle/spatialGrid.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{

using spatialGridDeviceView = spatialGridContainer::const_device_type;

namespace detail
{

/** Read-only rigid-particle arrays required by narrow-phase detection. */
struct particleDeviceView {
    const math::Vec3* position_{nullptr};          ///< World positions.
    const math::Quaternion* orientation_{nullptr}; ///< Unit orientations.
    const math::Real* inverseMass_{nullptr};       ///< Reciprocal masses.
    const math::Real* radius_{nullptr};            ///< Bounding radii.
};

/** Rigid state plus shared level-set geometry index. */
struct levelSetParticleDeviceView {
    particleDeviceView body_;           ///< Common rigid-particle arrays.
    const int* geometryIndex_{nullptr}; ///< Shared geometry descriptor indices.
};

/** Read-only level-set geometry arrays used by interpolation and patch detection. */
struct levelSetGeometryDeviceView {
    LSGeometryDescriptor::const_device_type descriptors_; ///< Per-geometry metadata.
    LSGridNode::const_device_type gridNodes_;             ///< Signed-distance samples.
    LSSurfaceNode::const_device_type surfaceNodes_;       ///< Local surface positions and areas.
};

/** Expanded LS surface-node mapping arrays. */
struct surfaceNodeMappingDeviceView {
    const int* geometrySurfaceNodeIndex_{nullptr}; ///< Combined geometry-node indices.
    const int* ownerParticleIndex_{nullptr};       ///< Owning LS-particle indices.
    const int* particleSurfaceNodeIndex_{nullptr}; ///< Owner-local surface-node indices.
};

/** Writable contact fields produced by detection kernels. */
struct contactOutputDeviceView {
    int* masterParticleIndex_{nullptr};               ///< Master indices.
    int* slaveParticleIndex_{nullptr};                ///< Slave indices.
    int* particleSurfaceNodeIndex_{nullptr};          ///< Owner-local LS node indices.
    math::Vec3* point_{nullptr};                      ///< World contact points.
    math::Vec3* normal_{nullptr};                     ///< Slave-to-master unit normals.
    math::Real* overlap_{nullptr};                    ///< Positive penetrations.
    math::Real* area_{nullptr};                       ///< Contact or patch areas.
    math::Real* effectiveMass_{nullptr};              ///< Damping effective masses.
    math::Real* effectiveRadius_{nullptr};            ///< Moment effective radii.
    math::Real* normalForceMagnitude_{nullptr};       ///< Reset normal-force magnitudes.
    math::Vec3* slidingSpringDeformation_{nullptr};   ///< Restored sliding histories.
    math::Vec3* rollingSpringDeformation_{nullptr};   ///< Restored rolling histories.
    math::Vec3* torsionalSpringDeformation_{nullptr}; ///< Restored torsional histories.
};

/** Read-only prior-step contact history and range-prefix arrays. */
struct contactHistoryDeviceView {
    const int* slaveParticleIndex_{nullptr};                ///< Prior slave indices.
    const math::Vec3* slidingSpringDeformation_{nullptr};   ///< Prior sliding states.
    const math::Vec3* rollingSpringDeformation_{nullptr};   ///< Prior rolling states.
    const math::Vec3* torsionalSpringDeformation_{nullptr}; ///< Prior torsional states.
    const int* prefixSum_{nullptr};                         ///< Prior owner-range ends.
};

} // namespace detail

/** Detects sphere-sphere contacts and updates interaction storage on @p stream. */
void launchSphereSphereContactDetection(interactionContainer& interactions, const particleContainer& spheres, const spatialGridDeviceView& spatialGrid, cudaStream_t stream = nullptr);

/** Detects sphere-LS contacts with LS particles ordered as slaves. */
void launchSphereLevelSetContactDetection(interactionContainer& interactions,
                                          const particleContainer& masterSpheres,
                                          const LSParticleContainer& slaveLSParticles,
                                          const LSGeometryContainer& geometries,
                                          const spatialGridDeviceView& slaveSpatialGrid,
                                          cudaStream_t stream = nullptr);

/** Detects LS-LS contacts using expanded surface-node owners. */
void launchLevelSetLevelSetContactDetection(interactionContainer& interactions,
                                            const LSParticleContainer& LSParticles,
                                            const LSGeometryContainer& geometries,
                                            const spatialGridDeviceView& spatialGrid,
                                            cudaStream_t stream = nullptr);

} // namespace fundem::cuda
