/**
 * @file surfaceNodeMapping.h
 * @brief Maps expanded LS-particle surface nodes to geometry nodes and owners.
 */
#pragma once

#include "particle/LSParticle.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace fundem
{

/** Maps one expanded particle surface node to its LSGeometryContainer node and owner. */
class surfaceNodeMapping
{
public:
    int geometrySurfaceNodeIndex_{-1}; ///< Index in the combined geometry surface-node array.
    int ownerParticleIndex_{-1};       ///< Owning LS-particle index.
    int particleSurfaceNodeIndex_{-1}; ///< Node index local to the expanded particle surface.

    struct geometrySurfaceNodeIndexField {
        inline static constexpr auto member = &surfaceNodeMapping::geometrySurfaceNodeIndex_;
    };
    struct ownerParticleIndexField {
        inline static constexpr auto member = &surfaceNodeMapping::ownerParticleIndex_;
    };
    struct particleSurfaceNodeIndexField {
        inline static constexpr auto member = &surfaceNodeMapping::particleSurfaceNodeIndex_;
    };

    using DeviceLayout = deviceLayout<geometrySurfaceNodeIndexField, ownerParticleIndexField, particleSurfaceNodeIndexField>;
};

using surfaceNodeMappingStorage = hostAoSDeviceSoA<surfaceNodeMapping, surfaceNodeMapping::DeviceLayout>;

/**
 * Append-order mapping storage generated from an LS-particle container.
 * Each mapping index is also the owner range used by LS contact detection.
 */
class surfaceNodeMappingContainer : public surfaceNodeMappingStorage
{
public:
    surfaceNodeMappingContainer() = default;
    /** Generates the complete append-order mapping at construction. */
    surfaceNodeMappingContainer(const LSParticleContainer& particles, const LSGeometryContainer& geometries) { generate(particles, geometries); }

    /**
     * Rebuilds mappings from host LS-particle geometry bindings.
     * @throws std::invalid_argument if any particle lacks valid geometry.
     * Existing device storage is released because mapping order may change.
     */
    void generate(const LSParticleContainer& particles, const LSGeometryContainer& geometries)
    {
        std::vector<surfaceNodeMapping> newMappings;
        const auto& particleHost = particles.host();
        const auto& descriptors = geometries.hostDescriptors();
        const int geometryCount = static_cast<int>(descriptors.size());

        for (int particleIndex = 0; particleIndex < static_cast<int>(particleHost.size()); ++particleIndex)
        {
            const LSParticle& particle = particleHost[particleIndex];
            const int geometryIndex = particle.geometryIndex();
            if (geometryIndex < 0 || geometryIndex >= geometryCount)
            {
                throw std::invalid_argument("Every LSParticle must have valid level-set geometry before generating its surface-node mapping.");
            }

            const int surfaceNodeBegin = descriptors[geometryIndex].surfaceNodeOffset_;
            const int surfaceNodeEnd = geometryIndex + 1 < geometryCount ? descriptors[geometryIndex + 1].surfaceNodeOffset_ : static_cast<int>(geometries.hostSurfaceNodes().size());
            for (int surfaceNodeIndex = surfaceNodeBegin; surfaceNodeIndex < surfaceNodeEnd; ++surfaceNodeIndex)
            {
                newMappings.push_back({surfaceNodeIndex, particleIndex, surfaceNodeIndex - surfaceNodeBegin});
            }
        }

        resetDevice();
        host() = std::move(newMappings);
    }

    /** Enqueues a complete mapping upload on @p stream. */
    void copyHostToDeviceAsync(cudaStream_t stream = nullptr) { surfaceNodeMappingStorage::copyHostToDeviceAsync(stream); }

    /** Uploads mappings and waits for @p stream. */
    void copyHostToDevice(cudaStream_t stream = nullptr)
    {
        copyHostToDeviceAsync(stream);
        host_device_detail::synchronize(stream);
    }

    int nodeCount() const noexcept { return static_cast<int>(hostSize()); }

    /**
     * Finds the flattened mapping index for one owner-local node pair.
     * @return Mapping index, or `-1` when no mapping matches.
     * @note This linear host lookup is used during initialization, not inner
     * contact loops.
     */
    int mappingIndex(int ownerParticleIndex, int particleSurfaceNodeIndex) const noexcept
    {
        if (ownerParticleIndex < 0 || particleSurfaceNodeIndex < 0)
        {
            return -1;
        }

        for (int index = 0; index < static_cast<int>(hostSize()); ++index)
        {
            const surfaceNodeMapping& mapping = host()[index];
            if (mapping.ownerParticleIndex_ == ownerParticleIndex && mapping.particleSurfaceNodeIndex_ == particleSurfaceNodeIndex)
            {
                return index;
            }
        }
        return -1;
    }
};

} // namespace fundem
