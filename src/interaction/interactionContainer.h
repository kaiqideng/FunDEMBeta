/**
 * @file interactionContainer.h
 * @brief Aggregates contact, bond, neighbor, history, and mapping storage.
 */
#pragma once

#include "bond.h"
#include "contact.h"
#include "neighborRange.h"
#include "surfaceNodeMapping.h"
#include "particle/LSParticle.h"
#include "particle/particle.h"

#include <vector>

namespace fundem
{

/**
 * Composite owner for contacts, bonds, search ranges, candidate neighbors,
 * persistent histories, and expanded LS surface-node mappings.
 */
class interactionContainer
{
public:
    interactionContainer() = default;
    interactionContainer(const interactionContainer&) = delete;
    interactionContainer& operator=(const interactionContainer&) = delete;
    interactionContainer(interactionContainer&&) noexcept = default;
    interactionContainer& operator=(interactionContainer&&) noexcept = default;

    contactContainer& contacts() noexcept { return contacts_; }
    const contactContainer& contacts() const noexcept { return contacts_; }

    bondContainer& bonds() noexcept { return bonds_; }
    const bondContainer& bonds() const noexcept { return bonds_; }

    const surfaceNodeMappingContainer& surfaceNodeMappings() const noexcept { return surfaceNodeMappings_; }

    neighborRange::device_type particleNeighborRanges() noexcept
    {
        return {particleNeighborRanges_.device<neighborRange::countField>(), particleNeighborRanges_.device<neighborRange::prefixSumField>()};
    }

    neighborRange::device_type surfaceNodeNeighborRanges() noexcept
    {
        return {surfaceNodeNeighborRanges_.device<neighborRange::countField>(), surfaceNodeNeighborRanges_.device<neighborRange::prefixSumField>()};
    }

    particleNeighbor::device_type particleNeighborList() noexcept
    {
        return {particleNeighbors_.device<particleNeighbor::particleIndexField>(), particleNeighbors_.device<particleNeighbor::patchAreaField>()};
    }

    contactHistory::device_type contactHistories() noexcept
    {
        return {contactHistories_.device<contactHistory::slaveParticleIndexField>(),
                contactHistories_.device<contactHistory::slidingSpringDeformationField>(),
                contactHistories_.device<contactHistory::rollingSpringDeformationField>(),
                contactHistories_.device<contactHistory::torsionalSpringDeformationField>()};
    }

    neighborRangeHistory::device_type neighborRangeHistories() noexcept { return {neighborRangeHistories_.device<neighborRangeHistory::prefixSumField>()}; }

    /** Returns allocated device memory across every owned subcontainer in GiB. */
    double deviceMemoryGB() const noexcept
    {
        return contacts_.deviceMemoryGB() + bonds_.deviceMemoryGB() + surfaceNodeMappings_.deviceMemoryGB() + particleNeighborRanges_.deviceMemoryGB() + surfaceNodeNeighborRanges_.deviceMemoryGB() +
               particleNeighbors_.deviceMemoryGB() + contactHistories_.deviceMemoryGB() + neighborRangeHistories_.deviceMemoryGB();
    }

    /**
     * Prepares sphere interaction storage for a new device solve.
     * Existing contacts and histories are preserved when particle count does not
     * decrease; ranges are resized and bonds are uploaded.
     * @param particles Sphere storage defining the required range count.
     * @param stream CUDA stream used by every allocation copy and synchronization.
     */
    void initializeDevice(const particleContainer& particles, cudaStream_t stream = nullptr)
    {
        const int particleCount = static_cast<int>(particles.hostSize());
        initializeStorage(particleCount, 0, false, stream);
        surfaceNodeMappings_ = surfaceNodeMappingContainer{};
    }

    /**
     * Prepares LS interaction storage, rebuilds expanded surface mappings, and
     * uploads all device-side ranges and bonds.
     * @param particles LS storage defining particles and expanded surface nodes.
     * @param geometries Geometry container referenced by @p particles.
     * @param stream CUDA stream used by every device operation.
     */
    void initializeDevice(const LSParticleContainer& particles, const LSGeometryContainer& geometries, cudaStream_t stream = nullptr)
    {
        const int particleCount = static_cast<int>(particles.hostSize());
        surfaceNodeMappings_.generate(particles, geometries);
        surfaceNodeMappings_.copyHostToDeviceAsync(stream);
        initializeStorage(particleCount, surfaceNodeMappings_.nodeCount(), true, stream);
    }

    /**
     * Prefix-sums the active range counts, reads the final total, and resizes the
     * logical contact array to the exact number produced by detection.
     * @param stream CUDA stream used for prefix sum and scalar transfer.
     * @return Exact current contact count.
     */
    int contactCount(cudaStream_t stream = nullptr)
    {
#if FUNDEM_HAS_CUDA
        neighborRangeContainer& ranges = activeNeighborRanges();
        const int appliedCount = static_cast<int>(ranges.deviceSize());
        if (appliedCount == 0)
        {
            contacts_.resizeDevice(0);
            return 0;
        }

        buildPrefixSum(ranges, stream);
        int count = 0;
        host_device_detail::checkCuda(cudaMemcpyAsync(&count, ranges.device<neighborRange::prefixSumField>() + appliedCount - 1, sizeof(int), cudaMemcpyDeviceToHost, stream),
                                      "cudaMemcpyAsync contact count");
        host_device_detail::synchronize(stream);
        contacts_.resizeDevice(count);
        return count;
#else
        (void)stream;
        host_device_detail::throwCUDAUnavailable();
#endif
    }

    /**
     * Prefix-sums particle-neighbor counts and resizes candidate-neighbor
     * storage to the resulting exact total.
     * @return Exact candidate-neighbor count.
     */
    int resizeParticleNeighborList(cudaStream_t stream = nullptr)
    {
#if FUNDEM_HAS_CUDA
        const int particleCount = static_cast<int>(particleNeighborRanges_.deviceSize());
        if (particleCount == 0)
        {
            particleNeighbors_.resizeDevice(0);
            return 0;
        }

        buildPrefixSum(particleNeighborRanges_, stream);
        int neighborCount = 0;
        host_device_detail::checkCuda(cudaMemcpyAsync(&neighborCount, particleNeighborRanges_.device<neighborRange::prefixSumField>() + particleCount - 1, sizeof(int), cudaMemcpyDeviceToHost, stream),
                                      "cudaMemcpyAsync particle neighbor count");
        host_device_detail::synchronize(stream);
        particleNeighbors_.resizeDevice(neighborCount);
        return neighborCount;
#else
        (void)stream;
        host_device_detail::throwCUDAUnavailable();
#endif
    }

    /**
     * Copies current contact identities, spring states, and active range prefix
     * sums into next-step history storage.
     * The first call after host-history restoration intentionally preserves that
     * restored state.
     */
    void saveCurrentStepToHistory(cudaStream_t stream = nullptr)
    {
        if (preserveInitializedContactHistory_)
        {
            preserveInitializedContactHistory_ = false;
            return;
        }

        const neighborRangeContainer& ranges = activeNeighborRanges();
        const int appliedCount = static_cast<int>(ranges.deviceSize());
        const int currentContactCount = static_cast<int>(contacts_.deviceSize());

        contactHistories_.resizeDevice(currentContactCount);
        neighborRangeHistories_.resizeDevice(appliedCount);

        copyDeviceArrayAsync(neighborRangeHistories_.device<fundem::neighborRangeHistory::prefixSumField>(), ranges.device<neighborRange::prefixSumField>(), appliedCount, stream);
        copyDeviceArrayAsync(contactHistories_.device<fundem::contactHistory::slaveParticleIndexField>(), contacts_.device<contact::slaveParticleIndexField>(), currentContactCount, stream);
        copyDeviceArrayAsync(contactHistories_.device<fundem::contactHistory::slidingSpringDeformationField>(),
                             contacts_.device<contact::slidingSpringDeformationField>(),
                             currentContactCount,
                             stream);
        copyDeviceArrayAsync(contactHistories_.device<fundem::contactHistory::rollingSpringDeformationField>(),
                             contacts_.device<contact::rollingSpringDeformationField>(),
                             currentContactCount,
                             stream);
        copyDeviceArrayAsync(contactHistories_.device<fundem::contactHistory::torsionalSpringDeformationField>(),
                             contacts_.device<contact::torsionalSpringDeformationField>(),
                             currentContactCount,
                             stream);
    }

    /** Enqueues device-to-host copies for contacts and bonds only. */
    void copyDeviceToHostAsync(cudaStream_t stream = nullptr)
    {
        contacts_.copyDeviceToHostAsync(stream);
        bonds_.copyDeviceToHostAsync(stream);
    }

    /** Copies contacts and bonds to the host and waits for @p stream. */
    void copyDeviceToHost(cudaStream_t stream = nullptr)
    {
        copyDeviceToHostAsync(stream);
        host_device_detail::synchronize(stream);
    }

private:
    /** Computes an in-place inclusive prefix sum for one device range array. */
    static void buildPrefixSum(neighborRangeContainer& ranges, cudaStream_t stream);

    /** Selects particle or expanded-node ranges for the active interaction type. */
    neighborRangeContainer& activeNeighborRanges() noexcept { return surfaceNodeNeighborRanges_.hostSize() > 0 ? surfaceNodeNeighborRanges_ : particleNeighborRanges_; }

    /** Enqueues a same-device copy used to preserve compact contact histories. */
    template <class Value> static void copyDeviceArrayAsync(Value* destination, const Value* source, int count, cudaStream_t stream)
    {
        if (count == 0)
        {
            return;
        }
#if FUNDEM_HAS_CUDA
        host_device_detail::checkCuda(cudaMemcpyAsync(destination, source, static_cast<std::size_t>(count) * sizeof(Value), cudaMemcpyDeviceToDevice, stream), "cudaMemcpyAsync contact history");
#else
        (void)destination;
        (void)source;
        (void)stream;
        host_device_detail::throwCUDAUnavailable();
#endif
    }

    /** Reconstructs compact device history from any contacts retained on the host. */
    void initializeContactHistoryAsync(int particleCount, bool usesSurfaceNodeRanges, cudaStream_t stream)
    {
        const int rangeCount = usesSurfaceNodeRanges ? surfaceNodeMappings_.nodeCount() : particleCount;
        const auto& currentContacts = contacts_.host();
        const int contactCount = static_cast<int>(currentContacts.size());
        std::vector<int> rangeIndices(contactCount, -1);
        std::vector<int> counts(rangeCount, 0);

        for (int contactIndex = 0; contactIndex < contactCount; ++contactIndex)
        {
            const contact& currentContact = currentContacts[contactIndex];
            const int masterParticleIndex = currentContact.masterParticleIndex();
            if (masterParticleIndex < 0 || masterParticleIndex >= particleCount)
            {
                continue;
            }

            int rangeIndex = masterParticleIndex;
            if (usesSurfaceNodeRanges)
            {
                rangeIndex = surfaceNodeMappings_.mappingIndex(masterParticleIndex, currentContact.particleSurfaceNodeIndex());
                if (rangeIndex < 0)
                {
                    continue;
                }
            }
            rangeIndices[contactIndex] = rangeIndex;
            ++counts[rangeIndex];
        }

        auto& rangeHistories = neighborRangeHistories_.host();
        rangeHistories.resize(rangeCount);
        std::vector<int> writeIndices(rangeCount, 0);
        int historyCount = 0;
        for (int rangeIndex = 0; rangeIndex < rangeCount; ++rangeIndex)
        {
            writeIndices[rangeIndex] = historyCount;
            historyCount += counts[rangeIndex];
            rangeHistories[rangeIndex] = fundem::neighborRangeHistory(historyCount);
        }

        auto& contactHistories = contactHistories_.host();
        contactHistories.resize(historyCount);
        for (int contactIndex = 0; contactIndex < contactCount; ++contactIndex)
        {
            const int rangeIndex = rangeIndices[contactIndex];
            if (rangeIndex < 0)
            {
                continue;
            }
            const contact& currentContact = currentContacts[contactIndex];
            contactHistories[writeIndices[rangeIndex]++] = fundem::contactHistory(currentContact.slaveParticleIndex(),
                                                                                  currentContact.slidingSpringDeformation(),
                                                                                  currentContact.rollingSpringDeformation(),
                                                                                  currentContact.torsionalSpringDeformation());
        }

        contactHistories_.copyHostToDeviceAsync(stream);
        neighborRangeHistories_.copyHostToDeviceAsync(stream);
        preserveInitializedContactHistory_ = historyCount > 0;
    }

    /** Applies grow-versus-reset policy and uploads every interaction dependency. */
    void initializeStorage(int particleCount, int surfaceNodeCount, bool usesSurfaceNodeRanges, cudaStream_t stream)
    {
        const bool sizeDecreased = particleCount < static_cast<int>(particleNeighborRanges_.hostSize()) || surfaceNodeCount < static_cast<int>(surfaceNodeNeighborRanges_.hostSize());
        if (sizeDecreased)
        {
            resetStorage();
        }

        particleNeighborRanges_.host().resize(particleCount);
        surfaceNodeNeighborRanges_.host().resize(surfaceNodeCount);

        bonds_.copyHostToDeviceAsync(stream);
        particleNeighborRanges_.copyHostToDeviceAsync(stream);
        surfaceNodeNeighborRanges_.copyHostToDeviceAsync(stream);

        initializeContactHistoryAsync(particleCount, usesSurfaceNodeRanges, stream);

        host_device_detail::synchronize(stream);
    }

    /** Drops contacts, histories, candidates, bonds, and range allocations together. */
    void resetStorage() noexcept
    {
        contacts_ = contactContainer{};
        bonds_ = bondContainer{};
        particleNeighborRanges_ = neighborRangeContainer{};
        surfaceNodeNeighborRanges_ = neighborRangeContainer{};
        particleNeighbors_ = particleNeighborContainer{};
        contactHistories_ = contactHistoryContainer{};
        neighborRangeHistories_ = neighborRangeHistoryContainer{};
        preserveInitializedContactHistory_ = false;
    }

    contactContainer contacts_;                            ///< Current contacts and their response state.
    bondContainer bonds_;                                  ///< User-defined persistent bonds.
    surfaceNodeMappingContainer surfaceNodeMappings_;      ///< Expanded LS node-to-owner mapping.
    neighborRangeContainer particleNeighborRanges_;        ///< Per-particle count and prefix sum.
    neighborRangeContainer surfaceNodeNeighborRanges_;     ///< Per-expanded-node count and prefix sum.
    particleNeighborContainer particleNeighbors_;          ///< Broad-phase candidate entries.
    contactHistoryContainer contactHistories_;             ///< Previous-step compact contact history.
    neighborRangeHistoryContainer neighborRangeHistories_; ///< Previous-step range boundaries.
    bool preserveInitializedContactHistory_{false};        ///< Protect restored host history for one step.
};

} // namespace fundem
