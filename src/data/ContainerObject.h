/**
 * @file ContainerObject.h
 * @brief Defines ownership metadata for objects inserted into solver containers.
 */
#pragma once

#include <atomic>

namespace fundem
{

template <class HostObject, class Layout> class hostAoSDeviceSoA;

/** Tracks the stable index assigned when a user-created object enters a container. */
class containerObject
{
public:
    containerObject() = default;
    /** Copies user state without copying membership in a container. */
    containerObject(const containerObject&) noexcept {}
    /** Assigns user state without changing either object's container membership. */
    containerObject& operator=(const containerObject&) noexcept { return *this; }
    /** Moves user state without transferring membership in a container. */
    containerObject(containerObject&&) noexcept {}
    /** Move-assigns user state without changing either object's container membership. */
    containerObject& operator=(containerObject&&) noexcept { return *this; }

    bool isInContainer() const noexcept { return containerIndex_.load(std::memory_order_relaxed) >= 0; }
    int containerIndex() const noexcept { return containerIndex_.load(std::memory_order_relaxed); }

private:
    template <class HostObject, class Layout> friend class hostAoSDeviceSoA;

    /** Binds the object to one stable host-container index. */
    void bindToContainer(int index) const noexcept { containerIndex_.store(index, std::memory_order_relaxed); }
    /** Tests the exact stable index expected by the owning container. */
    bool isBoundToContainerIndex(int index) const noexcept { return containerIndex_.load(std::memory_order_relaxed) == index; }

    mutable std::atomic<int> containerIndex_{-1}; ///< Bound container index, or `-1` before insertion.
};

} // namespace fundem
