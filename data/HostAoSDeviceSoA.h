/**
 * @file HostAoSDeviceSoA.h
 * @brief Implements synchronized host AoS and device SoA container storage.
 */
#pragma once

#include "ContainerObject.h"
#include "CudaTypes.h"
#include "DeviceLayout.h"

#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace fundem
{
namespace host_device_detail
{

#if FUNDEM_HAS_CUDA
/**
 * Throws a descriptive exception when a CUDA runtime operation fails.
 * @param error CUDA runtime status.
 * @param operation Name of the failed operation.
 */
inline void checkCuda(cudaError_t error, const char* operation)
{
    if (error != cudaSuccess)
    {
        throw std::runtime_error(std::string(operation) + " failed: " + cudaGetErrorString(error));
    }
}

/** Waits for all previously submitted work in @p stream. */
inline void synchronize(cudaStream_t stream = nullptr) { checkCuda(cudaStreamSynchronize(stream), "cudaStreamSynchronize"); }
#else
[[noreturn]] inline void throwCUDAUnavailable() { throw std::runtime_error("The CUDA backend is unavailable in this FunDEM build."); }

inline void synchronize(cudaStream_t = nullptr) noexcept {}
#endif

template <class Pointer> struct memberPointerTraits;

template <class Owner, class Value> struct memberPointerTraits<Value Owner::*> {
    using owner_type = Owner;
    using value_type = Value;
};

template <class Field> using fieldOwnerT = typename memberPointerTraits<std::remove_cv_t<decltype(Field::member)>>::owner_type;

template <class Field> using fieldValueT = typename memberPointerTraits<std::remove_cv_t<decltype(Field::member)>>::value_type;

template <class HostObject, class Field> inline constexpr bool fieldBelongsToHost = std::is_same_v<fieldOwnerT<Field>, HostObject> || std::is_base_of_v<fieldOwnerT<Field>, HostObject>;

template <class...> inline constexpr bool dependentFalse = false;

template <class Tag, class... Fields> struct findField;

template <class Tag, class... Rest> struct findField<Tag, Tag, Rest...> {
    using type = Tag;
};

template <class Tag, class First, class... Rest> struct findField<Tag, First, Rest...> : findField<Tag, Rest...> {};

template <class Tag> struct findField<Tag> {
    static_assert(dependentFalse<Tag>, "The field tag is not registered in DeviceLayout.");
};

template <class Tag, class... Fields> using findFieldT = typename findField<Tag, Fields...>::type;

template <class... Fields> struct haveUniqueTags;

template <> struct haveUniqueTags<> : std::true_type {};

template <class First, class... Rest> struct haveUniqueTags<First, Rest...> : std::bool_constant<((!std::is_same_v<First, Rest>) && ...) && haveUniqueTags<Rest...>::value> {};

/** Move-only owner for one raw device array. */
template <class Value> class deviceArray
{
public:
    static_assert(std::is_trivially_copyable_v<Value>, "Device fields must be trivially copyable.");
    static_assert(!std::is_const_v<Value> && !std::is_volatile_v<Value>, "Const and volatile Device fields are not supported.");
    static_assert(!std::is_array_v<Value>, "Raw array fields are not supported.");

    deviceArray() = default;

    ~deviceArray() { reset(); }

    deviceArray(const deviceArray&) = delete;
    deviceArray& operator=(const deviceArray&) = delete;

    deviceArray(deviceArray&& other) noexcept : data_(std::exchange(other.data_, nullptr)) {}

    deviceArray& operator=(deviceArray&& other) noexcept
    {
        if (this != &other)
        {
            reset();
            data_ = std::exchange(other.data_, nullptr);
        }
        return *this;
    }

    /**
     * Replaces the allocation with storage for exactly @p count values.
     * @param count Number of values to allocate; zero releases the allocation.
     * @throws std::length_error if the byte count overflows.
     */
    void allocate(std::size_t count)
    {
        if (count == 0)
        {
            reset();
            return;
        }

        if (count > std::numeric_limits<std::size_t>::max() / sizeof(Value))
        {
            throw std::length_error("CUDA allocation size overflow.");
        }

        Value* replacement = nullptr;
#if FUNDEM_HAS_CUDA
        checkCuda(cudaMalloc(reinterpret_cast<void**>(&replacement), count * sizeof(Value)), "cudaMalloc");
#else
        (void)replacement;
        throwCUDAUnavailable();
#endif

        reset();
        data_ = replacement;
    }

    /** Releases the device allocation, if present. */
    void reset() noexcept
    {
        if (data_ != nullptr)
        {
#if FUNDEM_HAS_CUDA
            // Destructors cannot report cudaFree failures.
            cudaFree(data_);
#endif
            data_ = nullptr;
        }
    }

    Value* data() noexcept { return data_; }

    const Value* data() const noexcept { return data_; }

private:
    Value* data_{nullptr}; ///< Owned device allocation, or null when empty.
};

/** Associates one registered field tag with its owned device allocation. */
template <class Field> struct deviceSlot {
    deviceArray<fieldValueT<Field>> array; ///< Device values registered by @c Field.
};

} // namespace host_device_detail

/**
 * Owns complete Host objects as AoS and registered Device fields as SoA.
 *
 * HostObject supplies a DeviceLayout. Full Host-to-Device copies update the
 * logical Device size and grow capacity when required. Full Device-to-Host
 * copies resize the Host vector to the logical Device size.
 *
 * Reallocation intentionally discards old Device values because it is only
 * performed immediately before a full Host-to-Device overwrite.
 */
template <class HostObject, class Layout = typename HostObject::DeviceLayout> class hostAoSDeviceSoA;

/** Concrete AoS/SoA storage generated for a registered field list. */
template <class HostObject, class... Fields> class hostAoSDeviceSoA<HostObject, deviceLayout<Fields...>>
{
    static_assert(sizeof...(Fields) > 0, "DeviceLayout must contain at least one field.");
    static_assert(std::is_default_constructible_v<HostObject>, "HostObject must be default constructible.");
    static_assert((std::is_member_object_pointer_v<std::remove_cv_t<decltype(Fields::member)>> && ...), "Every DeviceLayout field must provide a non-static data-member pointer.");
    static_assert((host_device_detail::fieldBelongsToHost<HostObject, Fields> && ...), "A registered field does not belong to HostObject or a base class.");
    static_assert((std::is_trivially_copyable_v<host_device_detail::fieldValueT<Fields>> && ...), "Device fields must be trivially copyable.");
    static_assert(host_device_detail::haveUniqueTags<Fields...>::value, "Device field tags must be unique.");

    template <class Field> using Slot = host_device_detail::deviceSlot<Field>;

    using DeviceTuple = std::tuple<Slot<Fields>...>;

public:
    using host_type = HostObject;
    using host_container_type = std::vector<HostObject>;
    using layout_type = deviceLayout<Fields...>;

    hostAoSDeviceSoA() = default;
    ~hostAoSDeviceSoA() = default;

    hostAoSDeviceSoA(const hostAoSDeviceSoA&) = delete;
    hostAoSDeviceSoA& operator=(const hostAoSDeviceSoA&) = delete;

    hostAoSDeviceSoA(hostAoSDeviceSoA&& other) noexcept
        : host_(std::move(other.host_)), device_(std::move(other.device_)), deviceSize_(std::exchange(other.deviceSize_, 0)), deviceCapacity_(std::exchange(other.deviceCapacity_, 0)),
          useDevice_(std::exchange(other.useDevice_, false))
    {}

    hostAoSDeviceSoA& operator=(hostAoSDeviceSoA&& other) noexcept
    {
        if (this != &other)
        {
            host_ = std::move(other.host_);
            device_ = std::move(other.device_);
            deviceSize_ = std::exchange(other.deviceSize_, 0);
            deviceCapacity_ = std::exchange(other.deviceCapacity_, 0);
            useDevice_ = std::exchange(other.useDevice_, false);
        }
        return *this;
    }

    bool usesDevice() const noexcept { return useDevice_; }

    void setUseDevice(bool value) noexcept { useDevice_ = value; }

    host_container_type& host() noexcept
    {
        bindHostObjects();
        return host_;
    }

    const host_container_type& host() const noexcept
    {
        bindHostObjects();
        return host_;
    }

    std::size_t hostSize() const noexcept { return host_.size(); }

    std::size_t deviceSize() const noexcept { return deviceSize_; }

    std::size_t deviceCapacity() const noexcept { return deviceCapacity_; }

    /** Returns the capacity of all registered device fields in GiB. */
    double deviceMemoryGB() const noexcept
    {
        constexpr std::size_t bytesPerElement = (sizeof(host_device_detail::fieldValueT<Fields>) + ...);
        return static_cast<double>(deviceCapacity_) * static_cast<double>(bytesPerElement) / (1024.0 * 1024.0 * 1024.0);
    }

    /**
     * Changes the logical device length without reallocating.
     * @param size New logical length.
     * @throws std::out_of_range if @p size exceeds device capacity.
     */
    void setDeviceSize(std::size_t size)
    {
        if (size > deviceCapacity_)
        {
            throw std::out_of_range("setDeviceSize(): size exceeds Device capacity.");
        }
        deviceSize_ = size;
    }

    /**
     * Ensures capacity for @p size objects and sets the logical device length.
     * Existing device values may be discarded when capacity grows.
     */
    void resizeDevice(std::size_t size)
    {
        ensureDeviceCapacity(size);
        deviceSize_ = size;
    }

    /** Releases every device field and resets logical size and capacity to zero. */
    void resetDevice() noexcept
    {
        device_ = DeviceTuple{};
        deviceSize_ = 0;
        deviceCapacity_ = 0;
    }

    template <class Tag> using field_type = host_device_detail::findFieldT<Tag, Fields...>;

    template <class Tag> using device_value_type = host_device_detail::fieldValueT<field_type<Tag>>;

    template <class Tag> device_value_type<Tag>* device() noexcept { return fieldArray<field_type<Tag>>().data(); }

    template <class Tag> const device_value_type<Tag>* device() const noexcept { return fieldArray<field_type<Tag>>().data(); }

    /**
     * Enqueues a complete host-AoS to device-SoA copy.
     * Device capacity grows as needed and logical device size becomes `hostSize()`.
     * @param stream CUDA stream receiving the copy operations.
     */
    void copyHostToDeviceAsync(cudaStream_t stream = nullptr)
    {
        const std::size_t count = host_.size();
        ensureDeviceCapacity(count);

        // Do not advertise partially enqueued data as valid if a copy throws.
        deviceSize_ = 0;
        (copyOneHostToDeviceAsync<Fields>(count, stream), ...);
        deviceSize_ = count;
    }

    /**
     * Performs a complete host-to-device copy and waits for @p stream.
     * @param stream CUDA stream receiving the copy operations.
     */
    void copyHostToDevice(cudaStream_t stream = nullptr)
    {
        copyHostToDeviceAsync(stream);
        synchronize(stream);
    }

    /**
     * Enqueues a complete device-SoA to host-AoS copy.
     * The host vector is resized to the logical device size before the copy.
     * @param stream CUDA stream receiving the copy operations.
     */
    void copyDeviceToHostAsync(cudaStream_t stream = nullptr)
    {
        host_.resize(deviceSize_);
        (copyOneDeviceToHostAsync<Fields>(deviceSize_, stream), ...);
    }

    /**
     * Performs a complete device-to-host copy and waits for @p stream.
     * @param stream CUDA stream receiving the copy operations.
     */
    void copyDeviceToHost(cudaStream_t stream = nullptr)
    {
        copyDeviceToHostAsync(stream);
        synchronize(stream);
    }

    /**
     * Enqueues one registered host field for copying to the corresponding device array.
     * @tparam Tag Registered field tag.
     * @throws std::logic_error if host and device logical sizes differ.
     */
    template <class Tag> void copyHostFieldToDeviceAsync(cudaStream_t stream = nullptr)
    {
        requireSameSize("copyHostFieldToDeviceAsync()");
        copyOneHostToDeviceAsync<field_type<Tag>>(host_.size(), stream);
    }

    /** Copies one registered host field to the device and waits for @p stream. */
    template <class Tag> void copyHostFieldToDevice(cudaStream_t stream = nullptr)
    {
        copyHostFieldToDeviceAsync<Tag>(stream);
        synchronize(stream);
    }

    /**
     * Enqueues one registered device field for copying to the host objects.
     * @tparam Tag Registered field tag.
     * @throws std::logic_error if host and device logical sizes differ.
     */
    template <class Tag> void copyDeviceFieldToHostAsync(cudaStream_t stream = nullptr)
    {
        requireSameSize("copyDeviceFieldToHostAsync()");
        copyOneDeviceToHostAsync<field_type<Tag>>(deviceSize_, stream);
    }

    /** Copies one registered device field to the host and waits for @p stream. */
    template <class Tag> void copyDeviceFieldToHost(cudaStream_t stream = nullptr)
    {
        copyDeviceFieldToHostAsync<Tag>(stream);
        synchronize(stream);
    }

    /** Waits for all previously submitted work in @p stream. */
    void synchronize(cudaStream_t stream = nullptr) const { host_device_detail::synchronize(stream); }

private:
    /** Refreshes insertion metadata after host-vector growth or relocation. */
    void bindHostObjects() const noexcept
    {
        if constexpr (std::is_base_of_v<containerObject, HostObject>)
        {
            if (host_.empty() || (host_.front().isBoundToContainerIndex(0) && host_.back().isBoundToContainerIndex(static_cast<int>(host_.size()) - 1)))
            {
                return;
            }
            for (std::size_t index = 0; index < host_.size(); ++index)
            {
                host_[index].bindToContainer(static_cast<int>(index));
            }
        }
    }

    /** Selects the device allocation associated with one field tag. */
    template <class Field> auto& fieldArray() noexcept { return std::get<Slot<Field>>(device_).array; }

    /** Const overload of `fieldArray`. */
    template <class Field> const auto& fieldArray() const noexcept { return std::get<Slot<Field>>(device_).array; }

    /** Reallocates every device field together when @p required exceeds capacity. */
    void ensureDeviceCapacity(std::size_t required)
    {
        if (required <= deviceCapacity_)
        {
            return;
        }

        DeviceTuple replacement;
        (std::get<Slot<Fields>>(replacement).array.allocate(required), ...);

        device_ = std::move(replacement);
        deviceCapacity_ = required;
        deviceSize_ = 0;
    }

    /** Enforces equal logical lengths before a single-field transfer. */
    void requireSameSize(const char* operation) const
    {
        if (host_.size() != deviceSize_)
        {
            throw std::logic_error(std::string(operation) + " requires hostSize() == deviceSize().");
        }
    }

    /** Enqueues the strided AoS-to-contiguous-SoA transfer for one field. */
    template <class Field> void copyOneHostToDeviceAsync(std::size_t count, cudaStream_t stream)
    {
        if (count == 0)
        {
            return;
        }

        using Value = host_device_detail::fieldValueT<Field>;
        const Value* source = std::addressof(host_.front().*(Field::member));

#if FUNDEM_HAS_CUDA
        host_device_detail::checkCuda(cudaMemcpy2DAsync(fieldArray<Field>().data(), sizeof(Value), source, sizeof(HostObject), sizeof(Value), count, cudaMemcpyHostToDevice, stream),
                                      "cudaMemcpy2DAsync Host AoS -> Device SoA");
#else
        (void)source;
        (void)stream;
        host_device_detail::throwCUDAUnavailable();
#endif
    }

    /** Enqueues the contiguous-SoA-to-strided-AoS transfer for one field. */
    template <class Field> void copyOneDeviceToHostAsync(std::size_t count, cudaStream_t stream)
    {
        if (count == 0)
        {
            return;
        }

        using Value = host_device_detail::fieldValueT<Field>;
        Value* destination = std::addressof(host_.front().*(Field::member));

#if FUNDEM_HAS_CUDA
        host_device_detail::checkCuda(cudaMemcpy2DAsync(destination, sizeof(HostObject), fieldArray<Field>().data(), sizeof(Value), sizeof(Value), count, cudaMemcpyDeviceToHost, stream),
                                      "cudaMemcpy2DAsync Device SoA -> Host AoS");
#else
        (void)destination;
        (void)stream;
        host_device_detail::throwCUDAUnavailable();
#endif
    }

    host_container_type host_;      ///< Complete host objects stored as an array of structures.
    DeviceTuple device_;            ///< Registered device fields stored as a structure of arrays.
    std::size_t deviceSize_{0};     ///< Logical number of valid device objects.
    std::size_t deviceCapacity_{0}; ///< Allocated element capacity of every device field.
    bool useDevice_{false};         ///< Whether the owning solver selected the device representation.
};

} // namespace fundem
