/**
 * @file DeviceLayout.h
 * @brief Declares compile-time field layouts for device Structure-of-Arrays storage.
 */
#pragma once

namespace fundem
{

/**
 * Lists field descriptors stored as independent Device arrays.
 *
 * Each descriptor is also its public field tag and provides a static constexpr
 * member pointer. Keeping that pointer out of the descriptor's template type
 * preserves private Host-object members with CUDA 13 and later.
 */
template <class... Fields> struct deviceLayout {};

/** Concatenates any number of device layouts at compile time. */
template <class... Layouts> struct concatDeviceLayouts;

template <> struct concatDeviceLayouts<> {
    using type = deviceLayout<>;
};

template <class... Fields> struct concatDeviceLayouts<deviceLayout<Fields...>> {
    using type = deviceLayout<Fields...>;
};

template <class... LeftFields, class... RightFields, class... Rest>
struct concatDeviceLayouts<deviceLayout<LeftFields...>, deviceLayout<RightFields...>, Rest...> : concatDeviceLayouts<deviceLayout<LeftFields..., RightFields...>, Rest...> {};
/** Resulting layout type produced by `concatDeviceLayouts`. */

template <class... Layouts> using concatDeviceLayoutsT = typename concatDeviceLayouts<Layouts...>::type;

template <class HostObject, class Layout> class hostAoSDeviceSoA;

} // namespace fundem
