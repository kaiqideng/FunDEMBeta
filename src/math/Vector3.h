/**
 * @file Vector3.h
 * @brief Defines three-dimensional vectors and host/device vector algebra.
 */
#pragma once

#include "Constants.h"

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace fundem::math
{

/** Trivially copyable three-dimensional vector shared by host and device code. */
struct Vec3 {
    Real x{0.0}; ///< x component.
    Real y{0.0}; ///< y component.
    Real z{0.0}; ///< z component.

    FUNDEM_MATH_HD constexpr Vec3() noexcept = default;

    FUNDEM_MATH_HD constexpr Vec3(Real xValue, Real yValue, Real zValue) noexcept : x(xValue), y(yValue), z(zValue) {}

    FUNDEM_MATH_HD explicit constexpr Vec3(Real value) noexcept : x(value), y(value), z(value) {}

    FUNDEM_MATH_HD static constexpr Vec3 zero() noexcept { return {}; }

    FUNDEM_MATH_HD static constexpr Vec3 one() noexcept { return {1.0, 1.0, 1.0}; }

    FUNDEM_MATH_HD static constexpr Vec3 unitX() noexcept { return {1.0, 0.0, 0.0}; }

    FUNDEM_MATH_HD static constexpr Vec3 unitY() noexcept { return {0.0, 1.0, 0.0}; }

    FUNDEM_MATH_HD static constexpr Vec3 unitZ() noexcept { return {0.0, 0.0, 1.0}; }

    FUNDEM_MATH_HD static constexpr std::size_t size() noexcept { return 3; }

    FUNDEM_MATH_HD constexpr Real& operator[](std::size_t index) noexcept
    {
        assert(index < size());
        return index == 0 ? x : (index == 1 ? y : z);
    }

    FUNDEM_MATH_HD constexpr const Real& operator[](std::size_t index) const noexcept
    {
        assert(index < size());
        return index == 0 ? x : (index == 1 ? y : z);
    }

    FUNDEM_MATH_HD constexpr Vec3& operator+=(const Vec3& rhs) noexcept
    {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    FUNDEM_MATH_HD constexpr Vec3& operator-=(const Vec3& rhs) noexcept
    {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    FUNDEM_MATH_HD constexpr Vec3& operator*=(Real scalar) noexcept
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    FUNDEM_MATH_HD constexpr Vec3& operator/=(Real scalar) noexcept
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }
};

static_assert(std::is_standard_layout_v<Vec3>);
static_assert(std::is_trivially_copyable_v<Vec3>);
static_assert(sizeof(Vec3) == 3 * sizeof(Real));

FUNDEM_MATH_HD constexpr bool operator==(const Vec3& lhs, const Vec3& rhs) noexcept { return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z; }

FUNDEM_MATH_HD constexpr bool operator!=(const Vec3& lhs, const Vec3& rhs) noexcept { return !(lhs == rhs); }

FUNDEM_MATH_HD constexpr Vec3 operator+(Vec3 lhs, const Vec3& rhs) noexcept { return lhs += rhs; }

FUNDEM_MATH_HD constexpr Vec3 operator-(Vec3 lhs, const Vec3& rhs) noexcept { return lhs -= rhs; }

FUNDEM_MATH_HD constexpr Vec3 operator-(const Vec3& value) noexcept { return {-value.x, -value.y, -value.z}; }

FUNDEM_MATH_HD constexpr Vec3 operator*(Vec3 value, Real scalar) noexcept { return value *= scalar; }

FUNDEM_MATH_HD constexpr Vec3 operator*(Real scalar, Vec3 value) noexcept { return value *= scalar; }

FUNDEM_MATH_HD constexpr Vec3 operator/(Vec3 value, Real scalar) noexcept { return value /= scalar; }

FUNDEM_MATH_HD constexpr Real dot(const Vec3& lhs, const Vec3& rhs) noexcept { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }

FUNDEM_MATH_HD constexpr Vec3 cross(const Vec3& lhs, const Vec3& rhs) noexcept { return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x}; }

FUNDEM_MATH_HD constexpr Vec3 hadamardProduct(const Vec3& lhs, const Vec3& rhs) noexcept { return {lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z}; }

FUNDEM_MATH_HD constexpr Vec3 componentMin(const Vec3& lhs, const Vec3& rhs) noexcept { return {lhs.x < rhs.x ? lhs.x : rhs.x, lhs.y < rhs.y ? lhs.y : rhs.y, lhs.z < rhs.z ? lhs.z : rhs.z}; }

FUNDEM_MATH_HD constexpr Vec3 componentMax(const Vec3& lhs, const Vec3& rhs) noexcept { return {lhs.x > rhs.x ? lhs.x : rhs.x, lhs.y > rhs.y ? lhs.y : rhs.y, lhs.z > rhs.z ? lhs.z : rhs.z}; }

FUNDEM_MATH_HD constexpr Real normSquared(const Vec3& value) noexcept { return dot(value, value); }

FUNDEM_MATH_HD inline Real norm(const Vec3& value) noexcept
{
    if (!isFinite(value.x) || !isFinite(value.y) || !isFinite(value.z))
    {
        return detail::sqrt(normSquared(value));
    }

    const Real scale = detail::maxAbs(value.x, value.y, value.z);
    if (scale == 0.0)
    {
        return 0.0;
    }

    const Vec3 scaled = value / scale;
    return scale * detail::sqrt(normSquared(scaled));
}

FUNDEM_MATH_HD constexpr Real distanceSquared(const Vec3& lhs, const Vec3& rhs) noexcept { return normSquared(lhs - rhs); }

FUNDEM_MATH_HD inline Real distance(const Vec3& lhs, const Vec3& rhs) noexcept { return norm(lhs - rhs); }

FUNDEM_MATH_HD inline bool isFinite(const Vec3& value) noexcept { return isFinite(value.x) && isFinite(value.y) && isFinite(value.z); }

FUNDEM_MATH_HD inline bool nearlyEqual(const Vec3& lhs, const Vec3& rhs, Real absoluteTolerance = defaultTolerance, Real relativeTolerance = defaultTolerance) noexcept
{
    return nearlyEqual(lhs.x, rhs.x, absoluteTolerance, relativeTolerance) && nearlyEqual(lhs.y, rhs.y, absoluteTolerance, relativeTolerance) &&
           nearlyEqual(lhs.z, rhs.z, absoluteTolerance, relativeTolerance);
}

/**
 * Normalizes @p value in place using scale-safe arithmetic.
 * @param value Vector to normalize; unchanged on failure.
 * @param tolerance Minimum accepted norm.
 * @return True when a finite unit vector was produced.
 */
FUNDEM_MATH_HD inline bool tryNormalize(Vec3& value, Real tolerance = defaultTolerance) noexcept
{
    if (!isFinite(value))
    {
        return false;
    }

    const Real scale = detail::maxAbs(value.x, value.y, value.z);
    const Real positiveTolerance = detail::abs(tolerance);
    if (scale == 0.0 || !isFinite(positiveTolerance))
    {
        return false;
    }

    const Vec3 scaled = value / scale;
    const Real scaledNorm = detail::sqrt(normSquared(scaled));
    if (!isFinite(scaledNorm) || scaledNorm == 0.0 || scale <= positiveTolerance / scaledNorm)
    {
        return false;
    }

    value = scaled / scaledNorm;
    return true;
}

/** Returns a normalized copy, or the zero vector when normalization is unsafe. */
FUNDEM_MATH_HD inline Vec3 normalizedOrZero(const Vec3& value, Real tolerance = defaultTolerance) noexcept
{
    Vec3 result = value;
    return tryNormalize(result, tolerance) ? result : Vec3::zero();
}

FUNDEM_MATH_HD constexpr Vec3 lerp(const Vec3& from, const Vec3& to, Real fraction) noexcept { return from + (to - from) * fraction; }

} // namespace fundem::math
