/**
 * @file Constants.h
 * @brief Defines scalar types, constants, tolerances, and host/device qualifiers.
 */
#pragma once

#include <cmath>

#ifndef FUNDEM_MATH_HD
#if defined(__CUDACC__)
#define FUNDEM_MATH_HD __host__ __device__
#else
#define FUNDEM_MATH_HD
#endif
#endif

namespace fundem::math
{

using Real = double;

inline constexpr Real pi = 3.141592653589793238462643383279502884;
inline constexpr Real twoPi = 2.0 * pi;
inline constexpr Real halfPi = 0.5 * pi;
inline constexpr Real defaultTolerance = 1.0e-12;
inline constexpr Real defaultToleranceSquared = defaultTolerance * defaultTolerance;

namespace detail
{

FUNDEM_MATH_HD inline Real abs(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::fabs(value);
#else
    return std::fabs(value);
#endif
}

FUNDEM_MATH_HD inline Real sqrt(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::sqrt(value);
#else
    return std::sqrt(value);
#endif
}

FUNDEM_MATH_HD inline Real sin(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::sin(value);
#else
    return std::sin(value);
#endif
}

FUNDEM_MATH_HD inline Real cos(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::cos(value);
#else
    return std::cos(value);
#endif
}

FUNDEM_MATH_HD inline Real log(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::log(value);
#else
    return std::log(value);
#endif
}

FUNDEM_MATH_HD inline Real floor(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::floor(value);
#else
    return std::floor(value);
#endif
}

FUNDEM_MATH_HD inline Real pow(Real base, Real exponent) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::pow(base, exponent);
#else
    return std::pow(base, exponent);
#endif
}

FUNDEM_MATH_HD inline bool finite(Real value) noexcept
{
#if defined(__CUDA_ARCH__)
    return ::isfinite(value) != 0;
#else
    return std::isfinite(value);
#endif
}

FUNDEM_MATH_HD constexpr Real max(Real lhs, Real rhs) noexcept { return lhs > rhs ? lhs : rhs; }

FUNDEM_MATH_HD inline Real maxAbs(Real x, Real y, Real z) noexcept { return max(abs(x), max(abs(y), abs(z))); }

FUNDEM_MATH_HD inline Real maxAbs(Real w, Real x, Real y, Real z) noexcept { return max(max(abs(w), abs(x)), max(abs(y), abs(z))); }

} // namespace detail

FUNDEM_MATH_HD constexpr Real degreesToRadians(Real degrees) noexcept { return degrees * (pi / 180.0); }

FUNDEM_MATH_HD constexpr Real radiansToDegrees(Real radians) noexcept { return radians * (180.0 / pi); }

FUNDEM_MATH_HD inline bool isFinite(Real value) noexcept { return detail::finite(value); }

FUNDEM_MATH_HD inline bool nearlyEqual(Real lhs, Real rhs, Real absoluteTolerance = defaultTolerance, Real relativeTolerance = defaultTolerance) noexcept
{
    if (lhs == rhs)
    {
        return true;
    }

    if (!isFinite(lhs) || !isFinite(rhs))
    {
        return false;
    }

    const Real absolute = detail::abs(absoluteTolerance);
    const Real relative = detail::abs(relativeTolerance);
    if (!isFinite(absolute) || !isFinite(relative))
    {
        return false;
    }

    const Real scale = detail::max(1.0, detail::max(detail::abs(lhs), detail::abs(rhs)));
    return detail::abs(lhs - rhs) <= absolute + relative * scale;
}

FUNDEM_MATH_HD inline bool nearlyZero(Real value, Real tolerance = defaultTolerance) noexcept
{
    const Real positiveTolerance = detail::abs(tolerance);
    return isFinite(value) && isFinite(positiveTolerance) && detail::abs(value) <= positiveTolerance;
}

} // namespace fundem::math
