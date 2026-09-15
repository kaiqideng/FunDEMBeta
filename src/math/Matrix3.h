/**
 * @file Matrix3.h
 * @brief Defines the row-major 3-by-3 matrix type and linear algebra operations.
 */
#pragma once

#include "Vector3.h"

#include <cassert>
#include <cstddef>
#include <type_traits>

namespace fundem::math
{

/** Row-major, trivially copyable 3-by-3 matrix shared by host and device code. */
struct Mat3 {
    Real values[9]{}; ///< Row-major matrix coefficients.

    FUNDEM_MATH_HD constexpr Mat3() noexcept = default;

    FUNDEM_MATH_HD constexpr Mat3(Real m00, Real m01, Real m02, Real m10, Real m11, Real m12, Real m20, Real m21, Real m22) noexcept : values{m00, m01, m02, m10, m11, m12, m20, m21, m22} {}

    FUNDEM_MATH_HD static constexpr Mat3 zero() noexcept { return {}; }

    FUNDEM_MATH_HD static constexpr Mat3 identity() noexcept { return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; }

    FUNDEM_MATH_HD static constexpr Mat3 diagonal(const Vec3& value) noexcept { return {value.x, 0.0, 0.0, 0.0, value.y, 0.0, 0.0, 0.0, value.z}; }

    FUNDEM_MATH_HD static constexpr std::size_t size() noexcept { return 9; }

    FUNDEM_MATH_HD constexpr Real* data() noexcept { return values; }

    FUNDEM_MATH_HD constexpr const Real* data() const noexcept { return values; }

    FUNDEM_MATH_HD constexpr Real& operator[](std::size_t index) noexcept
    {
        assert(index < size());
        return values[index];
    }

    FUNDEM_MATH_HD constexpr const Real& operator[](std::size_t index) const noexcept
    {
        assert(index < size());
        return values[index];
    }

    FUNDEM_MATH_HD constexpr Real& operator()(std::size_t row, std::size_t column) noexcept
    {
        assert(row < 3 && column < 3);
        return values[3 * row + column];
    }

    FUNDEM_MATH_HD constexpr const Real& operator()(std::size_t row, std::size_t column) const noexcept
    {
        assert(row < 3 && column < 3);
        return values[3 * row + column];
    }

    FUNDEM_MATH_HD constexpr Vec3 row(std::size_t index) const noexcept
    {
        assert(index < 3);
        return {(*this)(index, 0), (*this)(index, 1), (*this)(index, 2)};
    }

    FUNDEM_MATH_HD constexpr Vec3 column(std::size_t index) const noexcept
    {
        assert(index < 3);
        return {(*this)(0, index), (*this)(1, index), (*this)(2, index)};
    }

    FUNDEM_MATH_HD constexpr Mat3& operator+=(const Mat3& rhs) noexcept
    {
        for (std::size_t index = 0; index < size(); ++index)
        {
            values[index] += rhs.values[index];
        }
        return *this;
    }

    FUNDEM_MATH_HD constexpr Mat3& operator-=(const Mat3& rhs) noexcept
    {
        for (std::size_t index = 0; index < size(); ++index)
        {
            values[index] -= rhs.values[index];
        }
        return *this;
    }

    FUNDEM_MATH_HD constexpr Mat3& operator*=(Real scalar) noexcept
    {
        for (std::size_t index = 0; index < size(); ++index)
        {
            values[index] *= scalar;
        }
        return *this;
    }

    FUNDEM_MATH_HD constexpr Mat3& operator/=(Real scalar) noexcept
    {
        for (std::size_t index = 0; index < size(); ++index)
        {
            values[index] /= scalar;
        }
        return *this;
    }
};

static_assert(std::is_standard_layout_v<Mat3>);
static_assert(std::is_trivially_copyable_v<Mat3>);
static_assert(sizeof(Mat3) == 9 * sizeof(Real));

FUNDEM_MATH_HD constexpr bool operator==(const Mat3& lhs, const Mat3& rhs) noexcept
{
    for (std::size_t index = 0; index < Mat3::size(); ++index)
    {
        if (lhs[index] != rhs[index])
        {
            return false;
        }
    }
    return true;
}

FUNDEM_MATH_HD constexpr bool operator!=(const Mat3& lhs, const Mat3& rhs) noexcept { return !(lhs == rhs); }

FUNDEM_MATH_HD constexpr Mat3 operator+(Mat3 lhs, const Mat3& rhs) noexcept { return lhs += rhs; }

FUNDEM_MATH_HD constexpr Mat3 operator-(Mat3 lhs, const Mat3& rhs) noexcept { return lhs -= rhs; }

FUNDEM_MATH_HD constexpr Mat3 operator-(Mat3 value) noexcept { return value *= -1.0; }

FUNDEM_MATH_HD constexpr Mat3 operator*(Mat3 matrix, Real scalar) noexcept { return matrix *= scalar; }

FUNDEM_MATH_HD constexpr Mat3 operator*(Real scalar, Mat3 matrix) noexcept { return matrix *= scalar; }

FUNDEM_MATH_HD constexpr Mat3 operator/(Mat3 matrix, Real scalar) noexcept { return matrix /= scalar; }

FUNDEM_MATH_HD constexpr Vec3 operator*(const Mat3& matrix, const Vec3& vector) noexcept
{
    return {matrix(0, 0) * vector.x + matrix(0, 1) * vector.y + matrix(0, 2) * vector.z,
            matrix(1, 0) * vector.x + matrix(1, 1) * vector.y + matrix(1, 2) * vector.z,
            matrix(2, 0) * vector.x + matrix(2, 1) * vector.y + matrix(2, 2) * vector.z};
}

FUNDEM_MATH_HD constexpr Mat3 operator*(const Mat3& lhs, const Mat3& rhs) noexcept
{
    Mat3 result;
    for (std::size_t row = 0; row < 3; ++row)
    {
        for (std::size_t column = 0; column < 3; ++column)
        {
            result(row, column) = lhs(row, 0) * rhs(0, column) + lhs(row, 1) * rhs(1, column) + lhs(row, 2) * rhs(2, column);
        }
    }
    return result;
}

FUNDEM_MATH_HD constexpr Real trace(const Mat3& matrix) noexcept { return matrix(0, 0) + matrix(1, 1) + matrix(2, 2); }

FUNDEM_MATH_HD constexpr Mat3 transposed(const Mat3& matrix) noexcept
{
    return {matrix(0, 0), matrix(1, 0), matrix(2, 0), matrix(0, 1), matrix(1, 1), matrix(2, 1), matrix(0, 2), matrix(1, 2), matrix(2, 2)};
}

FUNDEM_MATH_HD constexpr Real determinant(const Mat3& matrix) noexcept
{
    return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)) - matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0)) +
           matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
}

namespace detail
{

FUNDEM_MATH_HD constexpr Mat3 adjugate(const Mat3& matrix) noexcept
{
    return {matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1),
            matrix(0, 2) * matrix(2, 1) - matrix(0, 1) * matrix(2, 2),
            matrix(0, 1) * matrix(1, 2) - matrix(0, 2) * matrix(1, 1),
            matrix(1, 2) * matrix(2, 0) - matrix(1, 0) * matrix(2, 2),
            matrix(0, 0) * matrix(2, 2) - matrix(0, 2) * matrix(2, 0),
            matrix(0, 2) * matrix(1, 0) - matrix(0, 0) * matrix(1, 2),
            matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0),
            matrix(0, 1) * matrix(2, 0) - matrix(0, 0) * matrix(2, 1),
            matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0)};
}

} // namespace detail

FUNDEM_MATH_HD inline bool isFinite(const Mat3& matrix) noexcept
{
    for (std::size_t index = 0; index < Mat3::size(); ++index)
    {
        if (!isFinite(matrix[index]))
        {
            return false;
        }
    }
    return true;
}

/**
 * Computes a scale-aware matrix inverse.
 * @param matrix Matrix to invert.
 * @param result Receives the inverse and remains unchanged on failure.
 * @param relativeTolerance Reciprocal-condition threshold.
 * @return False for non-finite or numerically singular matrices.
 */
FUNDEM_MATH_HD inline bool tryInverse(const Mat3& matrix, Mat3& result, Real relativeTolerance = defaultTolerance) noexcept
{
    Real scale = 0.0;
    for (std::size_t index = 0; index < Mat3::size(); ++index)
    {
        if (!isFinite(matrix[index]))
        {
            return false;
        }
        scale = detail::max(scale, detail::abs(matrix[index]));
    }

    const Real tolerance = detail::abs(relativeTolerance);
    if (scale == 0.0 || !isFinite(tolerance))
    {
        return false;
    }

    const Mat3 scaled = matrix / scale;
    const Real scaledDeterminant = determinant(scaled);
    if (!isFinite(scaledDeterminant))
    {
        return false;
    }

    const Mat3 adjugate = detail::adjugate(scaled);
    Real adjugateScale = 0.0;
    for (std::size_t index = 0; index < Mat3::size(); ++index)
    {
        adjugateScale = detail::max(adjugateScale, detail::abs(adjugate[index]));
    }

    // abs(det) / (max(abs(A)) * max(abs(adj(A)))) estimates the
    // reciprocal condition number. max(abs(scaled)) is one.
    if (adjugateScale == 0.0 || detail::abs(scaledDeterminant) <= tolerance * adjugateScale)
    {
        return false;
    }

    Mat3 candidate = adjugate / (scaledDeterminant * scale);
    if (!isFinite(candidate))
    {
        return false;
    }

    result = candidate;
    return true;
}

FUNDEM_MATH_HD constexpr Mat3 outerProduct(const Vec3& lhs, const Vec3& rhs) noexcept
{
    return {lhs.x * rhs.x, lhs.x * rhs.y, lhs.x * rhs.z, lhs.y * rhs.x, lhs.y * rhs.y, lhs.y * rhs.z, lhs.z * rhs.x, lhs.z * rhs.y, lhs.z * rhs.z};
}

FUNDEM_MATH_HD constexpr Mat3 crossProductMatrix(const Vec3& value) noexcept { return {0.0, -value.z, value.y, value.z, 0.0, -value.x, -value.y, value.x, 0.0}; }

FUNDEM_MATH_HD inline bool nearlyEqual(const Mat3& lhs, const Mat3& rhs, Real absoluteTolerance = defaultTolerance, Real relativeTolerance = defaultTolerance) noexcept
{
    for (std::size_t index = 0; index < Mat3::size(); ++index)
    {
        if (!nearlyEqual(lhs[index], rhs[index], absoluteTolerance, relativeTolerance))
        {
            return false;
        }
    }
    return true;
}

} // namespace fundem::math
