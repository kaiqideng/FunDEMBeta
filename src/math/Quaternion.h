/**
 * @file Quaternion.h
 * @brief Defines quaternions, unit rotations, and orientation utilities.
 */
#pragma once

#include "Matrix3.h"

#include <type_traits>

namespace fundem::math
{

/** Hamilton quaternion used for finite, unit-norm rigid-body orientations. */
struct Quaternion {
    Real w{1.0}; ///< Scalar component.
    Real x{0.0}; ///< x component of the vector part.
    Real y{0.0}; ///< y component of the vector part.
    Real z{0.0}; ///< z component of the vector part.

    FUNDEM_MATH_HD constexpr Quaternion() noexcept = default;

    FUNDEM_MATH_HD constexpr Quaternion(Real wValue, Real xValue, Real yValue, Real zValue) noexcept : w(wValue), x(xValue), y(yValue), z(zValue) {}

    FUNDEM_MATH_HD static constexpr Quaternion identity() noexcept { return {}; }

    FUNDEM_MATH_HD static constexpr Quaternion zero() noexcept { return {0.0, 0.0, 0.0, 0.0}; }

    /**
     * Fast path: unitAxis must have unit length and both inputs must be finite.
     */
    FUNDEM_MATH_HD static Quaternion fromUnitAxisAngle(const Vec3& unitAxis, Real angle) noexcept
    {
        const Real halfAngle = 0.5 * angle;
        const Real scale = detail::sin(halfAngle);
        return {detail::cos(halfAngle), scale * unitAxis.x, scale * unitAxis.y, scale * unitAxis.z};
    }

    /**
     * Constructs a unit quaternion from an arbitrary axis and angle.
     * @param axis Rotation axis.
     * @param angle Rotation angle in radians.
     * @param result Receives the unit quaternion on success.
     * @param tolerance Minimum usable axis length.
     * @return False when the inputs are non-finite or the axis is degenerate.
     */
    FUNDEM_MATH_HD static bool tryFromAxisAngle(const Vec3& axis, Real angle, Quaternion& result, Real tolerance = defaultTolerance) noexcept
    {
        if (!isFinite(angle))
        {
            return false;
        }

        Vec3 unitAxis = axis;
        if (!tryNormalize(unitAxis, tolerance))
        {
            return false;
        }

        const Quaternion candidate = fromUnitAxisAngle(unitAxis, angle);
        if (!isFinite(candidate.w) || !isFinite(candidate.x) || !isFinite(candidate.y) || !isFinite(candidate.z))
        {
            return false;
        }

        result = candidate;
        return true;
    }

    /** Converts a rotation vector to a unit quaternion, accepting a zero vector as identity. */
    FUNDEM_MATH_HD static bool tryFromRotationVector(const Vec3& rotationVector, Quaternion& result, Real tolerance = defaultTolerance) noexcept
    {
        if (!isFinite(rotationVector))
        {
            return false;
        }

        const Real angle = norm(rotationVector);
        const Real positiveTolerance = detail::abs(tolerance);
        if (!isFinite(angle) || !isFinite(positiveTolerance))
        {
            return false;
        }

        if (angle <= positiveTolerance)
        {
            result = identity();
            return true;
        }

        return tryFromAxisAngle(rotationVector / angle, angle, result, tolerance);
    }

    FUNDEM_MATH_HD constexpr Quaternion& operator+=(const Quaternion& rhs) noexcept
    {
        w += rhs.w;
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    FUNDEM_MATH_HD constexpr Quaternion& operator-=(const Quaternion& rhs) noexcept
    {
        w -= rhs.w;
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    FUNDEM_MATH_HD constexpr Quaternion& operator*=(Real scalar) noexcept
    {
        w *= scalar;
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    FUNDEM_MATH_HD constexpr Quaternion& operator/=(Real scalar) noexcept
    {
        w /= scalar;
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }
};

static_assert(std::is_standard_layout_v<Quaternion>);
static_assert(std::is_trivially_copyable_v<Quaternion>);
static_assert(sizeof(Quaternion) == 4 * sizeof(Real));

FUNDEM_MATH_HD constexpr bool operator==(const Quaternion& lhs, const Quaternion& rhs) noexcept { return lhs.w == rhs.w && lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z; }

FUNDEM_MATH_HD constexpr bool operator!=(const Quaternion& lhs, const Quaternion& rhs) noexcept { return !(lhs == rhs); }

FUNDEM_MATH_HD constexpr Quaternion operator+(Quaternion lhs, const Quaternion& rhs) noexcept { return lhs += rhs; }

FUNDEM_MATH_HD constexpr Quaternion operator-(Quaternion lhs, const Quaternion& rhs) noexcept { return lhs -= rhs; }

FUNDEM_MATH_HD constexpr Quaternion operator-(const Quaternion& value) noexcept { return {-value.w, -value.x, -value.y, -value.z}; }

FUNDEM_MATH_HD constexpr Quaternion operator*(Quaternion value, Real scalar) noexcept { return value *= scalar; }

FUNDEM_MATH_HD constexpr Quaternion operator*(Real scalar, Quaternion value) noexcept { return value *= scalar; }

FUNDEM_MATH_HD constexpr Quaternion operator/(Quaternion value, Real scalar) noexcept { return value /= scalar; }

FUNDEM_MATH_HD constexpr Quaternion operator*(const Quaternion& lhs, const Quaternion& rhs) noexcept
{
    return {lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z,
            lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x,
            lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w};
}

FUNDEM_MATH_HD constexpr Real dot(const Quaternion& lhs, const Quaternion& rhs) noexcept { return lhs.w * rhs.w + lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }

FUNDEM_MATH_HD constexpr Real normSquared(const Quaternion& value) noexcept { return dot(value, value); }

FUNDEM_MATH_HD inline Real norm(const Quaternion& value) noexcept
{
    if (!isFinite(value.w) || !isFinite(value.x) || !isFinite(value.y) || !isFinite(value.z))
    {
        return detail::sqrt(normSquared(value));
    }

    const Real scale = detail::maxAbs(value.w, value.x, value.y, value.z);
    if (scale == 0.0)
    {
        return 0.0;
    }

    const Quaternion scaled = value / scale;
    return scale * detail::sqrt(normSquared(scaled));
}

FUNDEM_MATH_HD constexpr Quaternion conjugate(const Quaternion& value) noexcept { return {value.w, -value.x, -value.y, -value.z}; }

FUNDEM_MATH_HD inline bool isFinite(const Quaternion& value) noexcept { return isFinite(value.w) && isFinite(value.x) && isFinite(value.y) && isFinite(value.z); }

FUNDEM_MATH_HD inline bool tryNormalize(Quaternion& value, Real tolerance = defaultTolerance) noexcept
{
    if (!isFinite(value))
    {
        return false;
    }

    const Real scale = detail::maxAbs(value.w, value.x, value.y, value.z);
    const Real positiveTolerance = detail::abs(tolerance);
    if (scale == 0.0 || !isFinite(positiveTolerance))
    {
        return false;
    }

    const Quaternion scaled = value / scale;
    const Real scaledNorm = detail::sqrt(normSquared(scaled));
    if (!isFinite(scaledNorm) || scaledNorm == 0.0 || scale <= positiveTolerance / scaledNorm)
    {
        return false;
    }

    value = scaled / scaledNorm;
    return true;
}

FUNDEM_MATH_HD inline Quaternion normalizedOrIdentity(const Quaternion& value, Real tolerance = defaultTolerance) noexcept
{
    Quaternion result = value;
    return tryNormalize(result, tolerance) ? result : Quaternion::identity();
}

/** Computes the inverse of a finite, nondegenerate quaternion. */
FUNDEM_MATH_HD inline bool tryInverse(const Quaternion& value, Quaternion& result, Real tolerance = defaultTolerance) noexcept
{
    if (!isFinite(value))
    {
        return false;
    }

    const Real scale = detail::maxAbs(value.w, value.x, value.y, value.z);
    const Real positiveTolerance = detail::abs(tolerance);
    if (scale == 0.0 || !isFinite(positiveTolerance))
    {
        return false;
    }

    const Quaternion scaled = value / scale;
    const Real scaledNormSquared = normSquared(scaled);
    const Real scaledNorm = detail::sqrt(scaledNormSquared);
    if (!isFinite(scaledNorm) || scaledNorm == 0.0 || scale <= positiveTolerance / scaledNorm)
    {
        return false;
    }

    const Quaternion candidate = conjugate(scaled) / (scaledNormSquared * scale);
    if (!isFinite(candidate))
    {
        return false;
    }

    result = candidate;
    return true;
}

/**
 * Fast path: orientation must be a finite unit quaternion.
 */
FUNDEM_MATH_HD inline Vec3 rotateUnit(const Quaternion& orientation, const Vec3& vector) noexcept
{
    const Vec3 imaginary{orientation.x, orientation.y, orientation.z};
    const Vec3 twiceCross = 2.0 * cross(imaginary, vector);
    return vector + orientation.w * twiceCross + cross(imaginary, twiceCross);
}

/**
 * Fast path: orientation must be a finite unit quaternion.
 */
FUNDEM_MATH_HD inline Vec3 inverseRotateUnit(const Quaternion& orientation, const Vec3& vector) noexcept { return rotateUnit(conjugate(orientation), vector); }

/** Normalizes @p orientation and rotates @p vector when both are finite. */
FUNDEM_MATH_HD inline bool tryRotate(const Quaternion& orientation, const Vec3& vector, Vec3& result, Real tolerance = defaultTolerance) noexcept
{
    if (!isFinite(vector))
    {
        return false;
    }

    Quaternion unitOrientation = orientation;
    if (!tryNormalize(unitOrientation, tolerance))
    {
        return false;
    }

    const Vec3 candidate = rotateUnit(unitOrientation, vector);
    if (!isFinite(candidate))
    {
        return false;
    }

    result = candidate;
    return true;
}

/** Safe inverse-rotation counterpart of `tryRotate`. */
FUNDEM_MATH_HD inline bool tryInverseRotate(const Quaternion& orientation, const Vec3& vector, Vec3& result, Real tolerance = defaultTolerance) noexcept
{
    return tryRotate(conjugate(orientation), vector, result, tolerance);
}

/**
 * Fast path: orientation must be a finite unit quaternion.
 */
FUNDEM_MATH_HD constexpr Mat3 rotationMatrixUnit(const Quaternion& orientation) noexcept
{
    const Real xx = orientation.x * orientation.x;
    const Real yy = orientation.y * orientation.y;
    const Real zz = orientation.z * orientation.z;
    const Real xy = orientation.x * orientation.y;
    const Real xz = orientation.x * orientation.z;
    const Real yz = orientation.y * orientation.z;
    const Real wx = orientation.w * orientation.x;
    const Real wy = orientation.w * orientation.y;
    const Real wz = orientation.w * orientation.z;

    return {1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy), 2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx), 2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)};
}

FUNDEM_MATH_HD inline bool tryRotationMatrix(const Quaternion& orientation, Mat3& result, Real tolerance = defaultTolerance) noexcept
{
    Quaternion unitOrientation = orientation;
    if (!tryNormalize(unitOrientation, tolerance))
    {
        return false;
    }

    const Mat3 candidate = rotationMatrixUnit(unitOrientation);
    if (!isFinite(candidate))
    {
        return false;
    }

    result = candidate;
    return true;
}

namespace detail
{

template <bool WorldFrame>
FUNDEM_MATH_HD inline bool tryIntegrateAngularVelocity(const Quaternion& orientation, const Vec3& angularVelocity, Real timeStep, Quaternion& result, Real tolerance) noexcept
{
    if (!isFinite(orientation) || !isFinite(angularVelocity) || !isFinite(timeStep))
    {
        return false;
    }

    const Quaternion omega{0.0, angularVelocity.x, angularVelocity.y, angularVelocity.z};
    const Quaternion derivative = WorldFrame ? omega * orientation : orientation * omega;
    Quaternion candidate = orientation + (0.5 * timeStep) * derivative;

    if (!tryNormalize(candidate, tolerance))
    {
        return false;
    }

    result = candidate;
    return true;
}

} // namespace detail

/** Integrates a world-frame angular velocity and renormalizes the orientation. */
FUNDEM_MATH_HD inline bool tryIntegrateWorldAngularVelocity(const Quaternion& orientation, const Vec3& angularVelocity, Real timeStep, Quaternion& result, Real tolerance = defaultTolerance) noexcept
{
    return detail::tryIntegrateAngularVelocity<true>(orientation, angularVelocity, timeStep, result, tolerance);
}

/** Integrates a body-frame angular velocity and renormalizes the orientation. */
FUNDEM_MATH_HD inline bool tryIntegrateBodyAngularVelocity(const Quaternion& orientation, const Vec3& angularVelocity, Real timeStep, Quaternion& result, Real tolerance = defaultTolerance) noexcept
{
    return detail::tryIntegrateAngularVelocity<false>(orientation, angularVelocity, timeStep, result, tolerance);
}

FUNDEM_MATH_HD inline bool nearlyEqual(const Quaternion& lhs, const Quaternion& rhs, Real absoluteTolerance = defaultTolerance, Real relativeTolerance = defaultTolerance) noexcept
{
    return nearlyEqual(lhs.w, rhs.w, absoluteTolerance, relativeTolerance) && nearlyEqual(lhs.x, rhs.x, absoluteTolerance, relativeTolerance) &&
           nearlyEqual(lhs.y, rhs.y, absoluteTolerance, relativeTolerance) && nearlyEqual(lhs.z, rhs.z, absoluteTolerance, relativeTolerance);
}

FUNDEM_MATH_HD inline bool representsSameRotation(const Quaternion& lhs, const Quaternion& rhs, Real tolerance = defaultTolerance) noexcept
{
    Quaternion normalizedLhs = lhs;
    Quaternion normalizedRhs = rhs;
    if (!tryNormalize(normalizedLhs, tolerance) || !tryNormalize(normalizedRhs, tolerance))
    {
        return false;
    }

    return nearlyEqual(normalizedLhs, normalizedRhs, tolerance, tolerance) || nearlyEqual(normalizedLhs, -normalizedRhs, tolerance, tolerance);
}

} // namespace fundem::math
