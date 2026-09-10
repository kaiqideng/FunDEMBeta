/**
 * @file SPHJet.h
 * @brief Defines the geometric and temporal parameters of an SPH jet.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem
{

/** Value type describing a finite-duration cylindrical SPH jet. */
class SPHJet
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    SPHJet() noexcept = default;

    /**
     * Creates a jet from its outlet, axis, radius, speed, and emission duration.
     * If any argument or the derived pipe length is invalid, the object remains
     * in its default invalid state.
     */
    SPHJet(const Vec3& outletCenter, const Vec3& direction, Real radius, Real speed, Real duration) noexcept
    {
        SPHJet candidate;
        if (candidate.setOutletCenter(outletCenter) && candidate.setDirection(direction) && candidate.setRadius(radius) && candidate.setSpeed(speed) && candidate.setDuration(duration))
        {
            *this = candidate;
        }
    }

    const Vec3& outletCenter() const noexcept { return outletCenter_; }
    const Vec3& direction() const noexcept { return direction_; }
    Real radius() const noexcept { return radius_; }
    Real speed() const noexcept { return speed_; }
    Real duration() const noexcept { return duration_; }
    Vec3 velocity() const noexcept { return speed_ * direction_; }
    Real pipeLength() const noexcept { return speed_ * duration_; }

    bool setOutletCenter(const Vec3& value) noexcept
    {
        if (!math::isFinite(value))
        {
            return false;
        }
        outletCenter_ = value;
        return true;
    }

    bool setDirection(const Vec3& value) noexcept
    {
        Vec3 normalized = value;
        if (!math::tryNormalize(normalized))
        {
            return false;
        }
        direction_ = normalized;
        return true;
    }

    bool setRadius(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= math::defaultTolerance)
        {
            return false;
        }
        radius_ = value;
        return true;
    }

    bool setSpeed(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= math::defaultTolerance)
        {
            return false;
        }
        speed_ = value;
        return true;
    }

    bool setDuration(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= math::defaultTolerance)
        {
            return false;
        }
        duration_ = value;
        return true;
    }

    /** Returns whether all required jet values have been set. */
    bool isValid() const noexcept { return direction_ != Vec3::zero() && radius_ > 0.0 && speed_ > 0.0 && duration_ > 0.0; }

private:
    Vec3 outletCenter_{Vec3::zero()};
    Vec3 direction_{Vec3::zero()};
    Real radius_{0.0};
    Real speed_{0.0};
    Real duration_{0.0};
};

} // namespace fundem
