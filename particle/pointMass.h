/**
 * @file pointMass.h
 * @brief Defines translational point-mass state shared by DEM and SPH particles.
 */
#pragma once

#include "data/ContainerObject.h"
#include "data/DeviceLayout.h"
#include "math/Vector3.h"

namespace fundem
{

/**
 * Translational state shared by rigid bodies, SPH particles, and internal
 * boundary samples.
 */
class pointMass : public containerObject
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    pointMass() = default;
    virtual ~pointMass() = default;
    pointMass(const pointMass&) = default;
    pointMass& operator=(const pointMass&) = default;
    pointMass(pointMass&&) noexcept = default;
    pointMass& operator=(pointMass&&) noexcept = default;

    /** Creates translational state at @p position with @p velocity and infinite mass. */
    pointMass(const Vec3& position, const Vec3& velocity) noexcept : position_(position), velocity_(velocity) {}

    const Vec3& position() const noexcept { return position_; }
    const Vec3& velocity() const noexcept { return velocity_; }
    const Vec3& force() const noexcept { return force_; }
    Real inverseMass() const noexcept { return inverseMass_; }

    void setPosition(const Vec3& value) noexcept { position_ = value; }
    void setVelocity(const Vec3& value) noexcept { velocity_ = value; }
    void setForce(const Vec3& value) noexcept { force_ = value; }
    void addForce(const Vec3& value) noexcept { force_ += value; }
    void setMass(Real mass) noexcept { inverseMass_ = math::isFinite(mass) && mass > math::defaultTolerance ? 1.0 / mass : 0.0; }
    void setInfiniteMass() noexcept { inverseMass_ = 0.0; }

    /** Advances position using the current velocity and @p timeStep. */
    void updatePosition(Real timeStep) noexcept;
    /** Advances velocity from accumulated force and gravity. */
    void updateVelocity(const Vec3& gravity, Real timeStep) noexcept;

    struct positionField;
    struct velocityField;
    struct forceField;
    struct inverseMassField;

    /** Mutable device structure-of-arrays view of translational state. */
    struct device_type {
        Vec3* position_{nullptr};    ///< World positions.
        Vec3* velocity_{nullptr};    ///< World linear velocities.
        Vec3* force_{nullptr};       ///< Accumulated world forces.
        Real* inverseMass_{nullptr}; ///< Inverse masses; zero denotes infinite mass.
        int size_{0};                ///< Logical particle count.
    };

private:
    Vec3 position_{Vec3::zero()}; ///< World position.
    Vec3 velocity_{Vec3::zero()}; ///< World linear velocity.
    Vec3 force_{Vec3::zero()};    ///< Accumulated world force.
    Real inverseMass_{0.0};       ///< Reciprocal mass; zero denotes infinite mass.

public:
    struct positionField {
        inline static constexpr auto member = &pointMass::position_;
    };
    struct velocityField {
        inline static constexpr auto member = &pointMass::velocity_;
    };
    struct forceField {
        inline static constexpr auto member = &pointMass::force_;
    };
    struct inverseMassField {
        inline static constexpr auto member = &pointMass::inverseMass_;
    };

    using DeviceLayout = deviceLayout<positionField, velocityField, forceField, inverseMassField>;
};

} // namespace fundem
