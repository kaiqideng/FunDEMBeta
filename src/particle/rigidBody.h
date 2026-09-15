/**
 * @file rigidBody.h
 * @brief Extends point-mass state with orientation, torque, and rotational inertia.
 */
#pragma once

#include "math/Quaternion.h"
#include "pointMass.h"

namespace fundem
{

/** Point mass extended with unit orientation and rotational dynamics. */
class rigidBody : public pointMass
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;
    using Mat3 = math::Mat3;
    using Quaternion = math::Quaternion;

    rigidBody() = default;
    ~rigidBody() override = default;
    rigidBody(const rigidBody&) = default;
    rigidBody& operator=(const rigidBody&) = default;
    rigidBody(rigidBody&&) noexcept = default;
    rigidBody& operator=(rigidBody&&) noexcept = default;

    /** Creates rigid kinematics with normalized orientation and zero mass/inertia. */
    rigidBody(const Vec3& position, const Quaternion& orientation, const Vec3& velocity, const Vec3& angularVelocity) : pointMass(position, velocity), angularVelocity_(angularVelocity)
    {
        setOrientation(orientation);
    }

    const Quaternion& orientation() const noexcept { return orientation_; }
    const Vec3& angularVelocity() const noexcept { return angularVelocity_; }
    const Vec3& torque() const noexcept { return torque_; }
    const Mat3& inertiaTensor() const noexcept { return inertiaTensor_; }
    const Mat3& inverseInertiaTensor() const noexcept { return inverseInertiaTensor_; }

    void setOrientation(const Quaternion& value) noexcept
    {
        Quaternion unitOrientation = value;
        if (math::tryNormalize(unitOrientation))
        {
            orientation_ = unitOrientation;
        }
    }
    void setAngularVelocity(const Vec3& value) noexcept { angularVelocity_ = value; }
    void setTorque(const Vec3& value) noexcept { torque_ = value; }
    void addTorque(const Vec3& value) noexcept { torque_ += value; }
    void setInfiniteMass() noexcept { setMassAndInertia(0.0, Mat3::zero()); }
    /**
     * Assigns finite mass and centroidal body-frame inertia together.
     * Invalid, non-positive, or singular input produces an infinite-mass body
     * with zero inertia tensors.
     * @param mass Body mass.
     * @param inertiaTensor Centroidal inertia tensor in the body frame.
     */
    void setMassAndInertia(const Real mass, const Mat3& inertiaTensor) noexcept
    {
        if (math::isFinite(mass) && mass > math::defaultTolerance && math::tryInverse(inertiaTensor, inverseInertiaTensor_))
        {
            pointMass::setMass(mass);
            inertiaTensor_ = inertiaTensor;
        }
        else
        {
            pointMass::setInfiniteMass();
            inertiaTensor_ = Mat3::zero();
            inverseInertiaTensor_ = Mat3::zero();
        }
    }

    /** Advances world position and unit orientation. */
    void updatePositionAndOrientation(Real timeStep) noexcept;
    /** Advances linear and angular velocity from accumulated loads. */
    void updateVelocityAndAngularVelocity(const Vec3& gravity, Real timeStep) noexcept;

    struct orientationField;
    struct angularVelocityField;
    struct torqueField;
    struct inertiaTensorField;
    struct inverseInertiaTensorField;

    /** Mutable device view extending `pointMass::device_type` with rotation state. */
    struct device_type : pointMass::device_type {
        Quaternion* orientation_{nullptr};    ///< Unit world-from-body orientations.
        Vec3* angularVelocity_{nullptr};      ///< World angular velocities.
        Vec3* torque_{nullptr};               ///< Accumulated world torques.
        Mat3* inertiaTensor_{nullptr};        ///< Body-frame inertia tensors.
        Mat3* inverseInertiaTensor_{nullptr}; ///< Inverse body-frame inertia tensors.
    };

private:
    using pointMass::setMass;
    using pointMass::updatePosition;
    using pointMass::updateVelocity;

    Quaternion orientation_{Quaternion::identity()}; ///< Unit world-from-body orientation.
    Vec3 angularVelocity_{Vec3::zero()};             ///< World angular velocity.
    Vec3 torque_{Vec3::zero()};                      ///< Accumulated world torque.
    Mat3 inertiaTensor_{Mat3::zero()};               ///< Centroidal body-frame inertia.
    Mat3 inverseInertiaTensor_{Mat3::zero()};        ///< Inverse centroidal body-frame inertia.

public:
    struct orientationField {
        inline static constexpr auto member = &rigidBody::orientation_;
    };
    struct angularVelocityField {
        inline static constexpr auto member = &rigidBody::angularVelocity_;
    };
    struct torqueField {
        inline static constexpr auto member = &rigidBody::torque_;
    };
    struct inertiaTensorField {
        inline static constexpr auto member = &rigidBody::inertiaTensor_;
    };
    struct inverseInertiaTensorField {
        inline static constexpr auto member = &rigidBody::inverseInertiaTensor_;
    };

    using DeviceLayout = concatDeviceLayoutsT<pointMass::DeviceLayout, deviceLayout<orientationField, angularVelocityField, torqueField, inertiaTensorField, inverseInertiaTensorField>>;
};

} // namespace fundem
