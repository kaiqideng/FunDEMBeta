/**
 * @file particle.h
 * @brief Defines spherical DEM particles and their host/device container.
 */
#pragma once

#include "rigidBody.h"
#include "material/material.h"

namespace fundem
{

/** Rigid spherical DEM particle with material-derived mass and inertia. */
class particle : public rigidBody
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;
    using Mat3 = math::Mat3;
    using Quaternion = math::Quaternion;

    particle() = default;
    ~particle() override = default;
    particle(const particle&) = default;
    particle& operator=(const particle&) = default;
    particle(particle&&) noexcept = default;
    particle& operator=(particle&&) noexcept = default;

    /** Creates sphere kinematics before radius and material are assigned. */
    particle(const Vec3& position, const Quaternion& orientation, const Vec3& velocity, const Vec3& angularVelocity) : rigidBody(position, orientation, velocity, angularVelocity) {}

    Real radius() const noexcept { return radius_; }
    int materialIndex() const noexcept { return materialIndex_; }
    /** Tests whether the stored stable material index belongs to @p container. */
    bool materialBelongsTo(const materialContainer& container) const noexcept { return materialContainer_ == &container; }

    const material& getMaterial() const noexcept { return material_; }

    /** Assigns radius and immediately refreshes mass and inertia when possible. */
    void setRadius(const Real value) noexcept
    {
        if (!setRadiusValue(value))
        {
            return;
        }
        setMassAndInertia();
    }

    /** Binds a stable material index and refreshes mass and inertia. */
    void setMaterial(const materialContainer& container, const int index) noexcept
    {
        if (index < 0 || index >= static_cast<int>(container.hostSize()))
            return;

        material_ = container.host()[index];
        materialContainer_ = &container;
        materialIndex_ = index;
        materialSet_ = true;
        setMassAndInertia();
    }

    /** Reports whether radius, material, and all kinematic values are ready for insertion. */
    bool isValid() const noexcept
    {
        return radiusSet_ && materialSet_ && math::isFinite(position()) && math::isFinite(orientation()) && math::isFinite(velocity()) && math::isFinite(angularVelocity());
    }

    struct radiusField;
    struct materialIndexField;

    /** Mutable device view extending rigid-body state with sphere metadata. */
    struct device_type : rigidBody::device_type {
        Real* radius_{nullptr};       ///< Sphere radii.
        int* materialIndex_{nullptr}; ///< Indices into the solver material container.
    };

protected:
    /** Accepts and stores a finite positive radius without recomputing mass properties. */
    bool setRadiusValue(const Real value) noexcept
    {
        if (!math::isFinite(value) || value <= 0.0)
        {
            return false;
        }
        radius_ = value;
        radiusSet_ = true;
        return true;
    }

private:
    /** Recomputes homogeneous-sphere mass and centroidal inertia from radius and material. */
    virtual void setMassAndInertia() noexcept
    {
        const Real mass = (4.0 / 3.0) * math::pi * radius_ * radius_ * radius_ * material_.density();
        const Real inertiaValue = (2.0 / 5.0) * mass * radius_ * radius_;
        const Mat3 inertiaTensor = Mat3::identity() * inertiaValue;
        rigidBody::setMassAndInertia(mass, inertiaTensor);
    }

    material material_;                                   ///< Host-side material snapshot.
    const materialContainer* materialContainer_{nullptr}; ///< Container that owns `materialIndex_`.
    bool radiusSet_{false};                               ///< Whether a positive radius was accepted.
    bool materialSet_{false};                             ///< Whether a material was bound.

    Real radius_{0.0};      ///< Bounding sphere radius.
    int materialIndex_{-1}; ///< Stable material-container index.

public:
    struct radiusField {
        inline static constexpr auto member = &particle::radius_;
    };
    struct materialIndexField {
        inline static constexpr auto member = &particle::materialIndex_;
    };

    using DeviceLayout = concatDeviceLayoutsT<rigidBody::DeviceLayout, deviceLayout<radiusField, materialIndexField>>;
};

using particleContainer = hostAoSDeviceSoA<particle, particle::DeviceLayout>;

/** Builds a mutable device view for any storage whose host type derives from `particle`. */
template <class ParticleStorage> particle::device_type particleDeviceFields(ParticleStorage& particles) noexcept
{
    particle::device_type result;
    result.position_ = particles.template device<pointMass::positionField>();
    result.velocity_ = particles.template device<pointMass::velocityField>();
    result.force_ = particles.template device<pointMass::forceField>();
    result.inverseMass_ = particles.template device<pointMass::inverseMassField>();
    result.orientation_ = particles.template device<rigidBody::orientationField>();
    result.angularVelocity_ = particles.template device<rigidBody::angularVelocityField>();
    result.torque_ = particles.template device<rigidBody::torqueField>();
    result.inertiaTensor_ = particles.template device<rigidBody::inertiaTensorField>();
    result.inverseInertiaTensor_ = particles.template device<rigidBody::inverseInertiaTensorField>();
    result.radius_ = particles.template device<particle::radiusField>();
    result.materialIndex_ = particles.template device<particle::materialIndexField>();
    result.size_ = static_cast<int>(particles.deviceSize());
    return result;
}

/** Builds the mutable device view of a spherical-particle container. */
inline particle::device_type deviceFields(particleContainer& particles) noexcept { return particleDeviceFields(particles); }

} // namespace fundem
