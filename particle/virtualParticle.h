/**
 * @file virtualParticle.h
 * @brief Defines internally generated SPH boundary samples and wall-reaction state.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "pointMass.h"

namespace fundem
{

/**
 * Internally generated boundary quadrature point attached to one level-set
 * rigid body for SPH wall interaction and reaction transfer.
 */
class virtualParticle : public pointMass
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    virtualParticle() = default;
    ~virtualParticle() override = default;
    virtualParticle(const virtualParticle&) = default;
    virtualParticle& operator=(const virtualParticle&) = default;
    virtualParticle(virtualParticle&&) noexcept = default;
    virtualParticle& operator=(virtualParticle&&) noexcept = default;

    /** Creates one owner-local boundary quadrature sample. */
    virtualParticle(const Vec3& localPosition, const Vec3& localNormal, Real volume, int ownerLSParticleIndex) noexcept { setLocalValues(localPosition, localNormal, volume, ownerLSParticleIndex); }

    const Vec3& localPosition() const noexcept { return localPosition_; }
    const Vec3& localNormal() const noexcept { return localNormal_; }
    Real volume() const noexcept { return volume_; }
    int ownerLSParticleIndex() const noexcept { return ownerLSParticleIndex_; }
    const Vec3& normal() const noexcept { return normal_; }
    const Vec3& acceleration() const noexcept { return acceleration_; }
    const Vec3& priorForce() const noexcept { return priorForce_; }

    /**
     * Assigns immutable owner-local geometry values.
     * @param localPosition Position in the owner's body frame.
     * @param localNormal Outward normal in the owner's body frame.
     * @param volume Boundary quadrature volume.
     * @param ownerLSParticleIndex Stable owner index.
     * @return True when every value is usable.
     */
    bool setLocalValues(const Vec3& localPosition, const Vec3& localNormal, Real volume, int ownerLSParticleIndex) noexcept
    {
        const Vec3 unitNormal = math::normalizedOrZero(localNormal);
        if (!math::isFinite(localPosition) || unitNormal == Vec3::zero() || !math::isFinite(volume) || volume <= 0.0 || ownerLSParticleIndex < 0)
        {
            return false;
        }
        localPosition_ = localPosition;
        localNormal_ = unitNormal;
        volume_ = volume;
        ownerLSParticleIndex_ = ownerLSParticleIndex;
        return true;
    }

    /** Normalizes and assigns the current world-frame wall normal. */
    bool setNormal(const Vec3& value) noexcept
    {
        const Vec3 unitNormal = math::normalizedOrZero(value);
        if (unitNormal == Vec3::zero())
        {
            return false;
        }
        normal_ = unitNormal;
        return true;
    }
    void setAcceleration(const Vec3& value) noexcept { acceleration_ = value; }
    void setPriorForce(const Vec3& value) noexcept { priorForce_ = value; }

    struct localPositionField;
    struct localNormalField;
    struct volumeField;
    struct ownerLSParticleIndexField;
    struct normalField;
    struct accelerationField;
    struct priorForceField;

    /** Mutable device view of generated wall-sample state. */
    struct device_type : pointMass::device_type {
        Vec3* localPosition_{nullptr};       ///< Body-frame sample positions.
        Vec3* localNormal_{nullptr};         ///< Body-frame outward normals.
        Real* volume_{nullptr};              ///< Boundary quadrature volumes.
        int* ownerLSParticleIndex_{nullptr}; ///< Owner rigid-body indices.
        Vec3* normal_{nullptr};              ///< Current world outward normals.
        Vec3* acceleration_{nullptr};        ///< Current wall accelerations.
        Vec3* priorForce_{nullptr};          ///< Previous fluid reaction forces.
    };

private:
    using pointMass::setInfiniteMass;
    using pointMass::setMass;
    using pointMass::updatePosition;
    using pointMass::updateVelocity;

    Vec3 localPosition_{Vec3::zero()}; ///< Body-frame sample position.
    Vec3 localNormal_{Vec3::zero()};   ///< Body-frame outward unit normal.
    Real volume_{0.0};                 ///< Boundary quadrature volume.
    int ownerLSParticleIndex_{-1};     ///< Stable owner rigid-body index.
    Vec3 normal_{Vec3::zero()};        ///< Current world outward unit normal.
    Vec3 acceleration_{Vec3::zero()};  ///< Current wall-point acceleration.
    Vec3 priorForce_{Vec3::zero()};    ///< Previous fluid reaction force.

public:
    struct localPositionField {
        inline static constexpr auto member = &virtualParticle::localPosition_;
    };
    struct localNormalField {
        inline static constexpr auto member = &virtualParticle::localNormal_;
    };
    struct volumeField {
        inline static constexpr auto member = &virtualParticle::volume_;
    };
    struct ownerLSParticleIndexField {
        inline static constexpr auto member = &virtualParticle::ownerLSParticleIndex_;
    };
    struct normalField {
        inline static constexpr auto member = &virtualParticle::normal_;
    };
    struct accelerationField {
        inline static constexpr auto member = &virtualParticle::acceleration_;
    };
    struct priorForceField {
        inline static constexpr auto member = &virtualParticle::priorForce_;
    };

    using DeviceLayout =
        concatDeviceLayoutsT<pointMass::DeviceLayout, deviceLayout<localPositionField, localNormalField, volumeField, ownerLSParticleIndexField, normalField, accelerationField, priorForceField>>;
};

using virtualParticleContainer = hostAoSDeviceSoA<virtualParticle, virtualParticle::DeviceLayout>;

/** Builds the mutable device view of a virtual-particle container. */
inline virtualParticle::device_type deviceFields(virtualParticleContainer& particles) noexcept
{
    virtualParticle::device_type result;
    result.position_ = particles.device<pointMass::positionField>();
    result.velocity_ = particles.device<pointMass::velocityField>();
    result.force_ = particles.device<pointMass::forceField>();
    result.inverseMass_ = particles.device<pointMass::inverseMassField>();
    result.localPosition_ = particles.device<virtualParticle::localPositionField>();
    result.localNormal_ = particles.device<virtualParticle::localNormalField>();
    result.volume_ = particles.device<virtualParticle::volumeField>();
    result.ownerLSParticleIndex_ = particles.device<virtualParticle::ownerLSParticleIndexField>();
    result.normal_ = particles.device<virtualParticle::normalField>();
    result.acceleration_ = particles.device<virtualParticle::accelerationField>();
    result.priorForce_ = particles.device<virtualParticle::priorForceField>();
    result.size_ = static_cast<int>(particles.deviceSize());
    return result;
}

} // namespace fundem
