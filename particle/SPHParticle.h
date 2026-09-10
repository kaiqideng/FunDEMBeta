/**
 * @file SPHParticle.h
 * @brief Defines WCSPH particle state and unified host/device storage.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "pointMass.h"

namespace fundem
{

/**
 * Weakly-compressible SPH particle. Uniform fluid properties remain solver
 * parameters; this type stores only per-particle evolving state.
 */
class SPHParticle : public pointMass
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    SPHParticle() = default;
    ~SPHParticle() override = default;
    SPHParticle(const SPHParticle&) = default;
    SPHParticle& operator=(const SPHParticle&) = default;
    SPHParticle(SPHParticle&&) noexcept = default;
    SPHParticle& operator=(SPHParticle&&) noexcept = default;

    /** Creates fluid kinematics; solver initialization assigns mass and thermodynamic state. */
    SPHParticle(const Vec3& position, const Vec3& velocity = Vec3::zero()) noexcept : pointMass(position, velocity) {}

    Real mass() const noexcept { return inverseMass() > 0.0 ? 1.0 / inverseMass() : 0.0; }
    Real density() const noexcept { return density_; }
    Real pressure() const noexcept { return pressure_; }
    Real densityRate() const noexcept { return densityRate_; }
    bool isFreeSurface() const noexcept { return freeSurface_ != 0; }
    bool isConstrained() const noexcept { return constrained_ != 0; }
    const Vec3& priorForce() const noexcept { return priorForce_; }

    bool setDensity(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= math::defaultTolerance)
        {
            return false;
        }
        density_ = value;
        return true;
    }

    bool setPressure(Real value) noexcept
    {
        if (!math::isFinite(value))
        {
            return false;
        }
        pressure_ = value;
        return true;
    }

    bool setDensityRate(Real value) noexcept
    {
        if (!math::isFinite(value))
        {
            return false;
        }
        densityRate_ = value;
        return true;
    }

    void setFreeSurface(bool value) noexcept { freeSurface_ = value ? 1 : 0; }
    void setConstrained(bool value) noexcept { constrained_ = value ? 1 : 0; }
    void setPriorForce(const Vec3& value) noexcept { priorForce_ = value; }
    /** Reports whether user-supplied kinematics and accumulated force are finite. */
    bool isValid() const noexcept { return math::isFinite(position()) && math::isFinite(velocity()) && math::isFinite(force()); }

    struct densityField;
    struct pressureField;
    struct densityRateField;
    struct freeSurfaceField;
    struct constrainedField;
    struct priorForceField;

    /** Mutable device view extending translational state with WCSPH variables. */
    struct device_type : pointMass::device_type {
        Real* density_{nullptr};     ///< Current densities.
        Real* pressure_{nullptr};    ///< Equation-of-state pressures.
        Real* densityRate_{nullptr}; ///< Material density rates.
        int* freeSurface_{nullptr};  ///< Nonzero for detected free-surface particles.
        int* constrained_{nullptr};  ///< Nonzero while the particle is subject to an inlet constraint.
        Vec3* priorForce_{nullptr};  ///< Force retained for split time integration.
    };

private:
    using pointMass::setInfiniteMass;

    Real density_{0.0};             ///< Current mass density.
    Real pressure_{0.0};            ///< Current pressure.
    Real densityRate_{0.0};         ///< Current material density derivative.
    int freeSurface_{0};            ///< Integer free-surface flag for device storage.
    int constrained_{0};            ///< Integer inlet-constraint flag for device storage.
    Vec3 priorForce_{Vec3::zero()}; ///< Force retained by the split integrator.

public:
    struct densityField {
        inline static constexpr auto member = &SPHParticle::density_;
    };
    struct pressureField {
        inline static constexpr auto member = &SPHParticle::pressure_;
    };
    struct densityRateField {
        inline static constexpr auto member = &SPHParticle::densityRate_;
    };
    struct freeSurfaceField {
        inline static constexpr auto member = &SPHParticle::freeSurface_;
    };
    struct constrainedField {
        inline static constexpr auto member = &SPHParticle::constrained_;
    };
    struct priorForceField {
        inline static constexpr auto member = &SPHParticle::priorForce_;
    };

    using DeviceLayout = concatDeviceLayoutsT<pointMass::DeviceLayout, deviceLayout<densityField, pressureField, densityRateField, freeSurfaceField, constrainedField, priorForceField>>;
};

using SPHParticleContainer = hostAoSDeviceSoA<SPHParticle, SPHParticle::DeviceLayout>;

/** Builds the mutable device view of an SPH-particle container. */
inline SPHParticle::device_type deviceFields(SPHParticleContainer& particles) noexcept
{
    SPHParticle::device_type result;
    result.position_ = particles.device<pointMass::positionField>();
    result.velocity_ = particles.device<pointMass::velocityField>();
    result.force_ = particles.device<pointMass::forceField>();
    result.inverseMass_ = particles.device<pointMass::inverseMassField>();
    result.density_ = particles.device<SPHParticle::densityField>();
    result.pressure_ = particles.device<SPHParticle::pressureField>();
    result.densityRate_ = particles.device<SPHParticle::densityRateField>();
    result.freeSurface_ = particles.device<SPHParticle::freeSurfaceField>();
    result.constrained_ = particles.device<SPHParticle::constrainedField>();
    result.priorForce_ = particles.device<SPHParticle::priorForceField>();
    result.size_ = static_cast<int>(particles.deviceSize());
    return result;
}

} // namespace fundem
