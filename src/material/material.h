/**
 * @file material.h
 * @brief Defines bulk elastic properties shared by all solid particle families.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "math/Constants.h"

#include <stdexcept>

namespace fundem
{

/** Bulk material properties; contact-specific coefficients belong to materialPair. */
class material
{
public:
    using Real = math::Real;

    material() = default;
    /** Constructs a complete material after validating all four properties. */
    material(Real youngModulus, Real poissonRatio, Real restitutionCoefficient, Real density) { setProperties(youngModulus, poissonRatio, restitutionCoefficient, density); }

    FUNDEM_MATH_HD Real youngModulus() const noexcept { return youngModulus_; }
    FUNDEM_MATH_HD Real poissonRatio() const noexcept { return poissonRatio_; }
    FUNDEM_MATH_HD Real restitutionCoefficient() const noexcept { return restitutionCoefficient_; }
    FUNDEM_MATH_HD Real density() const noexcept { return density_; }
    /** Reports whether the required positive elastic modulus and density were assigned. */
    FUNDEM_MATH_HD bool isValid() const noexcept { return youngModulus_ > 0.0 && density_ > 0.0; }

    void setYoungModulus(Real value)
    {
        validatePositive(value, "Young's modulus must be finite and positive.");
        youngModulus_ = value;
    }
    void setPoissonRatio(Real value)
    {
        validatePoissonRatio(value);
        poissonRatio_ = value;
    }
    void setRestitutionCoefficient(Real value)
    {
        validateRestitution(value);
        restitutionCoefficient_ = value;
    }
    void setDensity(Real value)
    {
        validatePositive(value, "Density must be finite and positive.");
        density_ = value;
    }

    /** Validates every value before replacing the complete definition. */
    void setProperties(Real youngModulus, Real poissonRatio, Real restitutionCoefficient, Real density)
    {
        validatePositive(youngModulus, "Young's modulus must be finite and positive.");
        validatePoissonRatio(poissonRatio);
        validateRestitution(restitutionCoefficient);
        validatePositive(density, "Density must be finite and positive.");
        youngModulus_ = youngModulus;
        poissonRatio_ = poissonRatio;
        restitutionCoefficient_ = restitutionCoefficient;
        density_ = density;
    }

    struct youngModulusField;
    struct poissonRatioField;
    struct restitutionCoefficientField;
    struct densityField;

private:
    /** Throws unless a required positive property is finite and positive. */
    static void validatePositive(Real value, const char* message)
    {
        if (!math::isFinite(value) || value <= 0.0)
        {
            throw std::invalid_argument(message);
        }
    }
    /** Enforces the isotropic elastic domain without excluding auxetic materials. */
    static void validatePoissonRatio(Real value)
    {
        if (!math::isFinite(value) || value <= -1.0 || value >= 0.5)
        {
            throw std::invalid_argument("Poisson ratio must be finite and within (-1, 0.5).");
        }
    }
    /** Enforces the closed interval used by restitution-based damping. */
    static void validateRestitution(Real value)
    {
        if (!math::isFinite(value) || value < 0.0 || value > 1.0)
        {
            throw std::invalid_argument("Restitution coefficient must be finite and within [0, 1].");
        }
    }

    Real youngModulus_{0.0};           ///< Young's modulus in Pa; zero denotes an incomplete material.
    Real poissonRatio_{0.0};           ///< Dimensionless Poisson ratio in (-1, 0.5).
    Real restitutionCoefficient_{0.0}; ///< Dimensionless coefficient of restitution in [0, 1].
    Real density_{0.0};                ///< Positive mass density in kg/m^3.

public:
    struct youngModulusField {
        inline static constexpr auto member = &material::youngModulus_;
    };
    struct poissonRatioField {
        inline static constexpr auto member = &material::poissonRatio_;
    };
    struct restitutionCoefficientField {
        inline static constexpr auto member = &material::restitutionCoefficient_;
    };
    struct densityField {
        inline static constexpr auto member = &material::density_;
    };

    using DeviceLayout = deviceLayout<youngModulusField, poissonRatioField, restitutionCoefficientField, densityField>;
};

using materialContainer = hostAoSDeviceSoA<material, material::DeviceLayout>;

} // namespace fundem
