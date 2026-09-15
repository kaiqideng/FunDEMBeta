/**
 * @file material.h
 * @brief Defines validated contact materials and unified host/device storage.
 */
#pragma once

#include "data/HostAoSDeviceSoA.h"
#include "math/Constants.h"

#include <stdexcept>
#include <string>

namespace fundem
{

/** Distinguishes ordinary stiffness values from level-set stiffness per unit area. */
enum class materialType : int
{
    standard,
    levelSet
};

/**
 * Validated contact parameters and density shared by spherical and level-set
 * particles. A level-set specialization changes semantics without adding storage.
 */
class material
{
public:
    using Real = math::Real;

    material() = default;

    /**
     * Constructs a complete standard material.
     * @throws std::invalid_argument when a coefficient is outside its physical
     * domain or is not finite.
     */
    material(Real normalStiffness,
             Real slidingStiffness,
             Real rollingStiffness,
             Real torsionalStiffness,
             Real slidingFrictionCoefficient,
             Real rollingFrictionCoefficient,
             Real torsionalFrictionCoefficient,
             Real restitutionCoefficient,
             Real density)
    {
        setProperties(normalStiffness,
                      slidingStiffness,
                      rollingStiffness,
                      torsionalStiffness,
                      slidingFrictionCoefficient,
                      rollingFrictionCoefficient,
                      torsionalFrictionCoefficient,
                      restitutionCoefficient,
                      density);
    }

    Real normalStiffness() const noexcept { return normalStiffness_; }

    Real slidingStiffness() const noexcept { return slidingStiffness_; }

    Real rollingStiffness() const noexcept { return rollingStiffness_; }

    Real torsionalStiffness() const noexcept { return torsionalStiffness_; }

    Real slidingFrictionCoefficient() const noexcept { return slidingFrictionCoefficient_; }

    Real rollingFrictionCoefficient() const noexcept { return rollingFrictionCoefficient_; }

    Real torsionalFrictionCoefficient() const noexcept { return torsionalFrictionCoefficient_; }

    Real restitutionCoefficient() const noexcept { return restitutionCoefficient_; }

    Real density() const noexcept { return density_; }

    materialType type() const noexcept { return type_; }

    bool isLevelSet() const noexcept { return type_ == materialType::levelSet; }

    /** Reports whether the required positive density has been assigned. */
    bool isValid() const noexcept { return densitySet_; }

    void setNormalStiffness(Real value)
    {
        validateNonNegative(value, "Normal stiffness");
        normalStiffness_ = value;
    }

    void setSlidingStiffness(Real value)
    {
        validateNonNegative(value, "Sliding stiffness");
        slidingStiffness_ = value;
    }

    void setRollingStiffness(Real value)
    {
        validateNonNegative(value, "Rolling stiffness");
        rollingStiffness_ = value;
    }

    void setTorsionalStiffness(Real value)
    {
        validateNonNegative(value, "Torsional stiffness");
        torsionalStiffness_ = value;
    }

    void setSlidingFrictionCoefficient(Real value)
    {
        validateNonNegative(value, "Sliding friction coefficient");
        slidingFrictionCoefficient_ = value;
    }

    void setRollingFrictionCoefficient(Real value)
    {
        validateNonNegative(value, "Rolling friction coefficient");
        rollingFrictionCoefficient_ = value;
    }

    void setTorsionalFrictionCoefficient(Real value)
    {
        validateNonNegative(value, "Torsional friction coefficient");
        torsionalFrictionCoefficient_ = value;
    }

    void setRestitutionCoefficient(Real value)
    {
        validateRestitution(value);
        restitutionCoefficient_ = value;
    }

    void setDensity(Real value)
    {
        validatePositive(value, "Density");
        density_ = value;
        densitySet_ = true;
    }

    /**
     * Validates and replaces the complete material definition.
     *
     * Stiffnesses and friction coefficients must be finite and non-negative,
     * restitution must lie in [0,1], and density must be finite and positive.
     */
    void setProperties(Real normalStiffness,
                       Real slidingStiffness,
                       Real rollingStiffness,
                       Real torsionalStiffness,
                       Real slidingFrictionCoefficient,
                       Real rollingFrictionCoefficient,
                       Real torsionalFrictionCoefficient,
                       Real restitutionCoefficient,
                       Real density)
    {
        validateNonNegative(normalStiffness, "Normal stiffness");
        validateNonNegative(slidingStiffness, "Sliding stiffness");
        validateNonNegative(rollingStiffness, "Rolling stiffness");
        validateNonNegative(torsionalStiffness, "Torsional stiffness");
        validateNonNegative(slidingFrictionCoefficient, "Sliding friction coefficient");
        validateNonNegative(rollingFrictionCoefficient, "Rolling friction coefficient");
        validateNonNegative(torsionalFrictionCoefficient, "Torsional friction coefficient");
        validateRestitution(restitutionCoefficient);
        validatePositive(density, "Density");

        normalStiffness_ = normalStiffness;
        slidingStiffness_ = slidingStiffness;
        rollingStiffness_ = rollingStiffness;
        torsionalStiffness_ = torsionalStiffness;
        slidingFrictionCoefficient_ = slidingFrictionCoefficient;
        rollingFrictionCoefficient_ = rollingFrictionCoefficient;
        torsionalFrictionCoefficient_ = torsionalFrictionCoefficient;
        restitutionCoefficient_ = restitutionCoefficient;
        density_ = density;
        densitySet_ = true;
    }

    struct normalStiffnessField;
    struct slidingStiffnessField;
    struct rollingStiffnessField;
    struct torsionalStiffnessField;
    struct slidingFrictionCoefficientField;
    struct rollingFrictionCoefficientField;
    struct torsionalFrictionCoefficientField;
    struct restitutionCoefficientField;
    struct densityField;
    struct typeField;

protected:
    /** Constructs an incomplete material carrying the requested runtime semantics. */
    explicit material(materialType type) noexcept : type_(type) {}

private:
    bool densitySet_{false}; ///< Whether a positive density has been assigned.

    /** Tests the common domain used by stiffness and friction coefficients. */
    static bool isNonNegativeFinite(Real value) noexcept { return math::isFinite(value) && value >= 0.0; }

    /** Throws unless @p value is finite and non-negative. */
    static void validateNonNegative(Real value, const char* name)
    {
        if (!isNonNegativeFinite(value))
        {
            throw std::invalid_argument(std::string(name) + " must be finite and non-negative.");
        }
    }

    /** Throws unless @p value is finite and strictly positive. */
    static void validatePositive(Real value, const char* name)
    {
        if (!math::isFinite(value) || value <= 0.0)
        {
            throw std::invalid_argument(std::string(name) + " must be finite and positive.");
        }
    }

    /** Throws unless restitution lies in the closed interval `[0, 1]`. */
    static void validateRestitution(Real value)
    {
        if (!math::isFinite(value) || value < 0.0 || value > 1.0)
        {
            throw std::invalid_argument("Restitution coefficient must be finite and within [0, 1].");
        }
    }

    Real normalStiffness_{0.0};                 ///< Normal stiffness, or stiffness per area for level-set material.
    Real slidingStiffness_{0.0};                ///< Sliding/shear stiffness in the material's native convention.
    Real rollingStiffness_{0.0};                ///< Rolling stiffness.
    Real torsionalStiffness_{0.0};              ///< Torsional stiffness.
    Real slidingFrictionCoefficient_{0.0};      ///< Coulomb sliding coefficient.
    Real rollingFrictionCoefficient_{0.0};      ///< Coulomb rolling coefficient.
    Real torsionalFrictionCoefficient_{0.0};    ///< Coulomb torsional coefficient.
    Real restitutionCoefficient_{0.0};          ///< Normal coefficient of restitution in [0,1].
    Real density_{0.0};                         ///< Mass density.
    materialType type_{materialType::standard}; ///< Runtime material semantics.

public:
    struct normalStiffnessField {
        inline static constexpr auto member = &material::normalStiffness_;
    };
    struct slidingStiffnessField {
        inline static constexpr auto member = &material::slidingStiffness_;
    };
    struct rollingStiffnessField {
        inline static constexpr auto member = &material::rollingStiffness_;
    };
    struct torsionalStiffnessField {
        inline static constexpr auto member = &material::torsionalStiffness_;
    };
    struct slidingFrictionCoefficientField {
        inline static constexpr auto member = &material::slidingFrictionCoefficient_;
    };
    struct rollingFrictionCoefficientField {
        inline static constexpr auto member = &material::rollingFrictionCoefficient_;
    };
    struct torsionalFrictionCoefficientField {
        inline static constexpr auto member = &material::torsionalFrictionCoefficient_;
    };
    struct restitutionCoefficientField {
        inline static constexpr auto member = &material::restitutionCoefficient_;
    };
    struct densityField {
        inline static constexpr auto member = &material::density_;
    };
    struct typeField {
        inline static constexpr auto member = &material::type_;
    };

    using DeviceLayout = deviceLayout<normalStiffnessField,
                                      slidingStiffnessField,
                                      rollingStiffnessField,
                                      torsionalStiffnessField,
                                      slidingFrictionCoefficientField,
                                      rollingFrictionCoefficientField,
                                      torsionalFrictionCoefficientField,
                                      restitutionCoefficientField,
                                      densityField,
                                      typeField>;
};

using materialContainer = hostAoSDeviceSoA<material, material::DeviceLayout>;

} // namespace fundem
