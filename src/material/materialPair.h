/**
 * @file materialPair.h
 * @brief Defines symmetric contact coefficients and append-stable material-pair storage.
 */
#pragma once

#include "material.h"

#include <limits>
#include <utility>

namespace fundem
{

/** Plain contact coefficients with one linear-law switch and one explicit-damping switch. */
struct materialPairCoefficients {
    math::Real normalStiffnessPerUnitArea_;   ///< LS normal stiffness per nodal area, N/m^3; -1 means unassigned.
    math::Real slidingStiffnessPerUnitArea_;  ///< LS sliding stiffness per nodal area, N/m^3; -1 means unassigned.
    math::Real normalStiffness_;              ///< Linear normal stiffness, N/m.
    math::Real slidingStiffness_;             ///< Linear sliding stiffness, N/m.
    math::Real rollingStiffness_;             ///< Length-scaled rolling stiffness, N/m.
    math::Real torsionalStiffness_;           ///< Length-scaled torsional stiffness, N/m.
    math::Real slidingFrictionCoefficient_;   ///< Dimensionless sliding friction.
    math::Real rollingFrictionCoefficient_;   ///< Dimensionless rolling friction.
    math::Real torsionalFrictionCoefficient_; ///< Dimensionless torsional friction.
    math::Real normalDampingCoefficient_;     ///< Normal dashpot coefficient, N s/m.
    math::Real slidingDampingCoefficient_;    ///< Sliding dashpot coefficient, N s/m.
    math::Real rollingDampingCoefficient_;    ///< Length-scaled rolling dashpot coefficient, N s/m.
    math::Real torsionalDampingCoefficient_;  ///< Length-scaled torsional dashpot coefficient, N s/m.
    bool linearStiffnessEnabled_;             ///< Uses all four configured linear stiffnesses instead of Hertz-Mindlin defaults.
    bool dampingEnabled_;                     ///< Uses all four configured dashpots instead of restitution-derived damping.
};

/** One unordered material pair, including equal-material self-pairs. */
class materialPair
{
public:
    using Real = math::Real;

    materialPair() = default;
    /** Stores nonnegative material indices in ascending order. */
    materialPair(int masterMaterialIndex, int slaveMaterialIndex)
    {
        if (masterMaterialIndex < 0 || slaveMaterialIndex < 0)
        {
            throw std::invalid_argument("Material-pair indices must be non-negative.");
        }
        masterMaterialIndex_ = masterMaterialIndex <= slaveMaterialIndex ? masterMaterialIndex : slaveMaterialIndex;
        slaveMaterialIndex_ = masterMaterialIndex <= slaveMaterialIndex ? slaveMaterialIndex : masterMaterialIndex;
    }

    FUNDEM_MATH_HD int masterMaterialIndex() const noexcept { return masterMaterialIndex_; }
    FUNDEM_MATH_HD int slaveMaterialIndex() const noexcept { return slaveMaterialIndex_; }
    FUNDEM_MATH_HD Real normalStiffnessPerUnitArea() const noexcept { return normalStiffnessPerUnitArea_; }
    FUNDEM_MATH_HD Real slidingStiffnessPerUnitArea() const noexcept { return slidingStiffnessPerUnitArea_; }
    FUNDEM_MATH_HD Real normalStiffness() const noexcept { return normalStiffness_; }
    FUNDEM_MATH_HD Real slidingStiffness() const noexcept { return slidingStiffness_; }
    FUNDEM_MATH_HD Real rollingStiffness() const noexcept { return rollingStiffness_; }
    FUNDEM_MATH_HD Real torsionalStiffness() const noexcept { return torsionalStiffness_; }
    FUNDEM_MATH_HD Real slidingFrictionCoefficient() const noexcept { return slidingFrictionCoefficient_; }
    FUNDEM_MATH_HD Real rollingFrictionCoefficient() const noexcept { return rollingFrictionCoefficient_; }
    FUNDEM_MATH_HD Real torsionalFrictionCoefficient() const noexcept { return torsionalFrictionCoefficient_; }
    FUNDEM_MATH_HD Real normalDampingCoefficient() const noexcept { return normalDampingCoefficient_; }
    FUNDEM_MATH_HD Real slidingDampingCoefficient() const noexcept { return slidingDampingCoefficient_; }
    FUNDEM_MATH_HD Real rollingDampingCoefficient() const noexcept { return rollingDampingCoefficient_; }
    FUNDEM_MATH_HD Real torsionalDampingCoefficient() const noexcept { return torsionalDampingCoefficient_; }
    FUNDEM_MATH_HD bool linearStiffnessEnabled() const noexcept { return linearStiffnessEnabled_; }
    FUNDEM_MATH_HD bool dampingEnabled() const noexcept { return dampingEnabled_; }
    void setLinearStiffnessEnabled(bool enabled) noexcept { linearStiffnessEnabled_ = enabled; }
    void setDampingEnabled(bool enabled) noexcept { dampingEnabled_ = enabled; }

    /** Reports whether both required LS stiffness-per-area values were assigned, including explicit zero. */
    FUNDEM_MATH_HD bool hasLevelSetStiffness() const noexcept { return normalStiffnessPerUnitArea_ >= 0.0 && slidingStiffnessPerUnitArea_ >= 0.0; }

    /** Sets required LS coefficients independently of the two contact-law switches. */
    void setLevelSetStiffness(Real normal, Real sliding)
    {
        validateNonNegative(normal);
        validateNonNegative(sliding);
        normalStiffnessPerUnitArea_ = normal;
        slidingStiffnessPerUnitArea_ = sliding;
    }
    /** Sets and enables all four effective linear stiffnesses. */
    void setLinearStiffness(Real normal, Real sliding, Real rolling = 0.0, Real torsional = 0.0)
    {
        validateNonNegative(normal);
        validateNonNegative(sliding);
        validateNonNegative(rolling);
        validateNonNegative(torsional);
        normalStiffness_ = normal;
        slidingStiffness_ = sliding;
        rollingStiffness_ = rolling;
        torsionalStiffness_ = torsional;
        linearStiffnessEnabled_ = true;
    }
    /** Replaces the three dimensionless friction coefficients without changing either switch. */
    void setFrictionCoefficients(Real sliding, Real rolling = 0.0, Real torsional = 0.0)
    {
        validateNonNegative(sliding);
        validateNonNegative(rolling);
        validateNonNegative(torsional);
        slidingFrictionCoefficient_ = sliding;
        rollingFrictionCoefficient_ = rolling;
        torsionalFrictionCoefficient_ = torsional;
    }
    /** Sets and enables all four directional dashpots; enabled zero explicitly suppresses damping. */
    void setDampingCoefficients(Real normal, Real sliding, Real rolling, Real torsional)
    {
        validateNonNegative(normal);
        validateNonNegative(sliding);
        validateNonNegative(rolling);
        validateNonNegative(torsional);
        normalDampingCoefficient_ = normal;
        slidingDampingCoefficient_ = sliding;
        rollingDampingCoefficient_ = rolling;
        torsionalDampingCoefficient_ = torsional;
        dampingEnabled_ = true;
    }
    void setNormalStiffnessPerUnitArea(Real value)
    {
        validateNonNegative(value);
        normalStiffnessPerUnitArea_ = value;
    }
    void setSlidingStiffnessPerUnitArea(Real value)
    {
        validateNonNegative(value);
        slidingStiffnessPerUnitArea_ = value;
    }
    void setNormalStiffness(Real value)
    {
        validateNonNegative(value);
        normalStiffness_ = value;
        linearStiffnessEnabled_ = true;
    }
    void setSlidingStiffness(Real value)
    {
        validateNonNegative(value);
        slidingStiffness_ = value;
        linearStiffnessEnabled_ = true;
    }
    void setRollingStiffness(Real value)
    {
        validateNonNegative(value);
        rollingStiffness_ = value;
        linearStiffnessEnabled_ = true;
    }
    void setTorsionalStiffness(Real value)
    {
        validateNonNegative(value);
        torsionalStiffness_ = value;
        linearStiffnessEnabled_ = true;
    }
    void setNormalDampingCoefficient(Real value)
    {
        validateNonNegative(value);
        normalDampingCoefficient_ = value;
        dampingEnabled_ = true;
    }
    void setSlidingDampingCoefficient(Real value)
    {
        validateNonNegative(value);
        slidingDampingCoefficient_ = value;
        dampingEnabled_ = true;
    }
    void setRollingDampingCoefficient(Real value)
    {
        validateNonNegative(value);
        rollingDampingCoefficient_ = value;
        dampingEnabled_ = true;
    }
    void setTorsionalDampingCoefficient(Real value)
    {
        validateNonNegative(value);
        torsionalDampingCoefficient_ = value;
        dampingEnabled_ = true;
    }

    /** Returns coefficients and both group switches, without material identity. */
    FUNDEM_MATH_HD materialPairCoefficients coefficients() const noexcept
    {
        return {normalStiffnessPerUnitArea_,
                slidingStiffnessPerUnitArea_,
                normalStiffness_,
                slidingStiffness_,
                rollingStiffness_,
                torsionalStiffness_,
                slidingFrictionCoefficient_,
                rollingFrictionCoefficient_,
                torsionalFrictionCoefficient_,
                normalDampingCoefficient_,
                slidingDampingCoefficient_,
                rollingDampingCoefficient_,
                torsionalDampingCoefficient_,
                linearStiffnessEnabled_,
                dampingEnabled_};
    }

    struct masterMaterialIndexField;
    struct slaveMaterialIndexField;
    struct normalStiffnessPerUnitAreaField;
    struct slidingStiffnessPerUnitAreaField;
    struct normalStiffnessField;
    struct slidingStiffnessField;
    struct rollingStiffnessField;
    struct torsionalStiffnessField;
    struct slidingFrictionCoefficientField;
    struct rollingFrictionCoefficientField;
    struct torsionalFrictionCoefficientField;
    struct normalDampingCoefficientField;
    struct slidingDampingCoefficientField;
    struct rollingDampingCoefficientField;
    struct torsionalDampingCoefficientField;
    struct linearStiffnessEnabledField;
    struct dampingEnabledField;

private:
    /** Rejects invalid input before a setter changes values or switches. */
    static void validateNonNegative(Real value)
    {
        if (!math::isFinite(value) || value < 0.0)
        {
            throw std::invalid_argument("Material-pair coefficients must be finite and non-negative.");
        }
    }

    int masterMaterialIndex_{-1};            ///< Canonical smaller material index.
    int slaveMaterialIndex_{-1};             ///< Canonical larger material index.
    Real normalStiffnessPerUnitArea_{-1.0};  ///< LS normal stiffness per nodal area, N/m^3; -1 means unassigned.
    Real slidingStiffnessPerUnitArea_{-1.0}; ///< LS sliding stiffness per nodal area, N/m^3; -1 means unassigned.
    Real normalStiffness_{0.0};              ///< Linear normal stiffness, N/m.
    Real slidingStiffness_{0.0};             ///< Linear sliding stiffness, N/m.
    Real rollingStiffness_{0.0};             ///< Length-scaled rolling stiffness, N/m.
    Real torsionalStiffness_{0.0};           ///< Length-scaled torsional stiffness, N/m.
    Real slidingFrictionCoefficient_{0.0};   ///< Dimensionless sliding friction.
    Real rollingFrictionCoefficient_{0.0};   ///< Dimensionless rolling friction.
    Real torsionalFrictionCoefficient_{0.0}; ///< Dimensionless torsional friction.
    Real normalDampingCoefficient_{0.0};     ///< Normal dashpot coefficient, N s/m.
    Real slidingDampingCoefficient_{0.0};    ///< Sliding dashpot coefficient, N s/m.
    Real rollingDampingCoefficient_{0.0};    ///< Length-scaled rolling dashpot coefficient, N s/m.
    Real torsionalDampingCoefficient_{0.0};  ///< Length-scaled torsional dashpot coefficient, N s/m.
    bool linearStiffnessEnabled_{false};     ///< One switch for the complete four-coefficient linear model.
    bool dampingEnabled_{false};             ///< One switch for all four explicit damping coefficients.

public:
    struct masterMaterialIndexField {
        inline static constexpr auto member = &materialPair::masterMaterialIndex_;
    };
    struct slaveMaterialIndexField {
        inline static constexpr auto member = &materialPair::slaveMaterialIndex_;
    };
    struct normalStiffnessPerUnitAreaField {
        inline static constexpr auto member = &materialPair::normalStiffnessPerUnitArea_;
    };
    struct slidingStiffnessPerUnitAreaField {
        inline static constexpr auto member = &materialPair::slidingStiffnessPerUnitArea_;
    };
    struct normalStiffnessField {
        inline static constexpr auto member = &materialPair::normalStiffness_;
    };
    struct slidingStiffnessField {
        inline static constexpr auto member = &materialPair::slidingStiffness_;
    };
    struct rollingStiffnessField {
        inline static constexpr auto member = &materialPair::rollingStiffness_;
    };
    struct torsionalStiffnessField {
        inline static constexpr auto member = &materialPair::torsionalStiffness_;
    };
    struct slidingFrictionCoefficientField {
        inline static constexpr auto member = &materialPair::slidingFrictionCoefficient_;
    };
    struct rollingFrictionCoefficientField {
        inline static constexpr auto member = &materialPair::rollingFrictionCoefficient_;
    };
    struct torsionalFrictionCoefficientField {
        inline static constexpr auto member = &materialPair::torsionalFrictionCoefficient_;
    };
    struct normalDampingCoefficientField {
        inline static constexpr auto member = &materialPair::normalDampingCoefficient_;
    };
    struct slidingDampingCoefficientField {
        inline static constexpr auto member = &materialPair::slidingDampingCoefficient_;
    };
    struct rollingDampingCoefficientField {
        inline static constexpr auto member = &materialPair::rollingDampingCoefficient_;
    };
    struct torsionalDampingCoefficientField {
        inline static constexpr auto member = &materialPair::torsionalDampingCoefficient_;
    };
    struct linearStiffnessEnabledField {
        inline static constexpr auto member = &materialPair::linearStiffnessEnabled_;
    };
    struct dampingEnabledField {
        inline static constexpr auto member = &materialPair::dampingEnabled_;
    };
    using DeviceLayout = deviceLayout<masterMaterialIndexField,
                                      slaveMaterialIndexField,
                                      normalStiffnessPerUnitAreaField,
                                      slidingStiffnessPerUnitAreaField,
                                      normalStiffnessField,
                                      slidingStiffnessField,
                                      rollingStiffnessField,
                                      torsionalStiffnessField,
                                      slidingFrictionCoefficientField,
                                      rollingFrictionCoefficientField,
                                      torsionalFrictionCoefficientField,
                                      normalDampingCoefficientField,
                                      slidingDampingCoefficientField,
                                      rollingDampingCoefficientField,
                                      torsionalDampingCoefficientField,
                                      linearStiffnessEnabledField,
                                      dampingEnabledField>;
};

using materialPairStorage = hostAoSDeviceSoA<materialPair, materialPair::DeviceLayout>;

/** Append-stable triangular table of unordered material pairs, including its diagonal. */
class materialPairContainer : public materialPairStorage
{
public:
    explicit materialPairContainer(int materialCount = 0) { resizeForMaterials(materialCount); }

    /** Transfers the table and leaves the source empty and reusable. */
    materialPairContainer(materialPairContainer&& other) noexcept : materialPairStorage(std::move(other)), materialCount_(std::exchange(other.materialCount_, 0)) {}

    /** Transfers table storage while preserving the moved-from object's empty-table invariant. */
    materialPairContainer& operator=(materialPairContainer&& other) noexcept
    {
        if (this != &other)
        {
            materialPairStorage::operator=(std::move(other));
            materialCount_ = std::exchange(other.materialCount_, 0);
        }
        return *this;
    }

    /** Extends the table for appended materials, retaining all existing overrides and pair indices. */
    void resizeForMaterials(int count)
    {
        if (count < materialCount_)
        {
            throw std::invalid_argument("Material-pair storage cannot remove existing materials.");
        }
        const auto total = static_cast<unsigned long long>(count) * (static_cast<unsigned long long>(count) + 1) / 2;
        if (total > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
        {
            throw std::length_error("The material-pair table exceeds the supported index range.");
        }
        auto& values = host();
        values.reserve(static_cast<std::size_t>(total));
        for (int larger = materialCount_; larger < count; ++larger)
        {
            for (int smaller = 0; smaller <= larger; ++smaller)
            {
                values.emplace_back(smaller, larger);
            }
        }
        materialCount_ = count;
    }

    materialPair& at(int first, int second)
    {
        checkIndices(first, second);
        return host()[index(first, second)];
    }
    const materialPair& at(int first, int second) const
    {
        checkIndices(first, second);
        return host()[index(first, second)];
    }
    int materialCount() const noexcept { return materialCount_; }

    /** Maps valid material IDs to a symmetric slot independent of the current material count. */
    FUNDEM_MATH_HD static int index(int first, int second) noexcept
    {
        const int larger = first >= second ? first : second;
        const int smaller = first >= second ? second : first;
        return static_cast<int>(static_cast<long long>(larger) * (larger + 1LL) / 2 + smaller);
    }

private:
    /** Checks configuration-time lookups; device kernels use the unchecked static index. */
    void checkIndices(int first, int second) const
    {
        if (first < 0 || second < 0 || first >= materialCount_ || second >= materialCount_)
        {
            throw std::out_of_range("Material-pair indices are outside this table.");
        }
    }
    int materialCount_{0}; ///< Number of materials represented by the table.
};

} // namespace fundem
