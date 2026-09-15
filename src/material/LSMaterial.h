/**
 * @file LSMaterial.h
 * @brief Defines the level-set-specific interface over unified material storage.
 */
#pragma once

#include "material.h"

namespace fundem
{

/** Level-set contact material using the storage defined by material. */
class LSMaterial : public material
{
public:
    LSMaterial() noexcept : material(materialType::levelSet) {}

    /**
     * Constructs a complete level-set material whose stiffness values are per
     * unit contact area.
     */
    LSMaterial(Real normalStiffnessPerUnitArea, Real shearStiffnessPerUnitArea, Real frictionCoefficient, Real restitutionCoefficient, Real density) : material(materialType::levelSet)
    {
        setProperties(normalStiffnessPerUnitArea, shearStiffnessPerUnitArea, frictionCoefficient, restitutionCoefficient, density);
    }

    Real normalStiffnessPerUnitArea() const noexcept { return material::normalStiffness(); }
    Real shearStiffnessPerUnitArea() const noexcept { return material::slidingStiffness(); }
    Real frictionCoefficient() const noexcept { return material::slidingFrictionCoefficient(); }

    using material::rollingFrictionCoefficient;
    using material::setRollingFrictionCoefficient;
    using material::setSlidingFrictionCoefficient;
    using material::setTorsionalFrictionCoefficient;
    using material::slidingFrictionCoefficient;
    using material::torsionalFrictionCoefficient;

    void setNormalStiffnessPerUnitArea(Real value) { material::setNormalStiffness(value); }
    void setShearStiffnessPerUnitArea(Real value) { material::setSlidingStiffness(value); }
    void setFrictionCoefficient(Real value) { material::setSlidingFrictionCoefficient(value); }

    /**
     * Replaces the validated level-set stiffness, friction, restitution, and
     * density properties. Rolling and torsional stiffness are kept at zero.
     */
    void setProperties(Real normalStiffnessPerUnitArea, Real shearStiffnessPerUnitArea, Real frictionCoefficient, Real restitutionCoefficient, Real density)
    {
        material::setProperties(normalStiffnessPerUnitArea, shearStiffnessPerUnitArea, 0.0, 0.0, frictionCoefficient, 0.0, 0.0, restitutionCoefficient, density);
    }

private:
    using material::normalStiffness;
    using material::rollingStiffness;
    using material::setNormalStiffness;
    using material::setProperties;
    using material::setRollingStiffness;
    using material::setSlidingStiffness;
    using material::setTorsionalStiffness;
    using material::slidingStiffness;
    using material::torsionalStiffness;
};

static_assert(sizeof(LSMaterial) == sizeof(material), "LSMaterial must not add storage to material.");

} // namespace fundem
