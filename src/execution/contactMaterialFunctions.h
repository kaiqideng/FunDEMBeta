/**
 * @file contactMaterialFunctions.h
 * @brief Resolves material-pair contact response identically on host and device.
 */
#pragma once

#include "execution/contactFunctions.h"
#include "material/materialPair.h"

#include <type_traits>

namespace fundem
{
class LSParticle;
}

namespace fundem::execution
{

/** Geometry family selected by the particle containers, never by a material flag. */
enum class contactType : int
{
    sphereSphere,
    sphereLevelSet,
    levelSetLevelSet
};

/** Selects the shared response policy before particle-container types are erased. */
template <class MasterParticle, class SlaveParticle> constexpr contactType contactTypeFor() noexcept
{
    if constexpr (std::is_base_of_v<LSParticle, MasterParticle> && std::is_base_of_v<LSParticle, SlaveParticle>)
        return contactType::levelSetLevelSet;
    if constexpr (std::is_base_of_v<LSParticle, MasterParticle> || std::is_base_of_v<LSParticle, SlaveParticle>)
        return contactType::sphereLevelSet;
    return contactType::sphereSphere;
}

/** Per-contact stiffness, physical dashpot coefficients and conservative normal energy. */
struct contactResponse {
    Real normalStiffness_{0.0};             ///< Normal secant stiffness, force divided by overlap.
    Real slidingStiffness_{0.0};            ///< Current tangential stiffness.
    Real rollingStiffness_{0.0};            ///< Length-scaled rolling spring stiffness.
    Real torsionalStiffness_{0.0};          ///< Length-scaled torsional spring stiffness.
    Real normalDampingCoefficient_{0.0};    ///< Actual normal dashpot coefficient.
    Real slidingDampingCoefficient_{0.0};   ///< Actual tangential dashpot coefficient.
    Real rollingDampingCoefficient_{0.0};   ///< Actual rolling pseudo-force dashpot coefficient.
    Real torsionalDampingCoefficient_{0.0}; ///< Actual torsional pseudo-force dashpot coefficient.
    Real normalElasticEnergy_{0.0};         ///< Integral of the normal elastic force over penetration.
};

/**
 * Recovers A_node / A_patch from the mass weighting already established by LS detection.
 * Since m_node = m_pair A_node/A_patch and 1/m_pair = invMass1+invMass2,
 * no second area reduction or additional contact/device field is needed.
 */
FUNDEM_MATH_HD inline Real nodalDampingWeight(Real nodalEffectiveMass, Real firstInverseMass, Real secondInverseMass) noexcept
{
    const Real weight = nodalEffectiveMass * (firstInverseMass + secondInverseMass);
    return weight < 0.0 ? 0.0 : weight > 1.0 ? 1.0 : weight;
}

/**
 * Resolves final pair values, Hertz-Mindlin defaults and restitution-derived damping.
 *
 * Hertz uses a=sqrt(R*delta), Fn=(4/3)E*a*delta and U=(2/5)Fn*delta;
 * its normal tangent stiffness is (3/2)Fn/delta, not the secant stiffness.
 * Mindlin uses Ks=8G*a; default length-scaled rolling/torsional stiffnesses
 * equal Ks. Hertz/Mindlin translational restitution damping is scaled by sqrt(5/6).
 * One linear switch selects all four sphere-related stiffnesses; LS per-area
 * stiffness does not use that switch. One damping switch selects all four
 * explicit dashpots; LS explicit values are distributed with dampingWeight.
 * Linear modes retain c=2*beta*sqrt(m*k).
 * Material-pair values are final: no mass-based stiffness or friction mixing occurs.
 * LS stiffness availability is checked once by solver initialization.
 */
FUNDEM_MATH_HD inline contactResponse resolveContactResponse(const materialPairCoefficients& pair,
                                                             contactType type,
                                                             Real masterYoungModulus,
                                                             Real masterPoissonRatio,
                                                             Real masterRestitution,
                                                             Real slaveYoungModulus,
                                                             Real slavePoissonRatio,
                                                             Real slaveRestitution,
                                                             Real effectiveMass,
                                                             Real effectiveRadius,
                                                             Real overlap,
                                                             Real area,
                                                             Real dampingWeight) noexcept
{
    contactResponse response;
    const bool nodal = type == contactType::levelSetLevelSet;
    if (nodal && (pair.normalStiffnessPerUnitArea_ < 0.0 || pair.slidingStiffnessPerUnitArea_ < 0.0))
        return response;
    const bool linear = nodal || pair.linearStiffnessEnabled_;
    if (nodal)
    {
        response.normalStiffness_ = pair.normalStiffnessPerUnitArea_ * area;
        response.slidingStiffness_ = pair.slidingStiffnessPerUnitArea_ * area;
    }
    else if (linear)
    {
        response.normalStiffness_ = pair.normalStiffness_;
        response.slidingStiffness_ = pair.slidingStiffness_;
        response.rollingStiffness_ = pair.rollingStiffness_;
        response.torsionalStiffness_ = pair.torsionalStiffness_;
    }
    else if (overlap > 0.0 && effectiveRadius > 0.0)
    {
        const Real contactRadius = math::detail::sqrt(effectiveRadius * overlap);
        const Real effectiveYoungModulus = seriesCombination(masterYoungModulus / (1.0 - masterPoissonRatio * masterPoissonRatio), slaveYoungModulus / (1.0 - slavePoissonRatio * slavePoissonRatio));
        const Real masterShearModulus = masterYoungModulus / (2.0 * (1.0 + masterPoissonRatio));
        const Real slaveShearModulus = slaveYoungModulus / (2.0 * (1.0 + slavePoissonRatio));
        const Real effectiveShearModulus = seriesCombination(masterShearModulus / (2.0 - masterPoissonRatio), slaveShearModulus / (2.0 - slavePoissonRatio));
        response.normalStiffness_ = (4.0 / 3.0) * effectiveYoungModulus * contactRadius;
        response.slidingStiffness_ = 8.0 * effectiveShearModulus * contactRadius;
        response.rollingStiffness_ = response.slidingStiffness_;
        response.torsionalStiffness_ = response.slidingStiffness_;
    }

    if (pair.dampingEnabled_)
    {
        const Real weight = nodal ? dampingWeight : 1.0;
        response.normalDampingCoefficient_ = pair.normalDampingCoefficient_ * weight;
        response.slidingDampingCoefficient_ = pair.slidingDampingCoefficient_ * weight;
        if (!nodal)
        {
            response.rollingDampingCoefficient_ = pair.rollingDampingCoefficient_;
            response.torsionalDampingCoefficient_ = pair.torsionalDampingCoefficient_;
        }
    }
    else
    {
        const Real dissipation = contactDissipationFactor(harmonicMean(masterRestitution, slaveRestitution));
        const Real scale = linear ? 1.0 : math::detail::sqrt(5.0 / 6.0);
        response.normalDampingCoefficient_ = scale * contactDampingCoefficient(dissipation, effectiveMass, (linear ? 1.0 : 1.5) * response.normalStiffness_);
        response.slidingDampingCoefficient_ = scale * contactDampingCoefficient(dissipation, effectiveMass, response.slidingStiffness_);
        if (!nodal)
        {
            response.rollingDampingCoefficient_ = contactDampingCoefficient(dissipation, effectiveMass, response.rollingStiffness_);
            response.torsionalDampingCoefficient_ = contactDampingCoefficient(dissipation, effectiveMass, response.torsionalStiffness_);
        }
    }
    response.normalElasticEnergy_ = overlap > 0.0 ? (linear ? 0.5 : 0.4) * response.normalStiffness_ * overlap * overlap : 0.0;
    return response;
}

} // namespace fundem::execution
