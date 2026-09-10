/**
 * @file contactFunctions.h
 * @brief Implements host/device contact forces, damping, friction, and elastic energy.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem::execution
{

using math::Real;
using math::Vec3;

/** Combines two positive coefficients acting in series. */
FUNDEM_MATH_HD inline Real seriesCombination(Real first, Real second) noexcept
{
    if (first <= 0.0 || second <= 0.0)
    {
        return 0.0;
    }
    return first <= second ? first / (1.0 + first / second) : second / (1.0 + second / first);
}

/** Combines pair stiffness while preserving the deformable side of a fixed/deformable pair. */
FUNDEM_MATH_HD inline Real effectiveStiffness(Real first, Real second, bool firstInfiniteMass, bool secondInfiniteMass) noexcept
{
    if (firstInfiniteMass && secondInfiniteMass)
    {
        return 0.0;
    }
    if (firstInfiniteMass)
    {
        return second;
    }
    if (secondInfiniteMass)
    {
        return first;
    }
    return seriesCombination(first, second);
}

/** Returns the harmonic mean of two positive values. */
FUNDEM_MATH_HD inline Real harmonicMean(Real first, Real second) noexcept
{
    if (first <= 0.0 || second <= 0.0)
    {
        return 0.0;
    }
    return first <= second ? first / (0.5 + 0.5 * first / second) : second / (0.5 + 0.5 * second / first);
}

/** Returns scalar linear-spring energy. */
FUNDEM_MATH_HD inline Real springElasticEnergy(Real stiffness, Real deformation) noexcept { return stiffness > 0.0 ? 0.5 * stiffness * deformation * deformation : 0.0; }

/** Returns vector linear-spring energy. */
FUNDEM_MATH_HD inline Real springElasticEnergy(Real stiffness, const Vec3& deformation) noexcept { return stiffness > 0.0 ? 0.5 * stiffness * math::normSquared(deformation) : 0.0; }

/** Distributes rigid-pair damping mass to one active level-set surface node by area. */
FUNDEM_MATH_HD inline Real nodalContactDampingMass(Real rigidPairEffectiveMass, Real nodeArea, Real activePatchArea) noexcept
{
    if (rigidPairEffectiveMass <= 0.0 || nodeArea <= 0.0 || activePatchArea <= 0.0)
    {
        return 0.0;
    }
    const Real areaFraction = nodeArea < activePatchArea ? nodeArea / activePatchArea : 1.0;
    return rigidPairEffectiveMass * areaFraction;
}

/** Returns rigid-body velocity at a world-space contact point. */
FUNDEM_MATH_HD inline Vec3 velocityAtContactPoint(const Vec3& position, const Vec3& velocity, const Vec3& angularVelocity, const Vec3& contactPoint) noexcept
{
    return velocity + math::cross(angularVelocity, contactPoint - position);
}

/** Returns first-body minus second-body velocity at a shared contact point. */
FUNDEM_MATH_HD inline Vec3 relativeVelocityAtContactPoint(const Vec3& firstPosition,
                                                          const Vec3& firstVelocity,
                                                          const Vec3& firstAngularVelocity,
                                                          const Vec3& secondPosition,
                                                          const Vec3& secondVelocity,
                                                          const Vec3& secondAngularVelocity,
                                                          const Vec3& contactPoint) noexcept
{
    return velocityAtContactPoint(firstPosition, firstVelocity, firstAngularVelocity, contactPoint) - velocityAtContactPoint(secondPosition, secondVelocity, secondAngularVelocity, contactPoint);
}

/** Converts restitution coefficient to the nondimensional linear dashpot factor. */
FUNDEM_MATH_HD inline Real contactDissipationFactor(Real restitutionCoefficient) noexcept
{
    if (restitutionCoefficient <= 0.0)
    {
        return 1.0;
    }
    const Real logarithm = math::detail::log(restitutionCoefficient);
    return -logarithm / math::detail::sqrt(logarithm * logarithm + math::pi * math::pi);
}

/** Returns the critical-damping-scaled dashpot coefficient. */
FUNDEM_MATH_HD inline Real contactDampingCoefficient(Real dissipation, Real dampingMass, Real stiffness) noexcept
{
    return dissipation > 0.0 && dampingMass > 0.0 && stiffness > 0.0 ? 2.0 * dissipation * math::detail::sqrt(dampingMass * stiffness) : 0.0;
}

/**
 * Scales only the damping force so the total tangential force stays inside the Coulomb disk.
 *
 * Let F(alpha) = Fe + alpha Fd, 0 <= alpha <= 1, and L = mu |Fn|. If the full trial force
 * lies outside ||F|| <= L, the boundary intersection follows from
 *
 *     a alpha^2 + 2 b alpha + c = 0,
 *     a = Fd.Fd,  b = Fe.Fd,  c = Fe.Fe - L^2.
 *
 * Fe has already undergone radial return, so c <= 0. The positive root is evaluated with
 * an anti-cancellation branch. This function never changes the elastic spring history.
 */
FUNDEM_MATH_HD inline Vec3 limitDampingWithinCoulombDisk(const Vec3& elasticForce, const Vec3& dampingForce, Real frictionLimit) noexcept
{
    if (frictionLimit <= 0.0)
    {
        return Vec3::zero();
    }

    const Real limitSquared = frictionLimit * frictionLimit;
    const Vec3 fullTrialForce = elasticForce + dampingForce;
    if (math::normSquared(fullTrialForce) <= limitSquared)
    {
        return fullTrialForce;
    }

    const Real a = math::normSquared(dampingForce);
    if (a == 0.0)
    {
        const Real elasticForceLength = math::norm(elasticForce);
        return elasticForceLength > frictionLimit ? (frictionLimit / elasticForceLength) * elasticForce : elasticForce;
    }

    const Real b = math::dot(elasticForce, dampingForce);
    const Real c = math::normSquared(elasticForce) - limitSquared;
    const Real discriminant = b * b - a * c;
    const Real squareRoot = math::detail::sqrt(discriminant > 0.0 ? discriminant : 0.0);

    Real positiveRoot = (-b + squareRoot) / a;
    if (b >= 0.0)
    {
        const Real denominator = b + squareRoot;
        positiveRoot = denominator > 0.0 ? -c / denominator : 0.0;
    }

    Real dampingScale = positiveRoot;
    if (dampingScale < 0.0)
    {
        dampingScale = 0.0;
    }
    else if (dampingScale > 1.0)
    {
        dampingScale = 1.0;
    }
    return elasticForce + dampingScale * dampingForce;
}

/** Integrates a tangential spring, performs Coulomb radial return, and adds admissible damping. */
FUNDEM_MATH_HD inline Vec3 integrateTangentialSpringForce(Vec3& springDeformation,
                                                          const Vec3& springVelocity,
                                                          const Vec3& contactNormal,
                                                          Real absoluteNormalForce,
                                                          Real frictionCoefficient,
                                                          Real stiffness,
                                                          Real dampingCoefficient,
                                                          Real timeStep) noexcept
{
    const Real frictionLimit = frictionCoefficient * absoluteNormalForce;
    if (frictionLimit <= 0.0 || stiffness <= 0.0)
    {
        springDeformation = Vec3::zero();
        return Vec3::zero();
    }

    springDeformation -= math::dot(springDeformation, contactNormal) * contactNormal;
    springDeformation += timeStep * springVelocity;
    Vec3 elasticForce = -stiffness * springDeformation;
    const Real elasticForceLength = math::norm(elasticForce);
    if (elasticForceLength > frictionLimit)
    {
        elasticForce *= frictionLimit / elasticForceLength;
        springDeformation = -elasticForce / stiffness;
    }

    return limitDampingWithinCoulombDisk(elasticForce, -dampingCoefficient * springVelocity, frictionLimit);
}

/** Integrates the normal-axis torsional spring subject to its friction limit. */
FUNDEM_MATH_HD inline Vec3 integrateTorsionalSpringForce(Vec3& springDeformation,
                                                         const Vec3& springVelocity,
                                                         const Vec3& contactNormal,
                                                         Real absoluteNormalForce,
                                                         Real frictionCoefficient,
                                                         Real stiffness,
                                                         Real dampingCoefficient,
                                                         Real timeStep) noexcept
{
    const Real frictionLimit = frictionCoefficient * absoluteNormalForce;
    if (frictionLimit <= 0.0 || stiffness <= 0.0)
    {
        springDeformation = Vec3::zero();
        return Vec3::zero();
    }

    springDeformation = math::dot(springDeformation + timeStep * springVelocity, contactNormal) * contactNormal;
    Vec3 elasticForce = -stiffness * springDeformation;
    const Real elasticForceLength = math::norm(elasticForce);
    if (elasticForceLength > frictionLimit)
    {
        elasticForce *= frictionLimit / elasticForceLength;
        springDeformation = -elasticForce / stiffness;
    }

    return limitDampingWithinCoulombDisk(elasticForce, -dampingCoefficient * springVelocity, frictionLimit);
}

/** Computes normal and sliding contact force while updating sliding history. */
FUNDEM_MATH_HD inline Vec3 calculateLinearContactForce(Real& normalForceMagnitude,
                                                       Vec3& slidingSpringDeformation,
                                                       const Vec3& relativeVelocity,
                                                       const Vec3& contactNormal,
                                                       Real overlap,
                                                       Real timeStep,
                                                       Real normalStiffness,
                                                       Real slidingStiffness,
                                                       Real slidingFrictionCoefficient,
                                                       Real restitutionCoefficient,
                                                       Real dampingMass) noexcept
{
    normalForceMagnitude = 0.0;
    if (overlap <= 0.0)
    {
        slidingSpringDeformation = Vec3::zero();
        return Vec3::zero();
    }

    const Real dissipation = contactDissipationFactor(restitutionCoefficient);
    const Real normalDampingCoefficient = contactDampingCoefficient(dissipation, dampingMass, normalStiffness);
    const Real normalVelocity = math::dot(relativeVelocity, contactNormal);
    const Vec3 normalRelativeVelocity = normalVelocity * contactNormal;
    const Real trialNormalForceMagnitude = normalStiffness * overlap - normalDampingCoefficient * normalVelocity;
    normalForceMagnitude = trialNormalForceMagnitude > 0.0 ? trialNormalForceMagnitude : 0.0;

    const Real slidingDampingCoefficient = contactDampingCoefficient(dissipation, dampingMass, slidingStiffness);
    const Vec3 slidingForce = integrateTangentialSpringForce(slidingSpringDeformation,
                                                             relativeVelocity - normalRelativeVelocity,
                                                             contactNormal,
                                                             normalForceMagnitude,
                                                             slidingFrictionCoefficient,
                                                             slidingStiffness,
                                                             slidingDampingCoefficient,
                                                             timeStep);
    return normalForceMagnitude * contactNormal + slidingForce;
}

/** Computes rolling and torsional resistance torque while updating both histories. */
FUNDEM_MATH_HD inline Vec3 calculateSphereRotationalContactTorque(Vec3& rollingSpringDeformation,
                                                                  Vec3& torsionalSpringDeformation,
                                                                  const Vec3& relativeAngularVelocity,
                                                                  const Vec3& contactNormal,
                                                                  Real timeStep,
                                                                  Real rollingStiffness,
                                                                  Real torsionalStiffness,
                                                                  Real rollingFrictionCoefficient,
                                                                  Real torsionalFrictionCoefficient,
                                                                  Real restitutionCoefficient,
                                                                  Real effectiveMass,
                                                                  Real effectiveRadius,
                                                                  Real absoluteNormalForce) noexcept
{
    const Real dissipation = contactDissipationFactor(restitutionCoefficient);
    const Real rollingDampingCoefficient = contactDampingCoefficient(dissipation, effectiveMass, rollingStiffness);
    const Vec3 rollingRelativeVelocity = -effectiveRadius * math::cross(contactNormal, relativeAngularVelocity);
    const Vec3 rollingForce = integrateTangentialSpringForce(rollingSpringDeformation,
                                                             rollingRelativeVelocity,
                                                             contactNormal,
                                                             absoluteNormalForce,
                                                             rollingFrictionCoefficient,
                                                             rollingStiffness,
                                                             rollingDampingCoefficient,
                                                             timeStep);

    const Real torsionalDampingCoefficient = contactDampingCoefficient(dissipation, effectiveMass, torsionalStiffness);
    const Real effectiveDiameter = 2.0 * effectiveRadius;
    const Vec3 torsionalRelativeVelocity = effectiveDiameter * math::dot(contactNormal, relativeAngularVelocity) * contactNormal;
    const Vec3 torsionalForce = integrateTorsionalSpringForce(torsionalSpringDeformation,
                                                              torsionalRelativeVelocity,
                                                              contactNormal,
                                                              absoluteNormalForce,
                                                              torsionalFrictionCoefficient,
                                                              torsionalStiffness,
                                                              torsionalDampingCoefficient,
                                                              timeStep);
    return effectiveRadius * math::cross(contactNormal, rollingForce) + effectiveDiameter * torsionalForce;
}

} // namespace fundem::execution
