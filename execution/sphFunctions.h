/**
 * @file sphFunctions.h
 * @brief Provides shared WCSPH kernels, Riemann states, and fluid-wall formulas.
 */
#pragma once

#include "math/Quaternion.h"

#include <cmath>

namespace fundem::execution
{

using math::Mat3;
using math::Quaternion;
using math::Real;
using math::Vec3;

inline constexpr Real minimumSPHSmoothingLengthToSpacingRatio = 1.0; ///< Smallest accepted smoothing-length/spacing ratio.
inline constexpr Real maximumSPHSmoothingLengthToSpacingRatio = 2.0; ///< Largest accepted smoothing-length/spacing ratio.
inline constexpr Real minimumSPHSoundSpeedToVelocityRatio = 10.0;    ///< Minimum weak-compressibility sound-speed ratio.
inline constexpr Real minimumSPHDensityRatio = 0.8;                  ///< Lower physical density-ratio diagnostic bound.
inline constexpr Real maximumSPHDensityRatio = 1.2;                  ///< Upper physical density-ratio diagnostic bound.
inline constexpr Real numericalSPHDensityFloorRatio = 0.1;           ///< Emergency density floor used by integration.
inline constexpr Real numericalSPHDensityCeilingRatio = 10.0;        ///< Emergency density ceiling used by integration.

/** Regularizes density denominators with the shared numerical tolerance. */
FUNDEM_MATH_HD inline Real nonZeroDensity(Real density) noexcept { return density > math::defaultTolerance ? density : math::defaultTolerance; }

/** Returns the three-dimensional Wendland C2 normalization coefficient. */
FUNDEM_MATH_HD inline Real wendlandKernelNormalization3D(Real smoothingLength) noexcept
{
    return smoothingLength > 0.0 ? 21.0 / (16.0 * math::pi * smoothingLength * smoothingLength * smoothingLength) : 0.0;
}

/** Evaluates the compact three-dimensional Wendland C2 kernel. */
FUNDEM_MATH_HD inline Real wendlandKernel3D(Real distance, Real smoothingLength) noexcept
{
    if (distance < 0.0 || smoothingLength <= 0.0)
    {
        return 0.0;
    }
    const Real normalizedDistance = distance / smoothingLength;
    if (normalizedDistance >= 2.0)
    {
        return 0.0;
    }
    const Real factor = 1.0 - 0.5 * normalizedDistance;
    return wendlandKernelNormalization3D(smoothingLength) * factor * factor * factor * factor * (2.0 * normalizedDistance + 1.0);
}

/** Evaluates the particle-position gradient of the Wendland C2 kernel. */
FUNDEM_MATH_HD inline Vec3 wendlandKernelGradient3D(const Vec3& separation, Real smoothingLength) noexcept
{
    const Real distanceSquared = math::normSquared(separation);
    if (distanceSquared <= math::defaultToleranceSquared || smoothingLength <= 0.0)
    {
        return Vec3::zero();
    }
    const Real distance = math::detail::sqrt(distanceSquared);
    const Real normalizedDistance = distance / smoothingLength;
    if (normalizedDistance >= 2.0)
    {
        return Vec3::zero();
    }
    const Real factor = 1.0 - 0.5 * normalizedDistance;
    const Real gradientScale = -5.0 * wendlandKernelNormalization3D(smoothingLength) * factor * factor * factor / (smoothingLength * smoothingLength);
    return gradientScale * separation;
}

/** Sums the kernel over an infinite-lattice stencil truncated by compact support. */
inline Real calculateLatticeKernelSum3D(Real particleSpacing, Real smoothingLength) noexcept
{
    if (particleSpacing <= 0.0 || smoothingLength <= 0.0)
    {
        return 0.0;
    }
    const Real cutoffRadius = 2.0 * smoothingLength;
    const int searchDepth = static_cast<int>(std::ceil(cutoffRadius / particleSpacing));
    Real sum = 0.0;
    for (int z = -searchDepth; z <= searchDepth; ++z)
    {
        for (int y = -searchDepth; y <= searchDepth; ++y)
        {
            for (int x = -searchDepth; x <= searchDepth; ++x)
            {
                const Vec3 separation = particleSpacing * Vec3{Real(x), Real(y), Real(z)};
                if (math::normSquared(separation) < cutoffRadius * cutoffRadius)
                {
                    sum += wendlandKernel3D(math::norm(separation), smoothingLength);
                }
            }
        }
    }
    return sum;
}

/** Solves the acoustic Riemann interface normal velocity. */
FUNDEM_MATH_HD inline Real computeRiemannInterfaceNormalVelocity(Real leftPressure,
                                                                 Real rightPressure,
                                                                 Real leftNormalVelocity,
                                                                 Real rightNormalVelocity,
                                                                 Real referenceDensity,
                                                                 Real soundSpeed) noexcept
{
    const Real acousticImpedance = nonZeroDensity(referenceDensity) * soundSpeed;
    return 0.5 * (leftNormalVelocity + rightNormalVelocity) + (leftPressure - rightPressure) / (2.0 * acousticImpedance);
}

/** Solves the dissipative Riemann interface pressure with limited compression speed. */
FUNDEM_MATH_HD inline Real computeRiemannInterfacePressure(Real leftPressure, Real rightPressure, Real leftNormalVelocity, Real rightNormalVelocity, Real referenceDensity, Real soundSpeed) noexcept
{
    const Real velocityDifference = leftNormalVelocity - rightNormalVelocity;
    const Real compressionSpeed = velocityDifference > 0.0 ? velocityDifference : 0.0;
    const Real limitedDissipationSpeed = 3.0 * compressionSpeed < soundSpeed ? 3.0 * compressionSpeed : soundSpeed;
    return 0.5 * (leftPressure + rightPressure) + 0.5 * referenceDensity * limitedDissipationSpeed * velocityDifference;
}

/** Reconstructs the full interface velocity from its normal solution and mean tangent. */
FUNDEM_MATH_HD inline Vec3 computeRiemannInterfaceVelocity(const Vec3& leftVelocity,
                                                           const Vec3& rightVelocity,
                                                           const Vec3& riemannPlaneNormal,
                                                           Real leftNormalVelocity,
                                                           Real rightNormalVelocity,
                                                           Real interfaceNormalVelocity) noexcept
{
    const Vec3 meanVelocity = 0.5 * (leftVelocity + rightVelocity);
    return interfaceNormalVelocity * riemannPlaneNormal + meanVelocity - 0.5 * (leftNormalVelocity + rightNormalVelocity) * riemannPlaneNormal;
}

/** Returns one neighbor's Riemann continuity-equation contribution. */
FUNDEM_MATH_HD inline Real computeDensityRateContribution(Real density, Real neighborVolume, const Vec3& velocity, const Vec3& interfaceVelocity, const Vec3& kernelGradient) noexcept
{
    return 2.0 * density * neighborVolume * math::dot(velocity - interfaceVelocity, kernelGradient);
}

/** Returns one symmetric pressure-acceleration contribution. */
FUNDEM_MATH_HD inline Vec3 computePressureAcceleration(Real neighborVolume, Real interfacePressure, Real density, const Vec3& kernelGradient) noexcept
{
    return -(2.0 * neighborVolume * interfacePressure / density) * kernelGradient;
}

/** Returns one pairwise dynamic-viscosity acceleration contribution. */
FUNDEM_MATH_HD inline Vec3 computeViscousAcceleration(const Vec3& velocity,
                                                      const Vec3& neighborVelocity,
                                                      Real dynamicViscosity,
                                                      Real neighborVolume,
                                                      Real density,
                                                      Real distance,
                                                      Real radialKernelDerivative) noexcept
{
    return (2.0 * dynamicViscosity * neighborVolume * radialKernelDerivative / (density * distance)) * (velocity - neighborVelocity);
}

/** Evaluates the linear weakly-compressible equation of state. */
FUNDEM_MATH_HD inline Real computePressureFromDensity(Real density, Real referenceDensity, Real soundSpeed) noexcept { return soundSpeed * soundSpeed * (density - referenceDensity); }

/** Inverts the linear weakly-compressible equation of state. */
FUNDEM_MATH_HD inline Real computeDensityFromPressure(Real pressure, Real referenceDensity, Real soundSpeed) noexcept { return referenceDensity + pressure / (soundSpeed * soundSpeed); }

/** Returns the acoustic CFL limit used for pressure-relaxation substeps. */
FUNDEM_MATH_HD inline Real calculateSPHAcousticTimeStep(Real smoothingLength, Real soundSpeed, Real maximumFluidVelocity) noexcept
{
    const Real denominator = soundSpeed + maximumFluidVelocity;
    return smoothingLength > 0.0 && soundSpeed > 0.0 && maximumFluidVelocity >= 0.0 ? 0.6 * smoothingLength / denominator : 0.0;
}

/** Returns the advection limit from pre-reduced global smoothing and viscosity scales. */
FUNDEM_MATH_HD inline Real calculateSPHAdvectionTimeStepFromLimits(Real minimumSmoothingLength,
                                                                   Real minimumViscousTimeScale,
                                                                   Real maximumFluidVelocity,
                                                                   Real maximumFluidAcceleration,
                                                                   Real configuredMaximumVelocity) noexcept;

/** Returns the advection limit for one uniform SPH parameter set. */
FUNDEM_MATH_HD inline Real calculateSPHAdvectionTimeStep(Real smoothingLength,
                                                         Real referenceDensity,
                                                         Real dynamicViscosity,
                                                         Real maximumFluidVelocity,
                                                         Real maximumFluidAcceleration,
                                                         Real configuredMaximumVelocity) noexcept
{
    if (referenceDensity <= 0.0 || dynamicViscosity < 0.0)
    {
        return 0.0;
    }
    const Real viscousTimeScale = dynamicViscosity > 0.0 ? referenceDensity * smoothingLength * smoothingLength / dynamicViscosity : 0.0;
    return calculateSPHAdvectionTimeStepFromLimits(smoothingLength, viscousTimeScale, maximumFluidVelocity, maximumFluidAcceleration, configuredMaximumVelocity);
}

/** Implements the velocity, acceleration, and viscous branches of the global advection limit. */
FUNDEM_MATH_HD inline Real calculateSPHAdvectionTimeStepFromLimits(Real minimumSmoothingLength,
                                                                   Real minimumViscousTimeScale,
                                                                   Real maximumFluidVelocity,
                                                                   Real maximumFluidAcceleration,
                                                                   Real configuredMaximumVelocity) noexcept
{
    if (minimumSmoothingLength <= 0.0 || minimumViscousTimeScale < 0.0 || maximumFluidVelocity < 0.0 || maximumFluidAcceleration < 0.0 || configuredMaximumVelocity <= 0.0)
    {
        return 0.0;
    }

    const Real accelerationVelocityScale = math::detail::sqrt(4.0 * minimumSmoothingLength * maximumFluidAcceleration);
    Real velocityScale = maximumFluidVelocity > configuredMaximumVelocity ? maximumFluidVelocity : configuredMaximumVelocity;
    velocityScale = accelerationVelocityScale > velocityScale ? accelerationVelocityScale : velocityScale;
    Real timeStep = minimumSmoothingLength / velocityScale;
    if (minimumViscousTimeScale > 0.0)
    {
        timeStep = minimumViscousTimeScale < timeStep ? minimumViscousTimeScale : timeStep;
    }
    return 0.25 * timeStep;
}

namespace detail
{

/** Precomputed pair geometry shared by pressure, viscosity, and density terms. */
struct SPHInteractionGeometry {
    Vec3 neighborDirection_;   ///< Unit direction from neighbor to particle.
    Vec3 kernelGradient_;      ///< Kernel gradient with respect to particle position.
    Real distance_;            ///< Pair separation.
    Real kernelGradientScale_; ///< Scalar multiplying the separation vector.
    bool valid_;               ///< Whether the pair lies in nondegenerate support.
};

/** Builds reusable compact-support geometry for a particle-neighbor pair. */
FUNDEM_MATH_HD inline SPHInteractionGeometry makeSPHInteractionGeometry(const Vec3& position, const Vec3& neighborPosition, Real smoothingLength) noexcept
{
    const Vec3 separation = position - neighborPosition;
    const Real distanceSquared = math::normSquared(separation);
    if (smoothingLength <= 0.0 || distanceSquared <= math::defaultToleranceSquared || distanceSquared >= 4.0 * smoothingLength * smoothingLength)
    {
        return {Vec3::zero(), Vec3::zero(), 0.0, 0.0, false};
    }

    const Real distance = math::detail::sqrt(distanceSquared);
    const Real normalizedDistance = distance / smoothingLength;
    const Real kernelFactor = 1.0 - 0.5 * normalizedDistance;
    const Real kernelGradientScale = -5.0 * wendlandKernelNormalization3D(smoothingLength) * kernelFactor * kernelFactor * kernelFactor / (smoothingLength * smoothingLength);
    return {separation / distance, kernelGradientScale * separation, distance, kernelGradientScale, true};
}

/** One-sided wall Riemann state and associated pair geometry. */
struct WallRiemannState {
    SPHInteractionGeometry geometry_; ///< Particle-wall pair geometry.
    Vec3 riemannPlaneNormal_;         ///< Normal used by the Riemann problem.
    Real density_;                    ///< Regularized fluid density.
    Real leftNormalVelocity_;         ///< Fluid-side normal velocity.
    Real rightNormalVelocity_;        ///< Reflected wall-side normal velocity.
    Real wallPressure_;               ///< Reconstructed wall pressure.
};

/** Builds the reflected one-sided state used by fluid-wall Riemann interactions. */
FUNDEM_MATH_HD inline WallRiemannState makeWallRiemannState(const Vec3& position,
                                                            const Vec3& velocity,
                                                            Real rawDensity,
                                                            Real pressure,
                                                            const Vec3& wallPosition,
                                                            const Vec3& wallNormal,
                                                            const Vec3& wallVelocity,
                                                            const Vec3& wallAcceleration,
                                                            Real smoothingLength,
                                                            const Vec3& gravity) noexcept
{
    const SPHInteractionGeometry geometry = makeSPHInteractionGeometry(position, wallPosition, smoothingLength);
    if (!geometry.valid_)
    {
        return {geometry, Vec3::zero(), 0.0, 0.0, 0.0, 0.0};
    }

    const Real density = nonZeroDensity(rawDensity);
    const Real wallAlignment = math::dot(geometry.neighborDirection_, wallNormal);
    const Real normalSign = wallAlignment < 0.0 ? -1.0 : (wallAlignment > 0.0 ? 1.0 : 0.0);
    const Vec3 riemannPlaneNormal = -normalSign * wallNormal;
    const Real leftNormalVelocity = math::dot(velocity, riemannPlaneNormal);
    const Real wallNormalVelocity = math::dot(wallVelocity, riemannPlaneNormal);
    const Real rightNormalVelocity = 2.0 * wallNormalVelocity - leftNormalVelocity;
    const Real externalAcceleration = math::dot(gravity - wallAcceleration, -geometry.neighborDirection_);
    const Real positiveExternalAcceleration = externalAcceleration > 0.0 ? externalAcceleration : 0.0;
    const Real wallPressure = pressure + density * geometry.distance_ * positiveExternalAcceleration;
    return {geometry, riemannPlaneNormal, density, leftNormalVelocity, rightNormalVelocity, wallPressure};
}

} // namespace detail

/** Computes the pressure part of one fluid-fluid acceleration pair. */
FUNDEM_MATH_HD inline Vec3 calculateFluidPressureAcceleration(const Vec3& position,
                                                              const Vec3& velocity,
                                                              Real rawDensity,
                                                              Real pressure,
                                                              const Vec3& neighborPosition,
                                                              const Vec3& neighborVelocity,
                                                              Real rawNeighborDensity,
                                                              Real neighborPressure,
                                                              Real neighborMass,
                                                              Real referenceDensity,
                                                              Real smoothingLength,
                                                              Real soundSpeed) noexcept
{
    const detail::SPHInteractionGeometry geometry = detail::makeSPHInteractionGeometry(position, neighborPosition, smoothingLength);
    if (!geometry.valid_)
    {
        return Vec3::zero();
    }

    const Vec3 riemannPlaneNormal = -geometry.neighborDirection_;
    const Real leftNormalVelocity = math::dot(velocity, riemannPlaneNormal);
    const Real rightNormalVelocity = math::dot(neighborVelocity, riemannPlaneNormal);
    const Real interfacePressure = computeRiemannInterfacePressure(pressure, neighborPressure, leftNormalVelocity, rightNormalVelocity, referenceDensity, soundSpeed);
    const Real neighborVolume = neighborMass / nonZeroDensity(rawNeighborDensity);
    return computePressureAcceleration(neighborVolume, interfacePressure, nonZeroDensity(rawDensity), geometry.kernelGradient_);
}

/** Computes the viscosity part of one fluid-fluid acceleration pair. */
FUNDEM_MATH_HD inline Vec3 calculateFluidViscousAcceleration(const Vec3& position,
                                                             const Vec3& velocity,
                                                             Real rawDensity,
                                                             const Vec3& neighborPosition,
                                                             const Vec3& neighborVelocity,
                                                             Real rawNeighborDensity,
                                                             Real neighborMass,
                                                             Real dynamicViscosity,
                                                             Real smoothingLength) noexcept
{
    const detail::SPHInteractionGeometry geometry = detail::makeSPHInteractionGeometry(position, neighborPosition, smoothingLength);
    if (!geometry.valid_)
    {
        return Vec3::zero();
    }

    const Real neighborVolume = neighborMass / nonZeroDensity(rawNeighborDensity);
    const Real radialKernelDerivative = geometry.kernelGradientScale_ * geometry.distance_;
    return computeViscousAcceleration(velocity, neighborVelocity, dynamicViscosity, neighborVolume, nonZeroDensity(rawDensity), geometry.distance_, radialKernelDerivative);
}

/** Computes one fluid neighbor's Riemann density-rate contribution. */
FUNDEM_MATH_HD inline Real calculateFluidDensityRate(const Vec3& position,
                                                     const Vec3& velocity,
                                                     Real rawDensity,
                                                     Real pressure,
                                                     const Vec3& neighborPosition,
                                                     const Vec3& neighborVelocity,
                                                     Real rawNeighborDensity,
                                                     Real neighborPressure,
                                                     Real neighborMass,
                                                     Real referenceDensity,
                                                     Real smoothingLength,
                                                     Real soundSpeed) noexcept
{
    const detail::SPHInteractionGeometry geometry = detail::makeSPHInteractionGeometry(position, neighborPosition, smoothingLength);
    if (!geometry.valid_)
    {
        return 0.0;
    }

    const Vec3 riemannPlaneNormal = -geometry.neighborDirection_;
    const Real leftNormalVelocity = math::dot(velocity, riemannPlaneNormal);
    const Real rightNormalVelocity = math::dot(neighborVelocity, riemannPlaneNormal);
    const Real interfaceNormalVelocity = computeRiemannInterfaceNormalVelocity(pressure, neighborPressure, leftNormalVelocity, rightNormalVelocity, referenceDensity, soundSpeed);
    const Vec3 interfaceVelocity = computeRiemannInterfaceVelocity(velocity, neighborVelocity, riemannPlaneNormal, leftNormalVelocity, rightNormalVelocity, interfaceNormalVelocity);
    const Real neighborVolume = neighborMass / nonZeroDensity(rawNeighborDensity);
    return computeDensityRateContribution(nonZeroDensity(rawDensity), neighborVolume, velocity, interfaceVelocity, geometry.kernelGradient_);
}

/** Computes the one-sided Riemann pressure acceleration from a moving wall sample. */
FUNDEM_MATH_HD inline Vec3 calculateWallPressureAcceleration(const Vec3& position,
                                                             const Vec3& velocity,
                                                             Real rawDensity,
                                                             Real pressure,
                                                             const Vec3& wallPosition,
                                                             const Vec3& wallNormal,
                                                             const Vec3& wallVelocity,
                                                             const Vec3& wallAcceleration,
                                                             Real boundaryVolume,
                                                             Real referenceDensity,
                                                             Real smoothingLength,
                                                             Real soundSpeed,
                                                             const Vec3& gravity) noexcept
{
    const detail::WallRiemannState wall = detail::makeWallRiemannState(position, velocity, rawDensity, pressure, wallPosition, wallNormal, wallVelocity, wallAcceleration, smoothingLength, gravity);
    if (!wall.geometry_.valid_)
    {
        return Vec3::zero();
    }

    const Real interfacePressure = computeRiemannInterfacePressure(pressure, wall.wallPressure_, wall.leftNormalVelocity_, wall.rightNormalVelocity_, referenceDensity, soundSpeed);
    return computePressureAcceleration(boundaryVolume, interfacePressure, wall.density_, wall.geometry_.kernelGradient_);
}

/** Computes the no-slip viscosity acceleration from a mirrored wall velocity. */
FUNDEM_MATH_HD inline Vec3 calculateWallViscousAcceleration(const Vec3& position,
                                                            const Vec3& velocity,
                                                            Real rawDensity,
                                                            const Vec3& wallPosition,
                                                            const Vec3& wallVelocity,
                                                            Real boundaryVolume,
                                                            Real dynamicViscosity,
                                                            Real smoothingLength) noexcept
{
    const detail::SPHInteractionGeometry geometry = detail::makeSPHInteractionGeometry(position, wallPosition, smoothingLength);
    if (!geometry.valid_)
    {
        return Vec3::zero();
    }

    const Vec3 mirroredWallVelocity = 2.0 * wallVelocity - velocity;
    const Real radialKernelDerivative = geometry.kernelGradientScale_ * geometry.distance_;
    return computeViscousAcceleration(velocity, mirroredWallVelocity, dynamicViscosity, boundaryVolume, nonZeroDensity(rawDensity), geometry.distance_, radialKernelDerivative);
}

/** Computes the one-sided wall contribution to fluid density rate. */
FUNDEM_MATH_HD inline Real calculateWallDensityRate(const Vec3& position,
                                                    const Vec3& velocity,
                                                    Real rawDensity,
                                                    Real pressure,
                                                    const Vec3& wallPosition,
                                                    const Vec3& wallNormal,
                                                    const Vec3& wallVelocity,
                                                    const Vec3& wallAcceleration,
                                                    Real boundaryVolume,
                                                    Real referenceDensity,
                                                    Real smoothingLength,
                                                    Real soundSpeed,
                                                    const Vec3& gravity) noexcept
{
    const detail::WallRiemannState wall = detail::makeWallRiemannState(position, velocity, rawDensity, pressure, wallPosition, wallNormal, wallVelocity, wallAcceleration, smoothingLength, gravity);
    if (!wall.geometry_.valid_)
    {
        return 0.0;
    }

    const Real interfaceNormalVelocity = computeRiemannInterfaceNormalVelocity(pressure, wall.wallPressure_, wall.leftNormalVelocity_, wall.rightNormalVelocity_, referenceDensity, soundSpeed);
    const Vec3 tangentialVelocity = velocity - wall.leftNormalVelocity_ * wall.riemannPlaneNormal_;
    const Vec3 interfaceVelocity = interfaceNormalVelocity * wall.riemannPlaneNormal_ + tangentialVelocity;
    return computeDensityRateContribution(wall.density_, boundaryVolume, velocity, interfaceVelocity, wall.geometry_.kernelGradient_);
}

/** Reaction force exerted on the wall by one fluid-particle acceleration contribution. */
FUNDEM_MATH_HD inline Vec3 calculateWallForce(Real particleMass, const Vec3& fluidAcceleration) noexcept { return -particleMass * fluidAcceleration; }

/** Returns one neighbor's contribution to the position-divergence free-surface indicator. */
FUNDEM_MATH_HD inline Real freeSurfacePositionDivergenceContribution(const Vec3& separation, Real neighborVolume, Real smoothingLength) noexcept
{
    return -neighborVolume * math::dot(wendlandKernelGradient3D(separation, smoothingLength), separation);
}

/** Classifies a particle as free surface from its position divergence. */
FUNDEM_MATH_HD inline bool isFreeSurface(Real positionDivergence) noexcept { return positionDivergence < 2.25; }

/** Returns the kernel-summation density correction without reducing density below reference. */
FUNDEM_MATH_HD inline Real reinitializedDensity(Real kernelSum, Real referenceKernelSum, Real referenceDensity) noexcept
{
    if (referenceKernelSum <= math::defaultTolerance)
    {
        return referenceDensity;
    }
    const Real summationDensity = referenceDensity * kernelSum / referenceKernelSum;
    return summationDensity > referenceDensity ? summationDensity : referenceDensity;
}

/** Advances density, applies emergency numerical bounds, and updates pressure. */
FUNDEM_MATH_HD inline void integrateDensityAndPressure(Real& density, Real& pressure, Real densityRate, Real referenceDensity, Real soundSpeed, Real timeStep) noexcept
{
    const Real updatedDensity = density + timeStep * densityRate;
    const Real densityFloor = numericalSPHDensityFloorRatio * referenceDensity;
    const Real densityCeiling = numericalSPHDensityCeilingRatio * referenceDensity;
    density = !math::isFinite(updatedDensity) || updatedDensity <= densityFloor ? densityFloor : (updatedDensity < densityCeiling ? updatedDensity : densityCeiling);
    pressure = computePressureFromDensity(density, referenceDensity, soundSpeed);
}

/** Transforms one boundary sample's local position and normal into world space. */
FUNDEM_MATH_HD inline void calculateVirtualParticlePositionAndNormal(Vec3& position,
                                                                     Vec3& normal,
                                                                     const Vec3& localPosition,
                                                                     const Vec3& localNormal,
                                                                     const Vec3& ownerPosition,
                                                                     const Quaternion& ownerOrientation) noexcept
{
    position = ownerPosition + math::rotateUnit(ownerOrientation, localPosition);
    normal = math::rotateUnit(ownerOrientation, localNormal);
}

/** Converts one virtual-particle reaction force into owner force and torque. */
FUNDEM_MATH_HD inline void calculateVirtualParticleForceAndTorque(Vec3& force, Vec3& torque, const Vec3& localPosition, const Quaternion& ownerOrientation, const Vec3& virtualParticleForce) noexcept
{
    force = virtualParticleForce;
    torque = math::cross(math::rotateUnit(ownerOrientation, localPosition), virtualParticleForce);
}

/** Accumulates boundary-sample kinematics implied by one owner rigid body. */
FUNDEM_MATH_HD inline void accumulateVirtualParticleVelocityAndAcceleration(Vec3& velocity,
                                                                            Vec3& acceleration,
                                                                            const Vec3& localPosition,
                                                                            const Quaternion& ownerOrientation,
                                                                            const Vec3& ownerVelocity,
                                                                            const Vec3& ownerAngularVelocity,
                                                                            const Vec3& ownerForce = Vec3::zero(),
                                                                            const Vec3& ownerTorque = Vec3::zero(),
                                                                            Real ownerInverseMass = 0.0,
                                                                            const Mat3& ownerInertiaTensor = Mat3::zero(),
                                                                            const Mat3& ownerInverseInertiaTensor = Mat3::zero(),
                                                                            const Vec3& gravity = Vec3{0.0, 0.0, -9.81}) noexcept
{
    const Vec3 leverArm = math::rotateUnit(ownerOrientation, localPosition);
    velocity += ownerVelocity + math::cross(ownerAngularVelocity, leverArm);
    acceleration += math::cross(ownerAngularVelocity, math::cross(ownerAngularVelocity, leverArm));
    if (ownerInverseMass <= 0.0)
    {
        return;
    }
    const Vec3 angularMomentum = math::rotateUnit(ownerOrientation, ownerInertiaTensor * math::inverseRotateUnit(ownerOrientation, ownerAngularVelocity));
    const Vec3 angularAcceleration =
        math::rotateUnit(ownerOrientation, ownerInverseInertiaTensor * math::inverseRotateUnit(ownerOrientation, ownerTorque - math::cross(ownerAngularVelocity, angularMomentum)));
    const Vec3 linearAcceleration = ownerInverseMass * ownerForce + gravity;
    acceleration += linearAcceleration + math::cross(angularAcceleration, leverArm);
}

} // namespace fundem::execution
