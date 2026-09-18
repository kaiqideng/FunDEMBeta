/**
 * @file SPHInteractionTypes.h
 * @brief Defines shared WCSPH interaction parameters and reduced statistics.
 */
#pragma once

#include "math/Vector3.h"

namespace fundem
{

/** Fluid-state extrema and validity counters produced by either backend. */
struct SPHStateStatistics {
    math::Real maximumVelocity_{0.0};     ///< Largest finite fluid speed.
    math::Real maximumAcceleration_{0.0}; ///< Largest finite acceleration including gravity.
    math::Real minimumDensity_{0.0};      ///< Smallest finite positive density.
    math::Real maximumDensity_{0.0};      ///< Largest finite density.
    int invalidValueCount_{0};            ///< Number of particles containing invalid state.
};

/** Virtual-boundary kinematic extrema and validity counters. */
struct SPHKinematicsStatistics {
    math::Real maximumVelocity_{0.0};     ///< Largest finite boundary speed.
    math::Real maximumAcceleration_{0.0}; ///< Largest finite boundary acceleration.
    int invalidValueCount_{0};            ///< Number of invalid boundary samples.
};

/** Immutable physical values shared by interaction stages. */
struct SPHInteractionParameters {
    math::Real smoothingLength_{0.0};        ///< Kernel smoothing length for this stage.
    math::Real particleMass_{0.0};           ///< Uniform fluid-particle mass.
    math::Real referenceDensity_{0.0};       ///< Equation-of-state reference density.
    math::Real latticeKernelSum3D_{0.0};     ///< Reference density reinitialization normalization.
    math::Real soundSpeed_{0.0};             ///< Artificial speed of sound.
    math::Real dynamicViscosity_{0.0};       ///< Uniform dynamic viscosity.
    math::Vec3 gravity_{math::Vec3::zero()}; ///< Uniform gravitational acceleration.
};

} // namespace fundem
