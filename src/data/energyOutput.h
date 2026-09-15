/**
 * @file energyOutput.h
 * @brief Defines solid-system energy records and DAT output assembly.
 */
#pragma once

#include "interaction/bond.h"
#include "interaction/contact.h"
#include "particle/LSParticle.h"
#include "particle/SPHParticle.h"
#include "particle/particle.h"

#include <string>

namespace fundem
{

/** Aggregated energy components of the solid DEM subsystem. */
struct energyRecord {
    using Real = math::Real;

    Real kineticEnergy_{0.0};                 ///< Translational plus rotational kinetic energy.
    Real gravitationalPotentialEnergy_{0.0};  ///< Gravitational potential energy of finite-mass solids.
    Real contactNormalElasticEnergy_{0.0};    ///< Elastic energy stored in normal contacts.
    Real contactSlidingElasticEnergy_{0.0};   ///< Elastic energy stored in sliding springs.
    Real contactRollingElasticEnergy_{0.0};   ///< Elastic energy stored in rolling springs.
    Real contactTorsionalElasticEnergy_{0.0}; ///< Elastic energy stored in torsional contact springs.
    Real bondNormalElasticEnergy_{0.0};       ///< Normal elastic energy stored in bonds.
    Real bondShearElasticEnergy_{0.0};        ///< Shear elastic energy stored in bonds.
    Real bondBendingElasticEnergy_{0.0};      ///< Bending elastic energy stored in bonds.
    Real bondTorsionalElasticEnergy_{0.0};    ///< Torsional elastic energy stored in bonds.

    /** Returns the sum of all stored components. */
    Real totalEnergy() const noexcept;
};

/** Adds spherical-particle mechanical energy to @p result. */
void addParticleEnergy(energyRecord& result, const particleContainer& particles, const math::Vec3& gravity) noexcept;
/** Adds finite-mass level-set-particle mechanical energy to @p result. */
void addParticleEnergy(energyRecord& result, const LSParticleContainer& particles, const math::Vec3& gravity) noexcept;
/** Adds SPH-particle mechanical energy to @p result when explicitly requested. */
void addParticleEnergy(energyRecord& result, const SPHParticleContainer& particles, const math::Vec3& gravity) noexcept;
/** Adds elastic contact energy to @p result. */
void addContactEnergy(energyRecord& result, const contactContainer& contacts) noexcept;
/** Adds intact-bond elastic energy to @p result. */
void addBondEnergy(energyRecord& result, const bondContainer& bonds) noexcept;
/** Appends one timestamped aggregate-energy row to a DAT file. */
void appendEnergyDAT(const std::string& fileName, math::Real time, const energyRecord& energy, bool appendExisting);
} // namespace fundem
