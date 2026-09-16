/**
 * @file vtuOutput.h
 * @brief Assembles default and user-selected VTU fields for simulation objects.
 */
#pragma once

#include "interaction/bond.h"
#include "interaction/contact.h"
#include "particle/LSParticle.h"
#include "particle/SPHParticle.h"
#include "particle/particle.h"
#include "vtuWriter.h"

#include <string>
#include <vector>

namespace fundem
{

/** Optional fields for spherical-particle VTU output. */
enum class sphereVTUField
{
    orientation,
    angularVelocity,
    force,
    torque,
    inverseMass,
    inertiaTensor,
    materialIndex,
    kineticEnergy,
    gravitationalPotentialEnergy
};

/** Optional fields for expanded level-set surface-node VTU output. */
enum class LSParticleVTUField
{
    particleIndex,
    distanceToCenter,
    orientation,
    angularVelocity,
    force,
    torque,
    inertiaTensor,
    boundingRadius,
    materialIndex,
    volume,
    geometryIndex,
    surfaceNodeArea,
    kineticEnergy,
    gravitationalPotentialEnergy
};

/** Optional fields for SPH-particle VTU output. */
enum class SPHParticleVTUField
{
    force,
    pressure,
    densityRate,
    freeSurface,
    kineticEnergy,
    gravitationalPotentialEnergy
};

/** Optional fields for contact-point VTU output. */
enum class contactVTUField
{
    force,
    torque,
    masterParticleIndex,
    slaveParticleIndex,
    overlap,
    area,
    effectiveMass,
    effectiveRadius,
    normalForceMagnitude,
    normalElasticEnergy,
    slidingElasticEnergy,
    rollingElasticEnergy,
    torsionalElasticEnergy,
    slidingSpringDeformation,
    rollingSpringDeformation,
    torsionalSpringDeformation
};
/** Optional fields for bond-point VTU output. */

enum class bondVTUField
{
    force,
    masterTorque,
    slaveTorque,
    masterParticleIndex,
    slaveParticleIndex,
    equivalentLength,
    normalElasticEnergy,
    shearElasticEnergy,
    bendingElasticEnergy,
    torsionalElasticEnergy,
    masterEndpointNormal,
    masterEndpointTangent1,
    masterEndpointTangent2,
    slaveEndpointNormal,
    slaveEndpointTangent1,
    slaveEndpointTangent2,
    crossSectionArea,
    damageFactor,
    modeICriticalEnergy,
    modeIICriticalEnergy,
    modeMixityExponent,
    damageInitiationRatio,
    maximumEnergyReleaseRatio
};

/** Creates a sphere writer with position, radius, and velocity plus selected fields. */
vtuWriter makeSphereVTUWriter(std::string fileName, const particleContainer& spheres, const math::Vec3& gravity, const std::vector<sphereVTUField>& fields = {});

/** Creates an expanded finite-mass LS-particle surface writer. */
vtuWriter makeLSParticleVTUWriter(std::string fileName, const LSParticleContainer& particles, const math::Vec3& gravity, const std::vector<LSParticleVTUField>& fields = {});
/** Creates an expanded infinite-mass LS-particle surface writer. */
vtuWriter makeFixedLSParticleVTUWriter(std::string fileName, const LSParticleContainer& particles, const math::Vec3& gravity, const std::vector<LSParticleVTUField>& fields = {});

/** Creates an SPH writer with position, velocity, mass, and density plus selected fields. */
vtuWriter makeSPHParticleVTUWriter(std::string fileName, const SPHParticleContainer& particles, const math::Vec3& gravity, const std::vector<SPHParticleVTUField>& fields = {});

/** Creates a contact writer with point and normal plus selected fields. */
vtuWriter makeContactVTUWriter(std::string fileName, const contactContainer& contacts, const std::vector<contactVTUField>& fields = {});

/** Creates a writer for sphere-sphere bonds. */
vtuWriter makeBondVTUWriter(std::string fileName, const bondContainer& bonds, const particleContainer& spheres, const std::vector<bondVTUField>& fields = {});

/** Creates a writer for LS-LS bonds. */
vtuWriter makeBondVTUWriter(std::string fileName, const bondContainer& bonds, const LSParticleContainer& LSParticles, const std::vector<bondVTUField>& fields = {});

/** Creates a writer for sphere-LS bonds. */
vtuWriter makeBondVTUWriter(std::string fileName,
                            const bondContainer& bonds,
                            const particleContainer& masterSpheres,
                            const LSParticleContainer& slaveLSParticles,
                            const std::vector<bondVTUField>& fields = {});

} // namespace fundem
