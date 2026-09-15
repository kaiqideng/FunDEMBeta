/**
 * @file contactForceKernel.cuh
 * @brief Declares contact-force device views and streamed CUDA launch interfaces.
 */
#pragma once

#include "interaction/contact.h"
#include "material/material.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{
namespace detail
{

/** Contact response arrays read and updated by the force kernel. */
struct contactDeviceView {
    const int* masterParticleIndex_{nullptr};         ///< Master particle indices.
    const int* slaveParticleIndex_{nullptr};          ///< Slave particle indices.
    math::Vec3* force_{nullptr};                      ///< Resultant forces on masters.
    math::Vec3* torque_{nullptr};                     ///< Resultant contact torques on masters.
    const math::Vec3* point_{nullptr};                ///< World contact points.
    const math::Vec3* normal_{nullptr};               ///< Slave-to-master normals.
    const math::Real* overlap_{nullptr};              ///< Normal penetrations.
    const math::Real* area_{nullptr};                 ///< Contact or patch areas.
    const math::Real* effectiveMass_{nullptr};        ///< Effective damping masses.
    const math::Real* effectiveRadius_{nullptr};      ///< Effective moment radii.
    math::Real* normalForceMagnitude_{nullptr};       ///< Computed normal-force magnitudes.
    math::Real* normalElasticEnergy_{nullptr};        ///< Normal elastic energies.
    math::Real* slidingElasticEnergy_{nullptr};       ///< Sliding elastic energies.
    math::Real* rollingElasticEnergy_{nullptr};       ///< Rolling elastic energies.
    math::Real* torsionalElasticEnergy_{nullptr};     ///< Torsional elastic energies.
    math::Vec3* slidingSpringDeformation_{nullptr};   ///< Sliding history state.
    math::Vec3* rollingSpringDeformation_{nullptr};   ///< Rolling history state.
    math::Vec3* torsionalSpringDeformation_{nullptr}; ///< Torsional history state.
};

/** Particle arrays required to evaluate and scatter contact response. */
struct contactParticleDeviceView {
    math::Vec3* force_{nullptr};                 ///< Accumulated particle forces.
    math::Vec3* torque_{nullptr};                ///< Accumulated particle torques.
    const math::Vec3* position_{nullptr};        ///< World positions.
    const math::Vec3* velocity_{nullptr};        ///< Linear velocities.
    const math::Vec3* angularVelocity_{nullptr}; ///< Angular velocities.
    const math::Real* inverseMass_{nullptr};     ///< Reciprocal masses.
    const int* materialIndex_{nullptr};          ///< Material-container indices.
};

/** Read-only material arrays consumed by the unified contact law. */
struct contactMaterialDeviceView {
    const math::Real* normalStiffness_{nullptr};              ///< Normal stiffness values.
    const math::Real* slidingStiffness_{nullptr};             ///< Sliding/shear stiffness values.
    const math::Real* rollingStiffness_{nullptr};             ///< Rolling stiffness values.
    const math::Real* torsionalStiffness_{nullptr};           ///< Torsional stiffness values.
    const math::Real* slidingFrictionCoefficient_{nullptr};   ///< Sliding friction coefficients.
    const math::Real* rollingFrictionCoefficient_{nullptr};   ///< Rolling friction coefficients.
    const math::Real* torsionalFrictionCoefficient_{nullptr}; ///< Torsional friction coefficients.
    const math::Real* restitutionCoefficient_{nullptr};       ///< Normal restitution coefficients.
    const materialType* type_{nullptr};                       ///< Runtime stiffness conventions.
};

/** Launches one contact-force kernel without synchronizing @p stream. */
void launchContactForce(contactDeviceView contacts,
                        contactParticleDeviceView masterParticles,
                        contactParticleDeviceView slaveParticles,
                        contactMaterialDeviceView materials,
                        math::Real timeStep,
                        int contactCount,
                        cudaStream_t stream);

inline contactDeviceView makeContactDeviceView(contactContainer& contacts) noexcept
{
    contactDeviceView view;
    view.masterParticleIndex_ = contacts.device<contact::masterParticleIndexField>();
    view.slaveParticleIndex_ = contacts.device<contact::slaveParticleIndexField>();
    view.force_ = contacts.device<contact::forceField>();
    view.torque_ = contacts.device<contact::torqueField>();
    view.point_ = contacts.device<contact::pointField>();
    view.normal_ = contacts.device<contact::normalField>();
    view.overlap_ = contacts.device<contact::overlapField>();
    view.area_ = contacts.device<contact::areaField>();
    view.effectiveMass_ = contacts.device<contact::effectiveMassField>();
    view.effectiveRadius_ = contacts.device<contact::effectiveRadiusField>();
    view.normalForceMagnitude_ = contacts.device<contact::normalForceMagnitudeField>();
    view.normalElasticEnergy_ = contacts.device<contact::normalElasticEnergyField>();
    view.slidingElasticEnergy_ = contacts.device<contact::slidingElasticEnergyField>();
    view.rollingElasticEnergy_ = contacts.device<contact::rollingElasticEnergyField>();
    view.torsionalElasticEnergy_ = contacts.device<contact::torsionalElasticEnergyField>();
    view.slidingSpringDeformation_ = contacts.device<contact::slidingSpringDeformationField>();
    view.rollingSpringDeformation_ = contacts.device<contact::rollingSpringDeformationField>();
    view.torsionalSpringDeformation_ = contacts.device<contact::torsionalSpringDeformationField>();
    return view;
}

template <class ParticleStorage> contactParticleDeviceView makeContactParticleDeviceView(ParticleStorage& particles) noexcept
{
    contactParticleDeviceView view;
    view.force_ = particles.template device<rigidBody::forceField>();
    view.torque_ = particles.template device<rigidBody::torqueField>();
    view.position_ = particles.template device<rigidBody::positionField>();
    view.velocity_ = particles.template device<rigidBody::velocityField>();
    view.angularVelocity_ = particles.template device<rigidBody::angularVelocityField>();
    view.inverseMass_ = particles.template device<rigidBody::inverseMassField>();
    view.materialIndex_ = particles.template device<particle::materialIndexField>();
    return view;
}

inline contactMaterialDeviceView makeContactMaterialDeviceView(const materialContainer& materials) noexcept
{
    contactMaterialDeviceView view;
    view.normalStiffness_ = materials.device<material::normalStiffnessField>();
    view.slidingStiffness_ = materials.device<material::slidingStiffnessField>();
    view.rollingStiffness_ = materials.device<material::rollingStiffnessField>();
    view.torsionalStiffness_ = materials.device<material::torsionalStiffnessField>();
    view.slidingFrictionCoefficient_ = materials.device<material::slidingFrictionCoefficientField>();
    view.rollingFrictionCoefficient_ = materials.device<material::rollingFrictionCoefficientField>();
    view.torsionalFrictionCoefficient_ = materials.device<material::torsionalFrictionCoefficientField>();
    view.restitutionCoefficient_ = materials.device<material::restitutionCoefficientField>();
    view.type_ = materials.device<material::typeField>();
    return view;
}

} // namespace detail

/**
 * Builds views and launches contact response for separate master/slave
 * containers. Every device operation is ordered on @p stream.
 */
template <class MasterParticleStorage, class SlaveParticleStorage>
void launchContactForce(contactContainer& contacts,
                        MasterParticleStorage& masterParticles,
                        SlaveParticleStorage& slaveParticles,
                        const materialContainer& materials,
                        math::Real timeStep,
                        cudaStream_t stream = nullptr)
{
    detail::launchContactForce(detail::makeContactDeviceView(contacts),
                               detail::makeContactParticleDeviceView(masterParticles),
                               detail::makeContactParticleDeviceView(slaveParticles),
                               detail::makeContactMaterialDeviceView(materials),
                               timeStep,
                               static_cast<int>(contacts.deviceSize()),
                               stream);
}

/** Same-container convenience wrapper for `launchContactForce`. */
template <class ParticleStorage> void launchContactForce(contactContainer& contacts, ParticleStorage& particles, const materialContainer& materials, math::Real timeStep, cudaStream_t stream = nullptr)
{
    launchContactForce(contacts, particles, particles, materials, timeStep, stream);
}

} // namespace fundem::cuda
