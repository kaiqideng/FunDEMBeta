/**
 * @file bondForceKernel.cuh
 * @brief Declares bond device views and streamed CUDA launch interfaces.
 */
#pragma once

#include "interaction/bond.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{
namespace detail
{

/** Bond state arrays read and updated by the cohesive-force kernel. */
struct bondDeviceView {
    const int* masterParticleIndex_{nullptr};                ///< Master particle indices.
    const int* slaveParticleIndex_{nullptr};                 ///< Slave particle indices.
    math::Vec3* force_{nullptr};                             ///< Forces applied to masters.
    math::Vec3* torque_{nullptr};                            ///< Bond-point bending/torsion moments.
    math::Vec3* point_{nullptr};                             ///< Current bond-center positions.
    math::Vec3* normal_{nullptr};                            ///< Current slave-to-master axes.
    const math::Real* equivalentLength_{nullptr};            ///< Reference bond lengths.
    const math::Real* coefficientB1_{nullptr};               ///< Axial response coefficients.
    const math::Real* coefficientB2_{nullptr};               ///< Shear response coefficients.
    const math::Real* coefficientB3_{nullptr};               ///< Bending response coefficients.
    const math::Real* coefficientB4_{nullptr};               ///< Torsional response coefficients.
    const math::Real* normalStiffness_{nullptr};             ///< Original undamaged normal stiffnesses.
    const math::Real* shearStiffness_{nullptr};              ///< Original undamaged shear stiffnesses.
    const math::Real* bendingStiffness_{nullptr};            ///< Original undamaged bending stiffnesses.
    math::Real* normalElasticEnergy_{nullptr};               ///< Axial elastic energies.
    math::Real* shearElasticEnergy_{nullptr};                ///< Shear elastic energies.
    math::Real* bendingElasticEnergy_{nullptr};              ///< Bending elastic energies.
    math::Real* torsionalElasticEnergy_{nullptr};            ///< Torsional elastic energies.
    const math::Vec3* masterEndpointLocalNormal_{nullptr};   ///< Master endpoint reference normal.
    const math::Vec3* masterEndpointLocalTangent1_{nullptr}; ///< Master endpoint first reference tangent.
    const math::Vec3* masterEndpointLocalTangent2_{nullptr}; ///< Master endpoint second reference tangent.
    const math::Vec3* masterEndpointLocalPosition_{nullptr}; ///< Master endpoint position in body coordinates.
    const math::Vec3* slaveEndpointLocalNormal_{nullptr};    ///< Slave endpoint reference normal.
    const math::Vec3* slaveEndpointLocalTangent1_{nullptr};  ///< Slave endpoint first reference tangent.
    const math::Vec3* slaveEndpointLocalTangent2_{nullptr};  ///< Slave endpoint second reference tangent.
    const math::Vec3* slaveEndpointLocalPosition_{nullptr};  ///< Slave endpoint position in body coordinates.
    const math::Real* crossSectionArea_{nullptr};            ///< Nominal cross-sectional areas; zero disables fracture.
    const fractureMode* fractureMode_{nullptr};              ///< Unbreakable, brittle, or progressive damage selection.
    math::Real* damageFactor_{nullptr};                      ///< Irreversible damage variables in `[0, 1]`.
    const math::Real* modeICriticalEnergy_{nullptr};         ///< Critical mode-I fracture energies.
    const math::Real* modeIICriticalEnergy_{nullptr};        ///< Critical mode-II fracture energies.
    const math::Real* modeMixityExponent_{nullptr};          ///< BK mixed-mode exponents.
    const math::Real* tensileStrength_{nullptr};             ///< Tensile strengths in Pa.
    const math::Real* cohesion_{nullptr};                    ///< Cohesive shear strengths in Pa.
    const math::Real* internalFrictionCoefficient_{nullptr}; ///< Internal friction coefficients.
    int* damageStarted_{nullptr};                           ///< Strength-onset history flags.
    math::Real* onsetEffectiveStress_{nullptr};              ///< Bond-section stress norms at onset, not energy-conjugate tractions, Pa.
    math::Real* onsetEffectiveDisplacement_{nullptr};        ///< Euclidean energy-equivalent separations at onset, m.
    math::Real* maximumEffectiveDisplacement_{nullptr};      ///< Historical maximum Euclidean separations, m.
    math::Real* failureEffectiveDisplacement_{nullptr};      ///< Onset-radial rupture separations, m.
    math::Real* effectiveFractureEnergy_{nullptr};           ///< Frozen effective mixed-mode fracture energies.
};

/** Particle arrays used to transform endpoints and scatter bond loads. */
struct bondParticleDeviceView {
    math::Vec3* force_{nullptr};                   ///< Accumulated particle forces.
    math::Vec3* torque_{nullptr};                  ///< Accumulated particle torques.
    const math::Vec3* position_{nullptr};          ///< World particle positions.
    const math::Quaternion* orientation_{nullptr}; ///< Body-to-world orientations.
};

/** Launches one bond-force kernel without synchronizing @p stream. */
void launchBondForce(bondDeviceView bonds, bondParticleDeviceView masterParticles, bondParticleDeviceView slaveParticles, int bondCount, cudaStream_t stream);

inline bondDeviceView makeBondDeviceView(bondContainer& bonds) noexcept
{
    bondDeviceView view;
    view.masterParticleIndex_ = bonds.device<bond::masterParticleIndexField>();
    view.slaveParticleIndex_ = bonds.device<bond::slaveParticleIndexField>();
    view.force_ = bonds.device<bond::forceField>();
    view.torque_ = bonds.device<bond::torqueField>();
    view.point_ = bonds.device<bond::pointField>();
    view.normal_ = bonds.device<bond::normalField>();
    view.equivalentLength_ = bonds.device<bond::equivalentLengthField>();
    view.coefficientB1_ = bonds.device<bond::coefficientB1Field>();
    view.coefficientB2_ = bonds.device<bond::coefficientB2Field>();
    view.coefficientB3_ = bonds.device<bond::coefficientB3Field>();
    view.coefficientB4_ = bonds.device<bond::coefficientB4Field>();
    view.normalStiffness_ = bonds.device<bond::normalStiffnessField>();
    view.shearStiffness_ = bonds.device<bond::shearStiffnessField>();
    view.bendingStiffness_ = bonds.device<bond::bendingStiffnessField>();
    view.normalElasticEnergy_ = bonds.device<bond::normalElasticEnergyField>();
    view.shearElasticEnergy_ = bonds.device<bond::shearElasticEnergyField>();
    view.bendingElasticEnergy_ = bonds.device<bond::bendingElasticEnergyField>();
    view.torsionalElasticEnergy_ = bonds.device<bond::torsionalElasticEnergyField>();
    view.masterEndpointLocalNormal_ = bonds.device<bond::masterEndpointLocalNormalField>();
    view.masterEndpointLocalTangent1_ = bonds.device<bond::masterEndpointLocalTangent1Field>();
    view.masterEndpointLocalTangent2_ = bonds.device<bond::masterEndpointLocalTangent2Field>();
    view.masterEndpointLocalPosition_ = bonds.device<bond::masterEndpointLocalPositionField>();
    view.slaveEndpointLocalNormal_ = bonds.device<bond::slaveEndpointLocalNormalField>();
    view.slaveEndpointLocalTangent1_ = bonds.device<bond::slaveEndpointLocalTangent1Field>();
    view.slaveEndpointLocalTangent2_ = bonds.device<bond::slaveEndpointLocalTangent2Field>();
    view.slaveEndpointLocalPosition_ = bonds.device<bond::slaveEndpointLocalPositionField>();
    view.crossSectionArea_ = bonds.device<bond::crossSectionAreaField>();
    view.fractureMode_ = bonds.device<bond::fractureModeField>();
    view.damageFactor_ = bonds.device<bond::damageFactorField>();
    view.modeICriticalEnergy_ = bonds.device<bond::modeICriticalEnergyField>();
    view.modeIICriticalEnergy_ = bonds.device<bond::modeIICriticalEnergyField>();
    view.modeMixityExponent_ = bonds.device<bond::modeMixityExponentField>();
    view.tensileStrength_ = bonds.device<bond::tensileStrengthField>();
    view.cohesion_ = bonds.device<bond::cohesionField>();
    view.internalFrictionCoefficient_ = bonds.device<bond::internalFrictionCoefficientField>();
    view.damageStarted_ = bonds.device<bond::damageStartedField>();
    view.onsetEffectiveStress_ = bonds.device<bond::onsetEffectiveStressField>();
    view.onsetEffectiveDisplacement_ = bonds.device<bond::onsetEffectiveDisplacementField>();
    view.maximumEffectiveDisplacement_ = bonds.device<bond::maximumEffectiveDisplacementField>();
    view.failureEffectiveDisplacement_ = bonds.device<bond::failureEffectiveDisplacementField>();
    view.effectiveFractureEnergy_ = bonds.device<bond::effectiveFractureEnergyField>();
    return view;
}

template <class ParticleStorage> bondParticleDeviceView makeBondParticleDeviceView(ParticleStorage& particles) noexcept
{
    bondParticleDeviceView view;
    view.force_ = particles.template device<rigidBody::forceField>();
    view.torque_ = particles.template device<rigidBody::torqueField>();
    view.position_ = particles.template device<rigidBody::positionField>();
    view.orientation_ = particles.template device<rigidBody::orientationField>();
    return view;
}

} // namespace detail

/**
 * Builds device views and evaluates bonds connecting two particle storages.
 * All work is ordered on @p stream.
 */
template <class MasterParticleStorage, class SlaveParticleStorage>
void launchBondForce(bondContainer& bonds, MasterParticleStorage& masterParticles, SlaveParticleStorage& slaveParticles, cudaStream_t stream = nullptr)
{
    detail::launchBondForce(detail::makeBondDeviceView(bonds),
                            detail::makeBondParticleDeviceView(masterParticles),
                            detail::makeBondParticleDeviceView(slaveParticles),
                            static_cast<int>(bonds.deviceSize()),
                            stream);
}

/** Same-container convenience wrapper for `launchBondForce`. */
template <class ParticleStorage> void launchBondForce(bondContainer& bonds, ParticleStorage& particles, cudaStream_t stream = nullptr) { launchBondForce(bonds, particles, particles, stream); }

} // namespace fundem::cuda
