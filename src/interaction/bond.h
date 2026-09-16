/**
 * @file bond.h
 * @brief Defines bonded interactions, endpoint state, damage, loads, and storage.
 */
#pragma once

#include "contact.h"
#include "data/HostAoSDeviceSoA.h"
#include "math/Vector3.h"
#include "particle/LSParticle.h"
#include "particle/rigidBody.h"

namespace fundem
{

/** Particle-container pairing represented by a bond. */
enum class bondType
{
    undefined,
    sphereSphere,
    sphereLSParticle,
    LSParticleLSParticle
};

/** One body-local endpoint frame required to continue a deformed bond. */
struct bondEndpointState {
    math::Vec3 normal_{math::Vec3::unitX()};
    math::Vec3 tangent1_{math::Vec3::unitY()};
    math::Vec3 tangent2_{math::Vec3::unitZ()};
    math::Vec3 position_{math::Vec3::zero()};
};

/** Complete body-local reference geometry required for lossless bond restart. */
struct bondReferenceState {
    bondEndpointState master_;
    bondEndpointState slave_;
};

/**
 * Cohesive bond with body-local endpoint frames, stiffness coefficients,
 * mixed-mode fracture state, particle snapshots, and resultant loads.
 */
class bond : public containerObject
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    bond() = default;
    /** Creates an unconnected bond with a validated reference length. */
    explicit bond(Real equivalentLength) noexcept { setEquivalentLength(equivalentLength); }
    ~bond() = default;
    bond(const bond&) = default;
    bond& operator=(const bond&) = default;
    bond(bond&&) noexcept = default;
    bond& operator=(bond&&) noexcept = default;

    const rigidBody& masterParticle() const noexcept { return masterParticle_; }
    const rigidBody& slaveParticle() const noexcept { return slaveParticle_; }
    /** Reports whether endpoint body snapshots were refreshed for the next force evaluation. */
    bool particleCopiesUpdated() const noexcept { return particleCopiesUpdated_; }
    const Vec3& force() const noexcept { return force_; }
    const Vec3& masterTorque() const noexcept { return masterTorque_; }
    const Vec3& slaveTorque() const noexcept { return slaveTorque_; }
    int masterParticleIndex() const noexcept { return masterParticleIndex_; }
    int slaveParticleIndex() const noexcept { return slaveParticleIndex_; }
    bondType type() const noexcept { return type_; }
    const Vec3& point() const noexcept { return point_; }
    const Vec3& normal() const noexcept { return normal_; }
    Real equivalentLength() const noexcept { return equivalentLength_; }
    Real coefficientB1() const noexcept { return coefficientB1_; }
    Real coefficientB2() const noexcept { return coefficientB2_; }
    Real coefficientB3() const noexcept { return coefficientB3_; }
    Real coefficientB4() const noexcept { return coefficientB4_; }
    Real normalElasticEnergy() const noexcept { return normalElasticEnergy_; }
    Real shearElasticEnergy() const noexcept { return shearElasticEnergy_; }
    Real bendingElasticEnergy() const noexcept { return bendingElasticEnergy_; }
    Real torsionalElasticEnergy() const noexcept { return torsionalElasticEnergy_; }
    const Vec3& masterEndpointLocalNormal() const noexcept { return masterEndpointLocalNormal_; }
    const Vec3& masterEndpointLocalTangent1() const noexcept { return masterEndpointLocalTangent1_; }
    const Vec3& masterEndpointLocalTangent2() const noexcept { return masterEndpointLocalTangent2_; }
    const Vec3& masterEndpointLocalPosition() const noexcept { return masterEndpointLocalPosition_; }
    const Vec3& slaveEndpointLocalNormal() const noexcept { return slaveEndpointLocalNormal_; }
    const Vec3& slaveEndpointLocalTangent1() const noexcept { return slaveEndpointLocalTangent1_; }
    const Vec3& slaveEndpointLocalTangent2() const noexcept { return slaveEndpointLocalTangent2_; }
    const Vec3& slaveEndpointLocalPosition() const noexcept { return slaveEndpointLocalPosition_; }
    Real crossSectionArea() const noexcept { return crossSectionArea_; }
    Real damageFactor() const noexcept { return damageFactor_; }
    Real modeICriticalEnergy() const noexcept { return modeICriticalEnergy_; }
    Real modeIICriticalEnergy() const noexcept { return modeIICriticalEnergy_; }
    Real modeMixityExponent() const noexcept { return modeMixityExponent_; }
    Real damageInitiationRatio() const noexcept { return damageInitiationRatio_; }
    Real maximumEnergyReleaseRatio() const noexcept { return maximumEnergyReleaseRatio_; }
    bondReferenceState referenceState() const noexcept
    {
        return {{masterEndpointLocalNormal_, masterEndpointLocalTangent1_, masterEndpointLocalTangent2_, masterEndpointLocalPosition_},
                {slaveEndpointLocalNormal_, slaveEndpointLocalTangent1_, slaveEndpointLocalTangent2_, slaveEndpointLocalPosition_}};
    }

    /**
     * Sets the positive reference length and invalidates connection and
     * stiffness values that depend on the old length.
     * @return True when @p value is finite and positive.
     */
    bool setEquivalentLength(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= 0.0)
        {
            return false;
        }
        if (equivalentLengthSet_ && value == equivalentLength_)
        {
            return true;
        }
        invalidateLengthDependentValues();
        equivalentLength_ = value;
        equivalentLengthSet_ = true;
        return true;
    }
    /**
     * Configures particle indices, reference point, normal, endpoint frames,
     * and container ownership. The overloads preserve the required particle
     * type ordering.
     */
    bool setConnection(const particleContainer& spheres, int masterParticleIndex, int slaveParticleIndex, const Vec3& normal) noexcept
    {
        return setConnectionForType(spheres, masterParticleIndex, slaveParticleIndex, normal, bondType::sphereSphere);
    }
    bool setConnection(const LSParticleContainer& LSParticles, int masterParticleIndex, int slaveParticleIndex, const Vec3& normal) noexcept
    {
        return setConnectionForType(LSParticles, masterParticleIndex, slaveParticleIndex, normal, bondType::LSParticleLSParticle);
    }
    bool setConnection(const particleContainer& masterSpheres, int masterParticleIndex, const LSParticleContainer& slaveLSParticles, int slaveParticleIndex, const Vec3& normal) noexcept
    {
        return setConnectionForType(masterSpheres, masterParticleIndex, slaveLSParticles, slaveParticleIndex, normal, bondType::sphereLSParticle);
    }
    bool setConnection(const particleContainer& spheres, const contact& singleContact) noexcept { return setConnectionForType(spheres, singleContact, bondType::sphereSphere); }
    bool setConnection(const LSParticleContainer& LSParticles, const contact& singleContact) noexcept { return setConnectionForType(LSParticles, singleContact, bondType::LSParticleLSParticle); }
    bool setConnection(const particleContainer& masterSpheres, const LSParticleContainer& slaveLSParticles, const contact& singleContact) noexcept
    {
        return setConnectionForType(masterSpheres, slaveLSParticles, singleContact, bondType::sphereLSParticle);
    }
    /**
     * Converts physical normal, shear, bending, and torsional stiffness values
     * into the four stored bond coefficients.
     * The derived bending coupling coefficient B3 may be negative.
     * @return False until equivalent length is set, when an input stiffness is
     * non-finite or negative, or when a derived coefficient is non-finite.
     */
    bool setStiffness(Real normalStiffness, Real shearStiffness, Real bendingStiffness, Real torsionalStiffness) noexcept
    {
        if (!hasEquivalentLength() || !isNonNegativeFinite(normalStiffness) || !isNonNegativeFinite(shearStiffness) || !isNonNegativeFinite(bendingStiffness) ||
            !isNonNegativeFinite(torsionalStiffness))
        {
            return false;
        }

        const Real coefficientB1 = normalStiffness;
        const Real coefficientB2 = shearStiffness * equivalentLength_ * equivalentLength_;
        const Real coefficientB4 = torsionalStiffness;
        const Real coefficientB3 = bendingStiffness - 0.25 * coefficientB2 - 0.5 * coefficientB4;
        if (!math::isFinite(coefficientB1) || !math::isFinite(coefficientB2) || !math::isFinite(coefficientB3) || !math::isFinite(coefficientB4))
        {
            return false;
        }

        coefficientB1_ = coefficientB1;
        coefficientB2_ = coefficientB2;
        coefficientB3_ = coefficientB3;
        coefficientB4_ = coefficientB4;
        stiffnessSet_ = true;
        return true;
    }
    bool setCrossSectionArea(Real value) noexcept { return setNonNegativeFinite(value, crossSectionArea_); }
    bool setModeICriticalEnergy(Real value) noexcept { return setNonNegativeFinite(value, modeICriticalEnergy_); }
    bool setModeIICriticalEnergy(Real value) noexcept { return setNonNegativeFinite(value, modeIICriticalEnergy_); }
    bool setMaximumEnergyReleaseRatio(Real value) noexcept { return setNonNegativeFinite(value, maximumEnergyReleaseRatio_); }
    bool setDamageFactor(Real value) noexcept
    {
        if (!math::isFinite(value) || value < 0.0 || value > 1.0)
        {
            return false;
        }
        damageFactor_ = value;
        return true;
    }
    bool setModeMixityExponent(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= 0.0)
        {
            return false;
        }
        modeMixityExponent_ = value;
        return true;
    }
    bool setDamageInitiationRatio(Real value) noexcept
    {
        if (!math::isFinite(value) || value <= 0.0 || value > 1.0)
        {
            return false;
        }
        damageInitiationRatio_ = value;
        return true;
    }

    /**
     * Restores the body-local endpoint geometry of a previously connected bond.
     * @return False when the bond has not been connected.
     */
    bool restoreReferenceState(const bondReferenceState& value) noexcept
    {
        if (type_ == bondType::undefined)
        {
            return false;
        }
        const math::Quaternion masterOrientation = math::normalizedOrIdentity(masterParticle_.orientation());
        const math::Quaternion slaveOrientation = math::normalizedOrIdentity(slaveParticle_.orientation());
        const Vec3 masterPosition = masterParticle_.position() + math::rotateUnit(masterOrientation, value.master_.position_);
        const Vec3 slavePosition = slaveParticle_.position() + math::rotateUnit(slaveOrientation, value.slave_.position_);
        const Vec3 restoredPoint = 0.5 * (masterPosition + slavePosition);
        Vec3 restoredNormal = math::normalizedOrZero(masterPosition - slavePosition);
        if (restoredNormal == Vec3::zero())
            restoredNormal = -math::rotateUnit(masterOrientation, value.master_.normal_);
        masterEndpointLocalNormal_ = value.master_.normal_;
        masterEndpointLocalTangent1_ = value.master_.tangent1_;
        masterEndpointLocalTangent2_ = value.master_.tangent2_;
        masterEndpointLocalPosition_ = value.master_.position_;
        slaveEndpointLocalNormal_ = value.slave_.normal_;
        slaveEndpointLocalTangent1_ = value.slave_.tangent1_;
        slaveEndpointLocalTangent2_ = value.slave_.tangent2_;
        slaveEndpointLocalPosition_ = value.slave_.position_;
        point_ = restoredPoint;
        normal_ = restoredNormal;
        particleCopiesUpdated_ = false;
        return true;
    }

    /** Reports whether a finite positive reference length has been accepted. */
    bool hasEquivalentLength() const noexcept { return equivalentLengthSet_; }
    /** Reports whether connection ownership, reference length, and stiffness are complete. */
    bool isValid() const noexcept { return type_ != bondType::undefined && equivalentLengthSet_ && stiffnessSet_ && masterParticleContainer_ != nullptr && slaveParticleContainer_ != nullptr; }

    /** Checks that both endpoints refer to @p particles. */
    template <class ParticleStorage> bool referencesParticleContainer(const ParticleStorage& particles) const noexcept
    {
        return masterParticleContainer_ == &particles && slaveParticleContainer_ == &particles;
    }

    /** Checks the stored master and slave container identities independently. */
    template <class MasterParticleStorage, class SlaveParticleStorage>
    bool referencesParticleContainers(const MasterParticleStorage& masterParticles, const SlaveParticleStorage& slaveParticles) const noexcept
    {
        return masterParticleContainer_ == &masterParticles && slaveParticleContainer_ == &slaveParticles;
    }

    /**
     * Refreshes both particle snapshots from their owning container and clears
     * previous resultant loads.
     * @return False before insertion or when stored indices are invalid.
     */
    template <class ParticleStorage> bool updateMasterSlaveParticle(const ParticleStorage& particles) noexcept
    {
        if (!isInContainer())
        {
            return false;
        }
        if (!setParticleCopies(particles, masterParticleIndex_, slaveParticleIndex_))
        {
            return false;
        }
        force_ = Vec3::zero();
        masterTorque_ = Vec3::zero();
        slaveTorque_ = Vec3::zero();
        particleCopiesUpdated_ = true;
        return true;
    }
    /**
     * Cross-container snapshot refresh counterpart of the single-container
     * overload.
     */
    template <class MasterParticleStorage, class SlaveParticleStorage> bool updateMasterSlaveParticle(const MasterParticleStorage& masterParticles, const SlaveParticleStorage& slaveParticles) noexcept
    {
        if (!isInContainer())
        {
            return false;
        }
        if (!setParticleCopies(masterParticles, masterParticleIndex_, slaveParticles, slaveParticleIndex_))
        {
            return false;
        }
        force_ = Vec3::zero();
        masterTorque_ = Vec3::zero();
        slaveTorque_ = Vec3::zero();
        particleCopiesUpdated_ = true;
        return true;
    }

    /**
     * Evaluates bond force, master/slave torque, elastic energy, and BK damage.
     * @return False unless the bond belongs to a container and its particle
     * snapshots were refreshed for the current assembly.
     */
    bool calculateForce() noexcept;

    struct masterParticleIndexField;
    struct slaveParticleIndexField;
    struct forceField;
    struct masterTorqueField;
    struct slaveTorqueField;
    struct pointField;
    struct normalField;
    struct equivalentLengthField;
    struct coefficientB1Field;
    struct coefficientB2Field;
    struct coefficientB3Field;
    struct coefficientB4Field;
    struct normalElasticEnergyField;
    struct shearElasticEnergyField;
    struct bendingElasticEnergyField;
    struct torsionalElasticEnergyField;
    struct masterEndpointLocalNormalField;
    struct masterEndpointLocalTangent1Field;
    struct masterEndpointLocalTangent2Field;
    struct masterEndpointLocalPositionField;
    struct slaveEndpointLocalNormalField;
    struct slaveEndpointLocalTangent1Field;
    struct slaveEndpointLocalTangent2Field;
    struct slaveEndpointLocalPositionField;
    struct crossSectionAreaField;
    struct damageFactorField;
    struct modeICriticalEnergyField;
    struct modeIICriticalEnergyField;
    struct modeMixityExponentField;
    struct damageInitiationRatioField;
    struct maximumEnergyReleaseRatioField;

private:
    /** Clears every connection and stiffness value derived from reference length. */
    void invalidateLengthDependentValues() noexcept
    {
        particleCopiesUpdated_ = false;
        stiffnessSet_ = false;
        masterParticleContainer_ = nullptr;
        slaveParticleContainer_ = nullptr;
        type_ = bondType::undefined;
        masterParticleIndex_ = -1;
        slaveParticleIndex_ = -1;
        coefficientB1_ = 0.0;
        coefficientB2_ = 0.0;
        coefficientB3_ = 0.0;
        coefficientB4_ = 0.0;
    }

    /** Validates two same-container indices and refreshes their body snapshots. */
    template <class ParticleStorage> bool setParticleCopies(const ParticleStorage& particles, int masterIndex, int slaveIndex) noexcept
    {
        particleCopiesUpdated_ = false;
        const int particleCount = static_cast<int>(particles.hostSize());
        if (masterIndex < 0 || masterIndex >= particleCount || slaveIndex < 0 || slaveIndex >= particleCount || masterIndex == slaveIndex)
        {
            return false;
        }

        masterParticle_ = particles.host()[masterIndex];
        slaveParticle_ = particles.host()[slaveIndex];
        masterParticleIndex_ = masterIndex;
        slaveParticleIndex_ = slaveIndex;
        return true;
    }

    /** Validates cross-container indices and refreshes both body snapshots. */
    template <class MasterParticleStorage, class SlaveParticleStorage>
    bool setParticleCopies(const MasterParticleStorage& masterParticles, int masterIndex, const SlaveParticleStorage& slaveParticles, int slaveIndex) noexcept
    {
        particleCopiesUpdated_ = false;
        const int masterParticleCount = static_cast<int>(masterParticles.hostSize());
        const int slaveParticleCount = static_cast<int>(slaveParticles.hostSize());
        if (masterIndex < 0 || masterIndex >= masterParticleCount || slaveIndex < 0 || slaveIndex >= slaveParticleCount)
        {
            return false;
        }

        masterParticle_ = masterParticles.host()[masterIndex];
        slaveParticle_ = slaveParticles.host()[slaveIndex];
        masterParticleIndex_ = masterIndex;
        slaveParticleIndex_ = slaveIndex;
        return true;
    }

    /** Establishes a typed same-container connection after validating both endpoints. */
    template <class ParticleStorage> bool setConnectionForType(const ParticleStorage& particles, int masterParticleIndex, int slaveParticleIndex, const Vec3& normal, bondType type) noexcept
    {
        type_ = bondType::undefined;
        masterParticleContainer_ = nullptr;
        slaveParticleContainer_ = nullptr;
        if (!setParticleCopies(particles, masterParticleIndex, slaveParticleIndex) || !setReferenceGeometry(0.5 * (masterParticle_.position() + slaveParticle_.position()), normal))
        {
            particleCopiesUpdated_ = false;
            return false;
        }
        masterParticleContainer_ = &particles;
        slaveParticleContainer_ = &particles;
        type_ = type;
        return true;
    }

    /** Establishes a typed cross-container connection after validating both endpoints. */
    template <class MasterParticleStorage, class SlaveParticleStorage>
    bool setConnectionForType(const MasterParticleStorage& masterParticles,
                              int masterParticleIndex,
                              const SlaveParticleStorage& slaveParticles,
                              int slaveParticleIndex,
                              const Vec3& normal,
                              bondType type) noexcept
    {
        type_ = bondType::undefined;
        masterParticleContainer_ = nullptr;
        slaveParticleContainer_ = nullptr;
        if (!setParticleCopies(masterParticles, masterParticleIndex, slaveParticles, slaveParticleIndex) ||
            !setReferenceGeometry(0.5 * (masterParticle_.position() + slaveParticle_.position()), normal))
        {
            particleCopiesUpdated_ = false;
            return false;
        }
        masterParticleContainer_ = &masterParticles;
        slaveParticleContainer_ = &slaveParticles;
        type_ = type;
        return true;
    }

    /** Establishes a typed same-container connection from existing contact geometry. */
    template <class ParticleStorage> bool setConnectionForType(const ParticleStorage& particles, const contact& singleContact, bondType type) noexcept
    {
        type_ = bondType::undefined;
        masterParticleContainer_ = nullptr;
        slaveParticleContainer_ = nullptr;
        if (!setParticleCopies(particles, singleContact.masterParticleIndex(), singleContact.slaveParticleIndex()) || !setReferenceGeometry(singleContact.point(), singleContact.normal()))
        {
            particleCopiesUpdated_ = false;
            return false;
        }
        masterParticleContainer_ = &particles;
        slaveParticleContainer_ = &particles;
        type_ = type;
        return true;
    }

    /** Establishes a typed cross-container connection from existing contact geometry. */
    template <class MasterParticleStorage, class SlaveParticleStorage>
    bool setConnectionForType(const MasterParticleStorage& masterParticles, const SlaveParticleStorage& slaveParticles, const contact& singleContact, bondType type) noexcept
    {
        type_ = bondType::undefined;
        masterParticleContainer_ = nullptr;
        slaveParticleContainer_ = nullptr;
        if (!setParticleCopies(masterParticles, singleContact.masterParticleIndex(), slaveParticles, singleContact.slaveParticleIndex()) ||
            !setReferenceGeometry(singleContact.point(), singleContact.normal()))
        {
            particleCopiesUpdated_ = false;
            return false;
        }
        masterParticleContainer_ = &masterParticles;
        slaveParticleContainer_ = &slaveParticles;
        type_ = type;
        return true;
    }

    /** Tests the common admissible domain of fracture and stiffness data. */
    static bool isNonNegativeFinite(Real value) noexcept { return math::isFinite(value) && value >= 0.0; }

    /** Assigns a finite non-negative value without modifying the destination on failure. */
    static bool setNonNegativeFinite(Real value, Real& destination) noexcept
    {
        if (!isNonNegativeFinite(value))
        {
            return false;
        }
        destination = value;
        return true;
    }

    /** Builds orthonormal endpoint frames and stores them in body coordinates. */
    bool setReferenceGeometry(const Vec3& point, const Vec3& normal) noexcept
    {
        const Vec3 unitNormal = math::normalizedOrZero(normal);
        if (unitNormal == Vec3::zero() || !hasEquivalentLength())
        {
            return false;
        }

        point_ = point;
        normal_ = unitNormal;

        const Vec3 masterNormal = -normal_;
        const Vec3 reference = masterNormal.x > -0.9 && masterNormal.x < 0.9 ? Vec3::unitX() : Vec3::unitY();
        const Vec3 tangent1 = math::normalizedOrZero(math::cross(reference, masterNormal));
        const Vec3 tangent2 = math::normalizedOrZero(math::cross(masterNormal, tangent1));
        const Vec3 masterPosition = point_ + 0.5 * equivalentLength_ * normal_;
        const Vec3 slavePosition = point_ - 0.5 * equivalentLength_ * normal_;
        const math::Quaternion masterOrientation = math::normalizedOrIdentity(masterParticle_.orientation());
        const math::Quaternion slaveOrientation = math::normalizedOrIdentity(slaveParticle_.orientation());
        masterEndpointLocalNormal_ = math::inverseRotateUnit(masterOrientation, masterNormal);
        masterEndpointLocalTangent1_ = math::inverseRotateUnit(masterOrientation, tangent1);
        masterEndpointLocalTangent2_ = math::inverseRotateUnit(masterOrientation, tangent2);
        masterEndpointLocalPosition_ = math::inverseRotateUnit(masterOrientation, masterPosition - masterParticle_.position());
        slaveEndpointLocalNormal_ = math::inverseRotateUnit(slaveOrientation, normal_);
        slaveEndpointLocalTangent1_ = math::inverseRotateUnit(slaveOrientation, tangent1);
        slaveEndpointLocalTangent2_ = math::inverseRotateUnit(slaveOrientation, tangent2);
        slaveEndpointLocalPosition_ = math::inverseRotateUnit(slaveOrientation, slavePosition - slaveParticle_.position());
        return true;
    }

    rigidBody masterParticle_;                     ///< Host snapshot of the master body.
    rigidBody slaveParticle_;                      ///< Host snapshot of the slave body.
    bool particleCopiesUpdated_{false};            ///< Force-evaluation readiness flag.
    bool equivalentLengthSet_{false};              ///< Whether reference length is usable.
    bool stiffnessSet_{false};                     ///< Whether all four coefficients are current.
    const void* masterParticleContainer_{nullptr}; ///< Owning master container identity.
    const void* slaveParticleContainer_{nullptr};  ///< Owning slave container identity.
    bondType type_{bondType::undefined};           ///< Stored particle-type pairing.

    int masterParticleIndex_{-1};                     ///< Master particle-container index.
    int slaveParticleIndex_{-1};                      ///< Slave particle-container index.
    Vec3 force_{Vec3::zero()};                        ///< Force applied to the master.
    Vec3 masterTorque_{Vec3::zero()};                 ///< Torque applied to the master.
    Vec3 slaveTorque_{Vec3::zero()};                  ///< Torque applied to the slave.
    Vec3 point_{Vec3::zero()};                        ///< Reference bond-center position.
    Vec3 normal_{Vec3::zero()};                       ///< Reference normal from slave to master.
    Real equivalentLength_{0.0};                      ///< Reference endpoint separation.
    Real coefficientB1_{0.0};                         ///< Normal stiffness coefficient.
    Real coefficientB2_{0.0};                         ///< Length-scaled shear coefficient.
    Real coefficientB3_{0.0};                         ///< Bending coupling coefficient.
    Real coefficientB4_{0.0};                         ///< Torsional stiffness coefficient.
    Real normalElasticEnergy_{0.0};                   ///< Stored normal elastic energy.
    Real shearElasticEnergy_{0.0};                    ///< Stored shear elastic energy.
    Real bendingElasticEnergy_{0.0};                  ///< Stored bending elastic energy.
    Real torsionalElasticEnergy_{0.0};                ///< Stored torsional elastic energy.
    Vec3 masterEndpointLocalNormal_{Vec3::unitX()};   ///< Master-frame endpoint normal.
    Vec3 masterEndpointLocalTangent1_{Vec3::unitY()}; ///< Master-frame first tangent.
    Vec3 masterEndpointLocalTangent2_{Vec3::unitZ()}; ///< Master-frame second tangent.
    Vec3 masterEndpointLocalPosition_{Vec3::zero()};  ///< Master-frame endpoint position.
    Vec3 slaveEndpointLocalNormal_{Vec3::unitX()};    ///< Slave-frame endpoint normal.
    Vec3 slaveEndpointLocalTangent1_{Vec3::unitY()};  ///< Slave-frame first tangent.
    Vec3 slaveEndpointLocalTangent2_{Vec3::unitZ()};  ///< Slave-frame second tangent.
    Vec3 slaveEndpointLocalPosition_{Vec3::zero()};   ///< Slave-frame endpoint position.
    Real crossSectionArea_{0.0};                      ///< Nominal cross-sectional area; zero disables fracture.
    Real damageFactor_{0.0};                          ///< Current damage in [0,1].
    Real modeICriticalEnergy_{0.0};                   ///< Mode-I critical energy release rate.
    Real modeIICriticalEnergy_{0.0};                  ///< Mode-II critical energy release rate.
    Real modeMixityExponent_{1.75};                   ///< BK mixed-mode exponent.
    Real damageInitiationRatio_{0.9};                 ///< Fraction of critical energy at damage onset.
    Real maximumEnergyReleaseRatio_{0.0};             ///< Largest normalized release reached.

public:
    struct masterParticleIndexField {
        inline static constexpr auto member = &bond::masterParticleIndex_;
    };
    struct slaveParticleIndexField {
        inline static constexpr auto member = &bond::slaveParticleIndex_;
    };
    struct forceField {
        inline static constexpr auto member = &bond::force_;
    };
    struct masterTorqueField {
        inline static constexpr auto member = &bond::masterTorque_;
    };
    struct slaveTorqueField {
        inline static constexpr auto member = &bond::slaveTorque_;
    };
    struct pointField {
        inline static constexpr auto member = &bond::point_;
    };
    struct normalField {
        inline static constexpr auto member = &bond::normal_;
    };
    struct equivalentLengthField {
        inline static constexpr auto member = &bond::equivalentLength_;
    };
    struct coefficientB1Field {
        inline static constexpr auto member = &bond::coefficientB1_;
    };
    struct coefficientB2Field {
        inline static constexpr auto member = &bond::coefficientB2_;
    };
    struct coefficientB3Field {
        inline static constexpr auto member = &bond::coefficientB3_;
    };
    struct coefficientB4Field {
        inline static constexpr auto member = &bond::coefficientB4_;
    };
    struct normalElasticEnergyField {
        inline static constexpr auto member = &bond::normalElasticEnergy_;
    };
    struct shearElasticEnergyField {
        inline static constexpr auto member = &bond::shearElasticEnergy_;
    };
    struct bendingElasticEnergyField {
        inline static constexpr auto member = &bond::bendingElasticEnergy_;
    };
    struct torsionalElasticEnergyField {
        inline static constexpr auto member = &bond::torsionalElasticEnergy_;
    };
    struct masterEndpointLocalNormalField {
        inline static constexpr auto member = &bond::masterEndpointLocalNormal_;
    };
    struct masterEndpointLocalTangent1Field {
        inline static constexpr auto member = &bond::masterEndpointLocalTangent1_;
    };
    struct masterEndpointLocalTangent2Field {
        inline static constexpr auto member = &bond::masterEndpointLocalTangent2_;
    };
    struct masterEndpointLocalPositionField {
        inline static constexpr auto member = &bond::masterEndpointLocalPosition_;
    };
    struct slaveEndpointLocalNormalField {
        inline static constexpr auto member = &bond::slaveEndpointLocalNormal_;
    };
    struct slaveEndpointLocalTangent1Field {
        inline static constexpr auto member = &bond::slaveEndpointLocalTangent1_;
    };
    struct slaveEndpointLocalTangent2Field {
        inline static constexpr auto member = &bond::slaveEndpointLocalTangent2_;
    };
    struct slaveEndpointLocalPositionField {
        inline static constexpr auto member = &bond::slaveEndpointLocalPosition_;
    };
    struct crossSectionAreaField {
        inline static constexpr auto member = &bond::crossSectionArea_;
    };
    struct damageFactorField {
        inline static constexpr auto member = &bond::damageFactor_;
    };
    struct modeICriticalEnergyField {
        inline static constexpr auto member = &bond::modeICriticalEnergy_;
    };
    struct modeIICriticalEnergyField {
        inline static constexpr auto member = &bond::modeIICriticalEnergy_;
    };
    struct modeMixityExponentField {
        inline static constexpr auto member = &bond::modeMixityExponent_;
    };
    struct damageInitiationRatioField {
        inline static constexpr auto member = &bond::damageInitiationRatio_;
    };
    struct maximumEnergyReleaseRatioField {
        inline static constexpr auto member = &bond::maximumEnergyReleaseRatio_;
    };

    using DeviceLayout = deviceLayout<masterParticleIndexField,
                                      slaveParticleIndexField,
                                      forceField,
                                      masterTorqueField,
                                      slaveTorqueField,
                                      pointField,
                                      normalField,
                                      equivalentLengthField,
                                      coefficientB1Field,
                                      coefficientB2Field,
                                      coefficientB3Field,
                                      coefficientB4Field,
                                      normalElasticEnergyField,
                                      shearElasticEnergyField,
                                      bendingElasticEnergyField,
                                      torsionalElasticEnergyField,
                                      masterEndpointLocalNormalField,
                                      masterEndpointLocalTangent1Field,
                                      masterEndpointLocalTangent2Field,
                                      masterEndpointLocalPositionField,
                                      slaveEndpointLocalNormalField,
                                      slaveEndpointLocalTangent1Field,
                                      slaveEndpointLocalTangent2Field,
                                      slaveEndpointLocalPositionField,
                                      crossSectionAreaField,
                                      damageFactorField,
                                      modeICriticalEnergyField,
                                      modeIICriticalEnergyField,
                                      modeMixityExponentField,
                                      damageInitiationRatioField,
                                      maximumEnergyReleaseRatioField>;
};

using bondContainer = hostAoSDeviceSoA<bond, bond::DeviceLayout>;

} // namespace fundem
