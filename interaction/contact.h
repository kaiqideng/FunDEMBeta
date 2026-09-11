/**
 * @file contact.h
 * @brief Defines contact state, material response data, loads, and device storage.
 */
#pragma once

#include "data/ContainerObject.h"
#include "data/HostAoSDeviceSoA.h"
#include "math/Vector3.h"
#include "particle/particle.h"

namespace fundem
{

/**
 * One ordered contact with host particle snapshots, geometric response state,
 * accumulated spring history, resultant force, and torque.
 */
class contact : public containerObject
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    contact() = default;
    ~contact() = default;
    contact(const contact&) = default;
    contact& operator=(const contact&) = default;
    contact(contact&&) noexcept = default;
    contact& operator=(contact&&) noexcept = default;

    const particle& masterParticle() const noexcept { return masterParticle_; }
    const particle& slaveParticle() const noexcept { return slaveParticle_; }
    const Vec3& force() const noexcept { return force_; }
    const Vec3& torque() const noexcept { return torque_; }
    int masterParticleIndex() const noexcept { return masterParticleIndex_; }
    int slaveParticleIndex() const noexcept { return slaveParticleIndex_; }
    int particleSurfaceNodeIndex() const noexcept { return particleSurfaceNodeIndex_; }
    const Vec3& point() const noexcept { return point_; }
    const Vec3& normal() const noexcept { return normal_; }
    Real overlap() const noexcept { return overlap_; }
    Real area() const noexcept { return area_; }
    Real effectiveMass() const noexcept { return effectiveMass_; }
    Real effectiveRadius() const noexcept { return effectiveRadius_; }
    Real normalForceMagnitude() const noexcept { return normalForceMagnitude_; }
    Real normalElasticEnergy() const noexcept { return normalElasticEnergy_; }
    Real slidingElasticEnergy() const noexcept { return slidingElasticEnergy_; }
    Real rollingElasticEnergy() const noexcept { return rollingElasticEnergy_; }
    Real torsionalElasticEnergy() const noexcept { return torsionalElasticEnergy_; }
    const Vec3& slidingSpringDeformation() const noexcept { return slidingSpringDeformation_; }
    const Vec3& rollingSpringDeformation() const noexcept { return rollingSpringDeformation_; }
    const Vec3& torsionalSpringDeformation() const noexcept { return torsionalSpringDeformation_; }

    /**
     * Refreshes master/slave snapshots from one particle container.
     * The update also clears the force-evaluation readiness flag until both
     * indices have been validated.
     * @return True when distinct valid indices were copied.
     */
    template <class ParticleStorage> bool setMasterSlaveParticle(const ParticleStorage& particles, int masterIndex, int slaveIndex) noexcept
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
        particleCopiesUpdated_ = true;
        return true;
    }

    /**
     * Refreshes master/slave snapshots from different particle containers.
     * Cross-type contacts permit identical numeric indices because the owners
     * differ.
     * @return True when both indices were valid.
     */
    template <class MasterParticleStorage, class SlaveParticleStorage>
    bool setMasterSlaveParticle(const MasterParticleStorage& masterParticles, int masterIndex, const SlaveParticleStorage& slaveParticles, int slaveIndex) noexcept
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
        particleCopiesUpdated_ = true;
        return true;
    }
    void setParticleSurfaceNodeIndex(int value) noexcept { particleSurfaceNodeIndex_ = value; }
    void setPoint(const Vec3& value) noexcept { point_ = value; }
    void setNormal(const Vec3& value) noexcept { normal_ = value; }
    void setOverlap(Real value) noexcept { overlap_ = value; }
    void setArea(Real value) noexcept { area_ = value; }
    void setEffectiveMass(Real value) noexcept { effectiveMass_ = value; }
    void setEffectiveRadius(Real value) noexcept { effectiveRadius_ = value; }
    void setSlidingSpringDeformation(const Vec3& value) noexcept { slidingSpringDeformation_ = value; }
    void setRollingSpringDeformation(const Vec3& value) noexcept { rollingSpringDeformation_ = value; }
    void setTorsionalSpringDeformation(const Vec3& value) noexcept { torsionalSpringDeformation_ = value; }

    /**
     * Evaluates the material contact law from the current snapshots and history.
     * @param timeStep Positive integration step used by damping and spring updates.
     * @return False unless this object belongs to a container and particle
     * snapshots were refreshed for the current force assembly.
     */
    bool calculateForce(Real timeStep) noexcept;

    struct masterParticleIndexField;
    struct slaveParticleIndexField;
    struct particleSurfaceNodeIndexField;
    struct forceField;
    struct torqueField;
    struct pointField;
    struct normalField;
    struct overlapField;
    struct areaField;
    struct effectiveMassField;
    struct effectiveRadiusField;
    struct normalForceMagnitudeField;
    struct normalElasticEnergyField;
    struct slidingElasticEnergyField;
    struct rollingElasticEnergyField;
    struct torsionalElasticEnergyField;
    struct slidingSpringDeformationField;
    struct rollingSpringDeformationField;
    struct torsionalSpringDeformationField;

private:
    particle masterParticle_;           ///< Host snapshot of the ordered master.
    particle slaveParticle_;            ///< Host snapshot of the ordered slave.
    bool particleCopiesUpdated_{false}; ///< Force-evaluation readiness flag.
    int particleSurfaceNodeIndex_{-1};  ///< Owner-local LS surface-node index, or `-1`.

    int masterParticleIndex_{-1};                   ///< Master particle-container index.
    int slaveParticleIndex_{-1};                    ///< Slave particle-container index.
    Vec3 force_{Vec3::zero()};                      ///< Force applied to the master.
    Vec3 torque_{Vec3::zero()};                     ///< Contact torque applied to the master.
    Vec3 point_{Vec3::zero()};                      ///< World contact point.
    Vec3 normal_{Vec3::zero()};                     ///< Unit normal directed from slave to master.
    Real overlap_{0.0};                             ///< Positive normal penetration.
    Real area_{0.0};                                ///< Contact or LS patch area.
    Real effectiveMass_{0.0};                       ///< Effective mass used by damping.
    Real effectiveRadius_{0.0};                     ///< Effective radius used by moment laws.
    Real normalForceMagnitude_{0.0};                ///< Absolute normal force magnitude.
    Real normalElasticEnergy_{0.0};                 ///< Stored normal elastic energy.
    Real slidingElasticEnergy_{0.0};                ///< Stored sliding elastic energy.
    Real rollingElasticEnergy_{0.0};                ///< Stored rolling elastic energy.
    Real torsionalElasticEnergy_{0.0};              ///< Stored torsional elastic energy.
    Vec3 slidingSpringDeformation_{Vec3::zero()};   ///< Tangential spring displacement.
    Vec3 rollingSpringDeformation_{Vec3::zero()};   ///< Rolling displacement history in length units, scaled by the effective radius.
    Vec3 torsionalSpringDeformation_{Vec3::zero()}; ///< Torsional displacement history in length units, scaled by twice the effective radius.

public:
    struct masterParticleIndexField {
        inline static constexpr auto member = &contact::masterParticleIndex_;
    };
    struct slaveParticleIndexField {
        inline static constexpr auto member = &contact::slaveParticleIndex_;
    };
    struct particleSurfaceNodeIndexField {
        inline static constexpr auto member = &contact::particleSurfaceNodeIndex_;
    };
    struct forceField {
        inline static constexpr auto member = &contact::force_;
    };
    struct torqueField {
        inline static constexpr auto member = &contact::torque_;
    };
    struct pointField {
        inline static constexpr auto member = &contact::point_;
    };
    struct normalField {
        inline static constexpr auto member = &contact::normal_;
    };
    struct overlapField {
        inline static constexpr auto member = &contact::overlap_;
    };
    struct areaField {
        inline static constexpr auto member = &contact::area_;
    };
    struct effectiveMassField {
        inline static constexpr auto member = &contact::effectiveMass_;
    };
    struct effectiveRadiusField {
        inline static constexpr auto member = &contact::effectiveRadius_;
    };
    struct normalForceMagnitudeField {
        inline static constexpr auto member = &contact::normalForceMagnitude_;
    };
    struct normalElasticEnergyField {
        inline static constexpr auto member = &contact::normalElasticEnergy_;
    };
    struct slidingElasticEnergyField {
        inline static constexpr auto member = &contact::slidingElasticEnergy_;
    };
    struct rollingElasticEnergyField {
        inline static constexpr auto member = &contact::rollingElasticEnergy_;
    };
    struct torsionalElasticEnergyField {
        inline static constexpr auto member = &contact::torsionalElasticEnergy_;
    };
    struct slidingSpringDeformationField {
        inline static constexpr auto member = &contact::slidingSpringDeformation_;
    };
    struct rollingSpringDeformationField {
        inline static constexpr auto member = &contact::rollingSpringDeformation_;
    };
    struct torsionalSpringDeformationField {
        inline static constexpr auto member = &contact::torsionalSpringDeformation_;
    };

    using DeviceLayout = deviceLayout<masterParticleIndexField,
                                      slaveParticleIndexField,
                                      particleSurfaceNodeIndexField,
                                      forceField,
                                      torqueField,
                                      pointField,
                                      normalField,
                                      overlapField,
                                      areaField,
                                      effectiveMassField,
                                      effectiveRadiusField,
                                      normalForceMagnitudeField,
                                      normalElasticEnergyField,
                                      slidingElasticEnergyField,
                                      rollingElasticEnergyField,
                                      torsionalElasticEnergyField,
                                      slidingSpringDeformationField,
                                      rollingSpringDeformationField,
                                      torsionalSpringDeformationField>;
};

using contactContainer = hostAoSDeviceSoA<contact, contact::DeviceLayout>;

/** Minimal previous-step state required to continue a persistent contact. */
class contactHistory
{
public:
    using Vec3 = math::Vec3;

    contactHistory() = default;
    /** Captures the slave identity and three persistent spring states. */
    contactHistory(int slaveParticleIndex, const Vec3& slidingSpringDeformation, const Vec3& rollingSpringDeformation, const Vec3& torsionalSpringDeformation) noexcept
        : slaveParticleIndex_(slaveParticleIndex), slidingSpringDeformation_(slidingSpringDeformation), rollingSpringDeformation_(rollingSpringDeformation),
          torsionalSpringDeformation_(torsionalSpringDeformation)
    {}

private:
    int slaveParticleIndex_{-1};                    ///< Slave identifier used for history matching.
    Vec3 slidingSpringDeformation_{Vec3::zero()};   ///< Previous sliding spring state.
    Vec3 rollingSpringDeformation_{Vec3::zero()};   ///< Previous rolling spring state.
    Vec3 torsionalSpringDeformation_{Vec3::zero()}; ///< Previous torsional spring state.

public:
    /** Mutable device view of compact contact-history arrays. */
    struct device_type {
        int* slaveParticleIndex_{nullptr};          ///< Previous slave indices.
        Vec3* slidingSpringDeformation_{nullptr};   ///< Previous sliding states.
        Vec3* rollingSpringDeformation_{nullptr};   ///< Previous rolling states.
        Vec3* torsionalSpringDeformation_{nullptr}; ///< Previous torsional states.
    };

    struct slaveParticleIndexField {
        inline static constexpr auto member = &contactHistory::slaveParticleIndex_;
    };
    struct slidingSpringDeformationField {
        inline static constexpr auto member = &contactHistory::slidingSpringDeformation_;
    };
    struct rollingSpringDeformationField {
        inline static constexpr auto member = &contactHistory::rollingSpringDeformation_;
    };
    struct torsionalSpringDeformationField {
        inline static constexpr auto member = &contactHistory::torsionalSpringDeformation_;
    };

    using DeviceLayout = deviceLayout<slaveParticleIndexField, slidingSpringDeformationField, rollingSpringDeformationField, torsionalSpringDeformationField>;
};

using contactHistoryContainer = hostAoSDeviceSoA<contactHistory, contactHistory::DeviceLayout>;

} // namespace fundem
