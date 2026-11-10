/**
 * @file periodicShear.h
 * @brief CPU-only ON-DEM bumpy-wall shear validation without changing default DEM laws.
 */
#pragma once

#include "periodicContactForce.h"
#include "solver/SphereDEM.h"

#include <vector>

namespace fundem::validation
{

using math::Real;
using math::Vec3;

/** One source-file sphere; insertion order and source IDs are retained exactly. */
struct initialSphere {
    int sourceId_{0};             ///< Source particle identifier, independent of the container index.
    Vec3 position_{Vec3::zero()}; ///< Initial center in the benchmark coordinates.
    Real radius_{0.5};            ///< Monodisperse sphere radius.
    int wall_{0};                 ///< -1: lower wall; +1: upper wall; 0: finite-mass grain.
};

/**
 * Direct external contact coefficients and prescribed wall speed.
 * The dashpots are force/relative-velocity coefficients, not restitution ratios.
 */
struct shearParameters : contactForceParameters {
    Real wallSpeed_{1.0}; ///< Wall speed magnitude; upper +X, lower -X.
};

/**
 * Dedicated CPU shear benchmark with X/Y periodicity and nonperiodic Z.
 *
 * Reuses FunDEM particle storage, rigid-body integration and output. Only this
 * validation uses its direct damping coefficients and periodic contact search.
 * Grain density is one; fixed wall grains retain their source order and follow
 * prescribed opposite X velocities. No side walls replace periodic boundaries.
 *
 * The Coulomb cap applies to the complete tangential spring-plus-dashpot force,
 * and the spring is returned consistently after sliding. The benchmark does not
 * establish that every participating code uses the same sliding convention.
 */
class periodicShear : public SphereDEM
{
public:
    /** Validates and inserts an ordered monodisperse source configuration. */
    periodicShear(const std::vector<initialSphere>& initial, const Vec3& minimum, const Vec3& maximum, shearParameters parameters = {});

    using SphereDEM::initialize;

    const std::vector<shearContact>& contacts() const noexcept { return contactForce_.contacts(); }
    const std::vector<int>& sourceIds() const noexcept { return sourceIds_; }
    const shearParameters& parameters() const noexcept { return parameters_; }

    /** Maps X/Y to the half-open periodic cell and leaves Z unchanged. */
    Vec3 wrapPosition(Vec3 position) const noexcept;
    /** Returns the nearest X/Y image displacement, leaving Z unchanged. */
    Vec3 minimumImage(Vec3 displacement) const noexcept;

protected:
    /** Initializes forces at the current configuration without advancing history. */
    void initialize(cpu::mode) override;
    /** Rejects a backend not implemented by this validation. */
    void initialize(gpu::mode, cudaStream_t stream) override;
    /** Rejects a backend not implemented by this validation. */
    void initialize(hybrid::mode, cudaStream_t stream) override;
    /** Advances one velocity-Verlet step using periodic contacts and prescribed walls. */
    void advanceCPU(Real timeStep) override;
    /** Rejects unsupported device advancement. */
    void advanceGPU(Real timeStep) override;
    /** Rejects unsupported mixed advancement. */
    void advanceHybrid(Real timeStep) override;
    /** Restricts this dedicated contact implementation to the CPU backend. */
    void validateExecutionMode(executionMode value) const override;
    /** Rejects unsupported changes to the benchmark population, gravity, or domain. */
    void validateSystemConfiguration() const override;
    /** Includes the validation's direct normal and tangential spring energies. */
    energyRecord systemEnergy() const override;
    /** Adds the standalone periodic contact force through the existing external-force hook. */
    void addSphereExternalForceAndTorque(particleContainer::host_container_type& particles) override;

private:
    /** Enforces wall kinematics and wraps all current particle centers. */
    void prepareKinematics();
    /** Clears loads once and invokes the additive external-force hook. */
    void assembleContacts(Real historyTimeStep);

    shearParameters parameters_;        ///< Benchmark-only force and wall-speed settings.
    periodicContactForce contactForce_; ///< Additive periodic contact law, search and spring histories.
    std::vector<int> sourceIds_;        ///< Original IDs in unchanged input order.
    std::vector<int> wallKinds_;        ///< Fixed-wall role aligned with stable sphere indices.
    std::vector<Real> wallZ_;           ///< Prescribed wall Z coordinates.
    Real historyTimeStep_{0.0};         ///< Current external-force history increment, zero at initialization.
};

} // namespace fundem::validation
