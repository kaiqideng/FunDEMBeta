/**
 * @file SPHDEM.h
 * @brief Declares CPU, GPU, and Hybrid WCSPH coupling to level-set solids.
 */
#pragma once

#include "particle/SPHParticle.h"
#include "particle/virtualParticle.h"
#include "solver/LSDEM.h"
#include "particle/SPHJet.h"
#include "interaction/SPHNeighborhood.h"
#include "interaction/SPHInteractionTypes.h"
#include "interaction/virtualParticleCoupling.h"
#include "execution/sphJetFunctions.h"

#include <memory>
#include <vector>

namespace fundem
{

namespace cpu
{
class SPHInteraction;
}

/**
 * Multirate WCSPH solver coupled to rigid level-set solids. SPH may run on CPU
 * or GPU; Hybrid mode keeps LS solids on CPU while SPH remains on the device.
 * Spherical DEM particles are intentionally unavailable.
 */
class SPHDEM : public LSDEM
{
public:
    /** Exposes the common lazy initialization entry point. */
    void initialize() { solver::initialize(); }

    /** Creates a CPU SPH/LS solver targeting @p device for later GPU use. */
    explicit SPHDEM(int device = 0);
    /** Creates an SPH/LS solver in the requested execution mode. */
    SPHDEM(executionMode mode, int device = 0);
    ~SPHDEM() override;

    /** Returns integration storage; use observeCurrentState() for a common-time view between acoustic endpoints. */
    const SPHParticleContainer& SPHParticles() const noexcept { return SPHParticles_; }
    const virtualParticleContainer& virtualParticles() const noexcept { return virtualParticles_; }

    /** Appends one valid SPH particle after uniform fluid properties are set. */
    int addSPHParticle(const SPHParticle& value);
    /** Appends a batch and returns the first stable index. */
    int addSPHParticles(std::vector<SPHParticle> values);
    /** Appends a cell-centered rectangular lattice of SPH particles. */
    int addSPHBlock(const Vec3& blockMinimum, const int3& particleCount, const Vec3& initialVelocity = Vec3::zero());
    /** Adds a finite downstream cylindrical particle column without an inlet constraint. */
    int addSPHJet(const Vec3& inletCenter, const Vec3& velocity, Real radius, Real length);
    /** Adds a short inlet buffer which emits append-only particle layers until its finite supply is exhausted. */
    int addSPHJet(const SPHJet& value);
    /** Returns this jet's currently created global indices, in axial-layer order; @p firstParticleIndex is the handle returned by addSPHJet. */
    const std::vector<int>& SPHJetParticleIndices(int firstParticleIndex) const;
    /** Reports whether every finite-duration jet has ended at the solver's visible time. */
    bool SPHJetsCompleted() const noexcept { return SPHJets_.completedAt(time()); }
    Real SPHMaximumVelocity() const noexcept { return SPHProperties_.maximumVelocity_; }
    Real SPHSoundSpeed() const noexcept { return SPHProperties_.soundSpeed_; }
    Real SPHSpacing() const noexcept { return SPHProperties_.spacing_; }
    Real SPHSmoothingLength() const noexcept { return SPHProperties_.smoothingLength_; }
    Real SPHReferenceDensity() const noexcept { return SPHProperties_.referenceDensity_; }
    Real SPHDynamicViscosity() const noexcept { return SPHProperties_.dynamicViscosity_; }
    Real SPHTimeStepLimit() const noexcept { return SPHStepState_.acoustic_.limit_; }
    Real SPHTimeStep() const noexcept { return SPHStepState_.acoustic_.timeStep_; }
    Real SPHAdvectionTimeStepLimit() const noexcept { return SPHStepState_.advection_.limit_; }
    Real SPHAdvectionTimeStep() const noexcept { return SPHStepState_.advection_.timeStep_; }
    int SPHStepInterval() const noexcept { return SPHStepState_.acoustic_.interval_; }
    int SPHAdvectionStepInterval() const noexcept { return SPHStepState_.advection_.interval_; }

    /** Sets uniform SPH discretization and fluid properties shared by all SPH particles. */
    void setSPHProperties(Real spacing, Real smoothingLength, Real referenceDensity, Real dynamicViscosity);
    /** Sets optional translational images for native fluid and virtual-boundary interactions. */
    void setSPHPeriodicDomain(const math::periodicDomain& value);
    const math::periodicDomain& SPHPeriodicDomain() const noexcept { return SPHPeriodicDomain_; }
    /** Sets the configured SPH velocity scale and resets the sound speed to ten times this value. */
    void setSPHMaximumVelocity(Real value);
    void setSPHSoundSpeed(Real value);
    void setSPHParticleVTUFields(std::vector<SPHParticleVTUField> fields);

protected:
    /** Mutable fluid storage for derived frontends restoring a saved host state. */
    SPHParticleContainer& mutableSPHParticles() noexcept { return SPHParticles_; }
    /** Declares scheduled boundary motion that may start after initialization. */
    void setSPHBoundaryMotionExpected(bool value) noexcept { SPHBoundaryMotionExpected_ = value; }
    /** Accepts CPU, GPU, and SPH-GPU/LS-CPU Hybrid modes. */
    void validateExecutionMode(executionMode value) const override;
    /** Updates fluid and LS storage selection after a mode change. */
    void executionModeChanged(executionMode value) override;
    /** Verifies uniform SPH properties and a usable fluid population. */
    void validateSystemConfiguration() const override;
    /** Initializes host SPH and LS state. */
    void initialize(cpu::mode) override;
    /** Initializes device SPH and LS state. */
    void initialize(gpu::mode, cudaStream_t stream) override;
    /** Initializes device SPH state and host LS boundary state. */
    void initialize(hybrid::mode, cudaStream_t stream) override;
    /** Advances one all-host multirate SPH/LS step. */
    void advanceCPU(Real timeStep) override;
    /** Prepares prescribed motion before boundary sampling and the first solid half kick. */
    virtual void beforeLSParticleDrift(cpu::mode) {}
    /** Prepares device prescribed motion before the leading fluid interval is evaluated. */
    virtual void beforeLSParticleDrift(gpu::mode, cudaStream_t) {}
    /** Advances one all-device multirate SPH/LS step. */
    void advanceGPU(Real timeStep) override;
    /** Advances one SPH-device/LS-host multirate step. */
    void advanceHybrid(Real timeStep) override;
    /** Exposes dense, current-time fluid output without changing the leading interval. */
    void beginOutputSnapshot() override;
    /** Restores leading host buffers after output without changing the integration. */
    void endOutputSnapshot() override;
    /** Retains a synchronized host view while preserving the leading interval for continuation. */
    void beginSolveCompletionSnapshot() override;
    /** Exposes current-time fluid/boundary state and copies inherited LS storage. */
    void copySystemDeviceToHost(cudaStream_t stream) override;
    /** Writes SPH output plus inherited LS and interaction groups. */
    void writeSystemVTU(int frameIndex) override;
    /** Returns combined SPH, boundary, and inherited LS device memory. */
    double systemDeviceMemoryGB() const noexcept override;
    /** Invalidates cached multirate phase and deferred fluid state. */
    void invalidateDeferredState() noexcept override;
    /** Enables restoration of a valid deferred SPH continuation. */
    void resumeContinuation() noexcept override;

private:
    /** Validated uniform fluid properties and their derived constants. */
    struct SPHProperties {
        /** Returns the uniform mass represented by one initially spaced particle. */
        Real particleMass() const noexcept { return referenceDensity_ * spacing_ * spacing_ * spacing_; }
        /** Creates stage inputs from current properties and gravity, without a second cached copy. */
        SPHInteractionParameters interactionParameters(const Vec3& gravity, const math::periodicDomain& periodicDomain) const noexcept
        {
            return {smoothingLength_, particleMass(), referenceDensity_, latticeKernelSum3D_, soundSpeed_, dynamicViscosity_, gravity, periodicDomain};
        }

        Real maximumVelocity_{1.0};    ///< Configured design velocity scale.
        Real soundSpeed_{10.0};        ///< Artificial speed of sound.
        Real spacing_{0.0};            ///< Uniform initial particle spacing.
        Real smoothingLength_{0.0};    ///< Kernel smoothing length.
        Real referenceDensity_{0.0};   ///< Equation-of-state reference density.
        Real dynamicViscosity_{0.0};   ///< Dynamic viscosity.
        Real latticeKernelSum3D_{0.0}; ///< Reference three-dimensional lattice kernel sum.
        Real viscousTimeScale_{0.0};   ///< Cached viscous stability scale.
        bool set_{false};              ///< Whether uniform properties were validated.
    };

    /** One stability limit and its integer-quantized integration interval. */
    struct SPHTimeInterval {
        Real limit_{0.0};    ///< Unquantized stability limit.
        Real timeStep_{0.0}; ///< Quantized physical duration.
        int interval_{1};    ///< DEM steps for acoustics, acoustic steps for advection.
    };

    /** Mutable multirate integration, stability, and neighborhood state. */
    struct SPHStepState {
        SPHTimeInterval acoustic_;       ///< Acoustic stability and DEM-step quantization.
        SPHTimeInterval advection_;      ///< Advection stability and acoustic-step quantization.
        SPHStateStatistics statistics_;  ///< Latest backend-independent fluid extrema.
        Real neighborSearchRadius_{0.0}; ///< Radius valid for the current neighborhood.
        Real timeInAdvectionStep_{0.0};  ///< Represented time since the last grid rebuild.
        int couplingDEMSteps_{0};        ///< Length of the current leading interval, independent of the next acoustic plan.
        int completedDEMSteps_{0};       ///< Completed solid substeps and their complete left-endpoint kinematic samples.

        bool couplingActive() const noexcept { return completedDEMSteps_ < couplingDEMSteps_; }

        /** Advances the advection phase and reports whether its interval has ended. */
        bool advanceAdvectionPhase(Real timeStep) noexcept
        {
            timeInAdvectionStep_ += timeStep;
            const bool completed = timeInAdvectionStep_ + math::defaultTolerance >= advection_.timeStep_;
            if (completed)
                timeInAdvectionStep_ = 0.0;
            return completed;
        }
    };

    /** Finite-jet configuration and its represented-time completion cache, independent of integration state. */
    struct SPHJetState {
        /** Finite supply, short inlet and stable ownership for one emitter. */
        struct constraint {
            SPHJet value_;                     ///< Outlet, direction, radius, speed and duration.
            Real startTime_{0.0};              ///< Absolute source registration time.
            int axialCount_{0};                ///< Total finite supply, including initial buffer layers.
            int bufferLayerCount_{0};          ///< Bounded upstream inlet depth in particle layers.
            int firstInletLayer_{0};           ///< Released chronological prefix; never re-enters the inlet.
            std::vector<Vec3> offsets_;        ///< World-space transverse lattice offsets for one layer.
            std::vector<int> particleIndices_; ///< Append-only ownership, independent of other sources.

            int layerSize() const noexcept { return static_cast<int>(offsets_.size()); }
            int layerCount() const noexcept { return static_cast<int>(particleIndices_.size()) / layerSize(); }
            int layerBegin(int layer) const noexcept { return particleIndices_[layer * layerSize()]; }
            Real endTime() const noexcept { return startTime_ + value_.duration(); }
            Real nextEmissionTime(Real spacing) const noexcept { return startTime_ + (layerCount() - bufferLayerCount_ + 1) * spacing / value_.speed(); }
            Real releaseTime(int axialIndex, Real spacing) const noexcept { return startTime_ + (axialIndex + 0.5) * spacing / value_.speed(); }
        };

        bool completedAt(Real time) const noexcept { return constraints_.empty() || time >= latestEndTime_; }

        std::vector<constraint> constraints_; ///< Emitters; observation never changes supply or ownership.
        Real latestEndTime_{0.0};             ///< Latest end time among all configured jets.
        bool completed_{true};                ///< Skip inlet checks after the represented state passes the final end.
    };

    /** Host buffers temporarily replaced by a current-time observation; device state is untouched. */
    struct outputSnapshot {
        SPHParticleContainer::host_container_type SPHParticles_;         ///< Saved fluid host state.
        virtualParticleContainer::host_container_type virtualParticles_; ///< Saved boundary host state.
        bool active_{false};                                             ///< Whether a snapshot is currently held.
    };
    bool outputObservationSnapshot_{false}; ///< Whether the observation scope owns the active snapshot.

    /** Initializes persistent SPH state and selects the mode-specific interaction engine. */
    void initializeSPHState();
    /** Creates host neighborhoods and initializes host fluid state. */
    void initializeSPH(cpu::mode);
    /** Shared host-solid KDK schedule for CPU and Hybrid fluid backends. */
    void advanceHostParticles(Real timeStep);
    /** Completes one host SPH and LS force assembly. */
    void assembleForceAndTorque(Real historyTimeStep, cpu::mode);
    /** Completes one device SPH and LS force assembly. */
    void assembleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream);
    /** Saves leading endpoint host buffers and exposes a dense synchronized view. */
    void captureSPHState(cudaStream_t stream);
    /** Restores leading host buffers; no force, history, or device state is changed. */
    void restoreSPHState() noexcept;
    /** Computes a complete leading acoustic interval and replaces only the held hydraulic load. */
    void beginSPHInterval(cudaStream_t stream);
    /** Appends due inlet layers at a common acoustic start without resetting any evolved particle. */
    void emitSPHJetParticles(cudaStream_t stream);
    /** Appends one source layer; registration and runtime injection use the same initialization. */
    void appendSPHJetLayer(SPHJetState::constraint& source);
    /** Clips the next interval to the first DEM node at or after an inlet-layer or source-end event. */
    int SPHJetInterval(int plannedSteps) const noexcept;
    /** Applies virtual-only inlet holes in the current world frame; reports a changed wall quadrature. */
    bool updateSPHJetOpenings(cpu::mode);
    bool updateSPHJetOpenings(gpu::mode, cudaStream_t stream);
    /** Allocates and uploads device fluid and boundary dependencies. */
    void initializeSPHDevice(cudaStream_t stream);
    /** Starts a host advection interval by rebuilding neighborhoods and density. */
    void prepareSPHAdvectionStep(cpu::mode);
    /** Starts a device advection interval by rebuilding grids and density. */
    void prepareSPHAdvectionStep(gpu::mode, cudaStream_t stream);
    /** Invalidates only search caches, without changing fluid phase or forces. */
    void invalidateSPHNeighborhood() noexcept { SPHNeighborhood_.invalidate(); }
    /** Refreshes searches, density and viscous forces after births or changed inlet holes. */
    void invalidateSPHAdvection() noexcept
    {
        invalidateSPHNeighborhood();
        SPHStepState_.timeInAdvectionStep_ = 0.0;
    }
    /** Rebuilds host neighbors only when motion has exhausted the search skin. */
    void ensureSPHNeighborhood(cpu::mode);
    /** Rebuilds device grids only when motion has exhausted the search skin. */
    void ensureSPHNeighborhood(gpu::mode, cudaStream_t stream);
    /** Advances one host acoustic substep without rebuilding the grid. */
    void advanceSPH(Real timeStep, cpu::mode);
    /** Advances one device acoustic substep without rebuilding the grid. */
    void advanceSPH(Real timeStep, gpu::mode, cudaStream_t stream);
    /** Completes jet timing and advances the common acoustic/advection planning state. */
    void finishSPHInterval(Real interval, cudaStream_t stream);
    /** Evaluates the leading host fluid interval against start-time boundary kinematics. */
    void advanceCouplingInterval(Real interval, cpu::mode);
    /** Evaluates the leading device fluid interval against start-time boundary kinematics. */
    void advanceCouplingInterval(Real interval, gpu::mode, cudaStream_t stream);
    /** Evaluates the leading fluid interval with host boundary owners. */
    void advanceCouplingInterval(Real interval, hybrid::mode, cudaStream_t stream);
    /** Updates the acoustic limit from host state statistics. */
    void updateSPHAcousticTimeStep(cpu::mode);
    /** Updates the acoustic limit from device-reduced state statistics. */
    void updateSPHAcousticTimeStep(gpu::mode, cudaStream_t stream);
    /** Checks reduced state and computes the common acoustic limit and quantized duration. */
    void updateSPHAcousticTimeStep(const SPHStateStatistics& statistics);
    /** Computes the shared advection search radius from fluid and boundary extrema. */
    void updateSPHNeighborSearchRadius(const SPHKinematicsStatistics& boundaryStatistics);
    /** Updates the advection limit from cached extrema and fluid properties. */
    void updateSPHAdvectionTimeStep();
    /** Quantizes an acoustic limit to an integer multiple of the DEM step. */
    void setSPHAcousticTimeStepFromLimit(Real maximumTimeStep);
    /** Enforces active virtual-pipe inlet velocities on host particles. */
    void applySPHJetVelocity(Real representedTime, cpu::mode);
    /** Enforces active virtual-pipe inlet velocities on device particles. */
    void applySPHJetVelocity(Real representedTime, gpu::mode, cudaStream_t stream);
    /** Writes one SPH-particle VTU frame. */
    void writeSPHVTU(int frameIndex);
    /** Returns fluid, boundary, and SPH-grid device capacity in GiB. */
    double SPHDeviceMemoryGB() const noexcept;

    using SPHIntervalStorage = hostAoSDeviceSoA<SPHParticle, deviceLayout<pointMass::positionField, pointMass::velocityField, SPHParticle::densityField>>;

    SPHParticleContainer SPHParticles_;                     ///< Fluid particles.
    SPHIntervalStorage SPHIntervalStart_;                   ///< Only position, velocity and density at acoustic start; device-local in CUDA modes.
    virtualParticleContainer virtualParticles_;             ///< Generated LS boundary samples.
    spatialGridContainer SPHSpatialGrid_;                   ///< Fluid background grid.
    SPHNeighborhood SPHNeighborhood_;                       ///< Actual-displacement validity of reused searches.
    spatialGridContainer virtualParticleSpatialGrid_;       ///< Boundary-sample background grid.
    std::vector<SPHParticleVTUField> SPHParticleVTUFields_; ///< Optional SPH output fields.
    virtualParticleCoupling virtualParticleCoupling_;       ///< Host aggregation by LS owner.
    std::unique_ptr<cpu::SPHInteraction> SPHInteractions_;  ///< CPU neighborhood and interaction engine.
    SPHJetState SPHJets_;                                   ///< Finite-jet definitions and completion cache.
    execution::SPHJetOpeningContainer SPHJetOpenings_;      ///< Immutable short-port geometry; wall slots remain stable.
    SPHProperties SPHProperties_;                           ///< Uniform fluid configuration.
    math::periodicDomain SPHPeriodicDomain_{};              ///< Neighbor images; integration keeps continuous physical coordinates.
    SPHStepState SPHStepState_;                             ///< Mutable multirate integration state.
    bool SPHBoundaryStatic_{false};                         ///< Whether all LS owners are stationary infinite-mass boundaries.
    bool SPHBoundaryMotionExpected_{false};                 ///< Whether a currently stationary boundary is scheduled to move.
    outputSnapshot outputSnapshot_;                         ///< Non-destructive synchronization snapshot.
};

} // namespace fundem
