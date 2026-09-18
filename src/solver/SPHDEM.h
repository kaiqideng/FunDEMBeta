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
    /** Adds a finite-duration jet backed by an upstream virtual pipe and returns its first particle index. */
    int addSPHJet(const SPHJet& value);
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
    /** Sets the configured SPH velocity scale and resets the sound speed to ten times this value. */
    void setSPHMaximumVelocity(Real value);
    void setSPHSoundSpeed(Real value);
    void setSPHParticleVTUFields(std::vector<SPHParticleVTUField> fields);

protected:
    /** Mutable fluid storage for derived frontends restoring a saved host state. */
    SPHParticleContainer& mutableSPHParticles() noexcept { return SPHParticles_; }
    /** Declares scheduled boundary motion that may start after initialization. */
    void setSPHBoundaryMotionExpected(bool value) noexcept { SPHBoundaryMotionExpected_ = value; }
    /** Returns inherited LS solid energy; fluid energy is intentionally excluded. */
    energyRecord systemEnergy() const override;
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
    /** Advances one all-device multirate SPH/LS step. */
    void advanceGPU(Real timeStep) override;
    /** Advances one SPH-device/LS-host multirate step. */
    void advanceHybrid(Real timeStep) override;
    /** Captures deferred SPH state before host-visible output. */
    void beginOutputSnapshot() override;
    /** Restores deferred SPH state after output. */
    void endOutputSnapshot() override;
    /** Flushes and captures deferred state before solve completion. */
    void beginSolveCompletionSnapshot() override;
    /** Restores continuation state after solve completion. */
    void endSolveCompletionSnapshot() override;
    /** Enqueues fluid, boundary, and inherited LS state to the host. */
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
        SPHInteractionParameters interactionParameters(const Vec3& gravity) const noexcept
        {
            return {smoothingLength_, particleMass(), referenceDensity_, latticeKernelSum3D_, soundSpeed_, dynamicViscosity_, gravity};
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
        Real representedTime_{0.0};      ///< Physical time represented by the current fluid state.
        int pendingDEMSteps_{0};         ///< Deferred DEM steps not yet flushed to SPH.
        Real couplingImpulseTime_{0.0};  ///< Endpoint matching interval; zero for observation-only flushes.
        bool resume_{false};             ///< Resume rather than reinitialize deferred state.

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

    /** Cached topology and boundary classification from the latest initialization. */
    struct SPHInitializationState {
        int fluidParticleCount_{-1};    ///< Fluid count used by current backend state.
        int boundaryParticleCount_{-1}; ///< LS solid count used by current boundary state.
        bool boundaryStatic_{false};    ///< Whether every LS owner is stationary with infinite mass.
    };

    /** Finite-jet configuration and its represented-time completion cache, independent of integration state. */
    struct SPHJetState {
        /** Immutable particle range and timing bound to one finite-duration jet. */
        struct constraint {
            SPHJet value_;              ///< Outlet, direction, radius, speed and duration.
            Real pipeLength_{0.0};      ///< Axial length of the upstream virtual pipe.
            Real endTime_{0.0};         ///< Absolute time at which the inlet constraint ends.
            int particleIndexBegin_{0}; ///< First SPH particle owned by this jet.
            int particleCount_{0};      ///< Number of particles owned by this jet.
        };

        bool completedAt(Real time) const noexcept { return constraints_.empty() || time >= latestEndTime_; }

        std::vector<constraint> constraints_; ///< Fixed particle ranges; unchanged during output snapshots.
        Real latestEndTime_{0.0};             ///< Latest end time among all configured jets.
        bool completed_{true};                ///< Skip inlet checks after the represented state passes the final end.
    };

    /** Complete deferred SPH and coupling state used for non-destructive output snapshots. */
    struct outputSnapshot {
        SPHParticleContainer::host_container_type SPHParticles_;         ///< Saved fluid host state.
        virtualParticleContainer::host_container_type virtualParticles_; ///< Saved boundary host state.
        SPHStepState stepState_;                                         ///< Saved multirate integration state.
        virtualParticleCoupling::state couplingState_;                   ///< Saved wall-coupling accumulators.
        bool jetsCompleted_{true};                                       ///< Saved jet completion; immutable jet configuration is shared.
        bool active_{false};                                             ///< Whether a snapshot is currently held.
        bool preserveCurrentHost_{false};                                ///< Whether snapshot capture retained host results.
    };
    bool outputObservationSnapshot_{false}; ///< Whether the observation scope owns the active snapshot.

    /** Initializes persistent SPH state and selects the mode-specific interaction engine. */
    void initializeSPHState();
    /** Creates host neighborhoods and initializes host fluid state. */
    void initializeSPH(cpu::mode);
    /** Uploads and initializes complete device fluid state. */
    void initializeSPH(gpu::mode, cudaStream_t stream);
    /** Uploads fluid while retaining LS boundary owners on the host. */
    void initializeSPH(hybrid::mode, cudaStream_t stream);
    /** Completes one host SPH and LS force assembly. */
    void assembleForceAndTorque(Real historyTimeStep, cpu::mode);
    /** Completes one device SPH and LS force assembly. */
    void assembleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream);
    /** Completes one device-fluid/host-solid force assembly. */
    void assembleForceAndTorque(Real historyTimeStep, hybrid::mode, cudaStream_t stream);
    /** Captures deferred multirate state around output or solve completion. */
    void captureSPHState(bool preserveCurrentHost, cudaStream_t stream);
    /** Restores the latest captured multirate state. */
    void restoreSPHState(cudaStream_t stream);
    /** Flushes deferred SPH work to the solver's represented current time. */
    void flushSPHToCurrentTime(cudaStream_t stream, bool commitCouplingImpulse = false);
    /** Allocates and uploads device fluid and boundary dependencies. */
    void initializeSPHDevice(cudaStream_t stream);
    /** Evaluates one accumulated SPH force interval on the host. */
    void calculateSPHForceAndTorque(Real DEMTimeStep, cpu::mode);
    /** Evaluates one accumulated SPH force interval on the device. */
    void calculateSPHForceAndTorque(Real DEMTimeStep, gpu::mode, cudaStream_t stream);
    /** Evaluates one device-fluid/host-solid coupling interval. */
    void calculateSPHForceAndTorque(Real DEMTimeStep, hybrid::mode, cudaStream_t stream);
    /** Accumulates a DEM interval and flushes fluid work when the acoustic interval is ready. */
    void advanceSPHIfReady(Real DEMTimeStep, cudaStream_t stream);
    /** Starts a host advection interval by rebuilding neighborhoods and density. */
    void prepareSPHAdvectionStep(cpu::mode);
    /** Starts a device advection interval by rebuilding grids and density. */
    void prepareSPHAdvectionStep(gpu::mode, cudaStream_t stream);
    /** Invalidates only search caches, without changing fluid phase or forces. */
    void invalidateSPHNeighborhood() noexcept { SPHNeighborhood_.invalidate(); }
    /** Rebuilds host neighbors only when motion has exhausted the search skin. */
    void ensureSPHNeighborhood(cpu::mode);
    /** Rebuilds device grids only when motion has exhausted the search skin. */
    void ensureSPHNeighborhood(gpu::mode, cudaStream_t stream);
    /** Advances one host acoustic substep without rebuilding the grid. */
    void advanceSPH(Real timeStep, cpu::mode);
    /** Advances one device acoustic substep without rebuilding the grid. */
    void advanceSPH(Real timeStep, gpu::mode, cudaStream_t stream);
    /** Flushes a partially represented host SPH interval. */
    void flushSPH(Real representedTime, cpu::mode);
    /** Flushes a partially represented device SPH interval. */
    void flushSPH(Real representedTime, gpu::mode, cudaStream_t stream);
    /** Flushes a partially represented hybrid SPH interval. */
    void flushSPH(Real representedTime, hybrid::mode, cudaStream_t stream);
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
    /** Enqueues current fluid state to the host. */
    void copySPHToHost(cudaStream_t stream);
    /** Writes one SPH-particle VTU frame. */
    void writeSPHVTU(int frameIndex);
    /** Returns fluid, boundary, and SPH-grid device capacity in GiB. */
    double SPHDeviceMemoryGB() const noexcept;

    SPHParticleContainer SPHParticles_;                     ///< Fluid particles.
    virtualParticleContainer virtualParticles_;             ///< Generated LS boundary samples.
    spatialGridContainer SPHSpatialGrid_;                   ///< Fluid background grid.
    SPHNeighborhood SPHNeighborhood_;                       ///< Actual-displacement validity of reused searches.
    spatialGridContainer virtualParticleSpatialGrid_;       ///< Boundary-sample background grid.
    std::vector<SPHParticleVTUField> SPHParticleVTUFields_; ///< Optional SPH output fields.
    virtualParticleCoupling virtualParticleCoupling_;       ///< Host aggregation by LS owner.
    std::unique_ptr<cpu::SPHInteraction> SPHInteractions_;  ///< CPU neighborhood and interaction engine.
    SPHJetState SPHJets_;                                   ///< Finite-jet definitions and completion cache.
    SPHProperties SPHProperties_;                           ///< Uniform fluid configuration.
    SPHStepState SPHStepState_;                             ///< Mutable multirate integration state.
    SPHInitializationState SPHInitializationState_;         ///< Cached initialized topology state.
    bool SPHBoundaryMotionExpected_{false};                 ///< Whether a currently stationary boundary is scheduled to move.
    outputSnapshot outputSnapshot_;                         ///< Non-destructive synchronization snapshot.
};

} // namespace fundem
