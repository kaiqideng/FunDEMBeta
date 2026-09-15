/**
 * @file solver.h
 * @brief Declares the common solver lifecycle, scheduling, backend, and output API.
 */
#pragma once

#include "data/CudaTypes.h"
#include "data/vtuWriter.h"
#include "math/Vector3.h"

#include <functional>
#include <memory>
#include <string>

namespace fundem
{

struct energyRecord;

namespace cpu
{
/** Compile-time tag selecting host-only solver implementations. */
struct mode {};
} // namespace cpu

namespace gpu
{
/** Compile-time tag selecting device-only solver implementations. */
struct mode {};
} // namespace gpu

namespace hybrid
{
/** Compile-time tag selecting GPU primary particles with CPU LS solids. */
struct mode {};
} // namespace hybrid

/** Runtime execution policy selected once during solver initialization. */
enum class executionMode
{
    CPU,
    GPU,
    Hybrid
};

/** Common lifecycle, time integration, output scheduling, and GPU context for all solvers. */
class solver
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    /**
     * Creates a solver targeting one CUDA device and execution policy.
     * CPU mode does not create a CUDA context.
     * @param device Zero-based CUDA device index.
     * @param mode Requested execution policy.
     */
    explicit solver(int device = 0, executionMode mode = executionMode::CPU);
    virtual ~solver();
    solver(const solver&) = delete;
    solver& operator=(const solver&) = delete;
    solver(solver&&) = delete;
    solver& operator=(solver&&) = delete;

    const Vec3& minimumBoundary() const noexcept { return minimumBoundary_; }
    const Vec3& maximumBoundary() const noexcept { return maximumBoundary_; }
    const Vec3& gravity() const noexcept { return gravity_; }
    Real timeStep() const noexcept { return timeStep_; }
    Real time() const noexcept { return time_; }
    int stepCount() const noexcept { return stepCount_; }
    int device() const noexcept { return device_; }
    executionMode mode() const noexcept { return executionMode_; }
    bool initialized() const noexcept { return initialized_; }
    const std::string& outputDirectory() const noexcept { return outputDirectory_; }
    int outputStepInterval() const noexcept { return outputStepInterval_; }
    int outputFrameCount() const noexcept { return outputFrameCount_; }
    vtuFormat vtuOutputFormat() const noexcept { return vtuOutputFormat_; }

    void setBoundary(const Vec3& minimumBoundary, const Vec3& maximumBoundary);
    void setGravity(const Vec3& value);
    void setTimeStep(Real value);
    void setExecutionMode(executionMode value);
    void setOutputDirectory(std::string directory);
    void setOutputStepInterval(int stepInterval);
    void setVTUOutputFormat(vtuFormat value) noexcept { vtuOutputFormat_ = value; }

    /** Writes one output frame immediately from a synchronized snapshot. */
    void writeOutput();
    /**
     * Calls @p observer with host state synchronized to time(), temporarily
     * completing deferred multirate work without changing its continuation.
     * Nested observations and writeOutput() share the same snapshot. The
     * observer must not advance or modify the solver; copied values remain
     * valid after return, but references to temporary host state do not.
     * Deferred state is restored even when the observer throws.
     * Reentrant initialize(), step(), and solve() calls are rejected.
     */
    void observeCurrentState(const std::function<void(const solver&)>& observer);
    /**
     * Advances @p numberOfSteps DEM steps, initializing or resuming as needed.
     * Repeated calls continue from the current time and permit append-only model
     * additions between solves.
     */
    void solve(int numberOfSteps);
    /** Validates and initializes the selected execution policy without advancing. */
    void initialize();
    /** Advances exactly one DEM time step. */
    void step();
    /** Waits for all queued solver work in the selected backend. */
    void synchronize();
    /** Synchronizes solver-owned device state back to host containers. */
    void copyDeviceToHost();
    /** Returns total currently allocated solver device memory in GiB. */
    double deviceMemoryGB() const noexcept;

protected:
    cudaStream_t gpuStream() const noexcept;
    /** Selects the configured CUDA device for the calling host thread. */
    void activateGPUDevice() const;
    /** Throws when the current build cannot execute the requested GPU mode. */
    void validateGPUEnable() const;
    /** Constructs an output path under the requested group and frame number. */
    std::string makeVTUFileName(const char* group, const char* name, int frameIndex) const;
    /** Marks deferred integration state invalid after user model changes. */
    void invalidateContinuation() noexcept;

private:
    using advanceTemplate = void (solver::*)(Real);

    bool usesDevice() const noexcept { return executionMode_ != executionMode::CPU; }
    /** Derived-class mode validation hook. */
    virtual void validateExecutionMode(executionMode value) const = 0;
    /** Allows derived containers to follow a newly selected execution mode. */
    virtual void executionModeChanged(executionMode) {}
    /** Validates derived model invariants immediately before initialization. */
    virtual void validateSystemConfiguration() const {}
    /** Initializes a derived system for host execution. */
    virtual void initialize(cpu::mode) = 0;
    /** Initializes a derived system for device execution on @p stream. */
    virtual void initialize(gpu::mode, cudaStream_t stream) = 0;
    /** Initializes a derived mixed host/device system on @p stream. */
    virtual void initialize(hybrid::mode, cudaStream_t stream) = 0;
    /** Advances a derived host system by one DEM step. */
    virtual void advanceCPU(Real timeStep) = 0;
    /** Advances a derived device system by one DEM step. */
    virtual void advanceGPU(Real timeStep) = 0;
    /** Advances a derived hybrid system by one DEM step. */
    virtual void advanceHybrid(Real timeStep) = 0;
    /** Captures deferred state before a non-destructive output snapshot. */
    virtual void beginOutputSnapshot() {}
    /** Restores deferred state after an output snapshot. */
    virtual void endOutputSnapshot() {}
    /** Captures deferred state before final solve synchronization. */
    virtual void beginSolveCompletionSnapshot() {}
    /** Restores or finalizes deferred state after solve completion. */
    virtual void endSolveCompletionSnapshot() {}
    /** Enqueues derived device-to-host state transfer. */
    virtual void copySystemDeviceToHost(cudaStream_t) {}
    /** Writes all derived particle and interaction VTU groups. */
    virtual void writeSystemVTU(int frameIndex) = 0;
    /** Returns the derived solid-system energy record. */
    virtual energyRecord systemEnergy() const = 0;
    /** Returns allocated device memory owned by the derived system. */
    virtual double systemDeviceMemoryGB() const noexcept { return 0.0; }
    /** Invalidates derived deferred/multirate continuation state. */
    virtual void invalidateDeferredState() noexcept {}
    /** Marks derived deferred state eligible for continuation. */
    virtual void resumeContinuation() noexcept {}

    /** Validates common time-step, boundary, output, and mode configuration. */
    void validateConfiguration() const;
    /** Rejects integration lifecycle changes while temporary state is visible. */
    void requireNoCurrentStateObservation() const;
    /** Lazily creates the configured CUDA-device stream. */
    void initializeDeviceContext();
    /** Resolves and creates the configured output-directory structure. */
    void ensureOutputDirectory() const;
    /** Removes VTU and DAT files once at the first solve. */
    void clearOutputFiles();
    /** Brings the selected backend to a coherent host-visible state. */
    void synchronizeHostState();
    /** Writes a scheduled frame after the current completed step. */
    void writeOutputIfDue();

    struct gpuState;

    std::unique_ptr<gpuState> gpu_;                   ///< Lazily created CUDA stream and device state.
    Vec3 minimumBoundary_{Vec3::zero()};              ///< Lower search-domain bound.
    Vec3 maximumBoundary_{Vec3{1.0, 1.0, 1.0}};       ///< Upper search-domain bound.
    Vec3 gravity_{Vec3{0.0, 0.0, -9.81}};             ///< Uniform gravitational acceleration.
    Real timeStep_{1.0e-6};                           ///< DEM integration step.
    Real time_{0.0};                                  ///< Current represented simulation time.
    int stepCount_{0};                                ///< Number of completed DEM steps.
    int device_{0};                                   ///< Zero-based CUDA device index.
    executionMode executionMode_{executionMode::CPU}; ///< Selected runtime policy.
    std::string outputDirectory_;                     ///< Root output directory.
    vtuFormat vtuOutputFormat_{vtuFormat::binary};    ///< VTU serialization format.
    int outputStepInterval_{0};                       ///< DEM steps between output frames.
    int nextOutputStep_{0};                           ///< Next scheduled output step.
    int outputFrameCount_{0};                         ///< Number of frames already written.
    bool initialized_{false};                         ///< Whether backend state matches the host model.
    bool firstSolve_{true};                           ///< Whether output cleanup is still pending.
    bool continuationValid_{false};                   ///< Whether deferred state may resume.
    bool observingCurrentState_{false};               ///< Whether a synchronized observation scope is active.
    advanceTemplate advanceTemplate_{nullptr};        ///< Cached mode-specific advance dispatch.
};

} // namespace fundem
