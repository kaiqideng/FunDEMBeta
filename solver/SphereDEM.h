/**
 * @file SphereDEM.h
 * @brief Declares sphere and sphere-level-set DEM solver extensions.
 */
#pragma once

#include "particle/particle.h"
#include "solver/LSDEM.h"

#include <vector>

namespace fundem
{

/**
 * DEM solver extending LSDEM with spherical particles, sphere-sphere contacts,
 * sphere-LS contacts, and all three supported execution policies.
 */
class SphereDEM : public LSDEM
{
public:
    /** Exposes the common lazy initialization entry point. */
    void initialize() { solver::initialize(); }

    /** Creates a CPU sphere/LS solver targeting @p device for later GPU use. */
    explicit SphereDEM(int device = 0) : LSDEM(device) {}
    /** Creates a sphere/LS solver in the requested execution mode. */
    SphereDEM(executionMode mode, int device = 0) : LSDEM(device) { setExecutionMode(mode); }
    ~SphereDEM() override = default;

    const particleContainer& spheres() const noexcept { return spheres_; }
    const interactionContainer& sphereInteractions() const noexcept { return sphereInteractions_; }
    const interactionContainer& sphereLSInteractions() const noexcept { return sphereLSInteractions_; }

    using LSDEM::addMaterial;
    /** Appends an ordinary sphere material and returns its stable index. */
    int addMaterial(const material& value) { return appendMaterial(value); }
    /** Appends a valid spherical particle and returns its stable index. */
    int addSphere(const particle& value);
    /** Replaces one sphere position before or between solve calls. */
    void setSpherePosition(int particleIndex, const Vec3& value);
    /** Replaces one sphere linear velocity before or between solve calls. */
    void setSphereVelocity(int particleIndex, const Vec3& value);
    /** Replaces one sphere angular velocity before or between solve calls. */
    void setSphereAngularVelocity(int particleIndex, const Vec3& value);
    /** Dispatches a valid bond to the sphere-sphere, sphere-LS, or LS-LS store. */
    int addBond(const bond& value) override;

    void setSphereVTUFields(std::vector<sphereVTUField> fields);

protected:
    /** Mutable sphere storage for derived frontends implementing prescribed kinematics. */
    particleContainer& mutableSpheres() noexcept { return spheres_; }
    /** Mutable sphere interaction storage for derived frontends restoring a saved host state. */
    interactionContainer& mutableSphereInteractions() noexcept { return sphereInteractions_; }
    /** Mutable sphere-LS interaction storage for derived frontends restoring a saved host state. */
    interactionContainer& mutableSphereLSInteractions() noexcept { return sphereLSInteractions_; }
    /** Aggregates sphere and inherited LS solid-system energy. */
    energyRecord systemEnergy() const override;
    /** Accepts CPU, GPU, and sphere-GPU/LS-CPU Hybrid modes. */
    void validateExecutionMode(executionMode value) const override;
    /** Updates sphere and LS storage selection after a mode change. */
    void executionModeChanged(executionMode value) override;
    /** Initializes both sphere and inherited LS host state. */
    void initialize(cpu::mode) override;
    /** Initializes both sphere and inherited LS device state. */
    void initialize(gpu::mode, cudaStream_t stream) override;
    /** Initializes sphere device state and LS host state. */
    void initialize(hybrid::mode, cudaStream_t stream) override;
    /** Advances one all-host sphere/LS step. */
    void advanceCPU(Real timeStep) override;
    /** Advances one all-device sphere/LS step. */
    void advanceGPU(Real timeStep) override;
    /** Advances one sphere-device/LS-host step. */
    void advanceHybrid(Real timeStep) override;
    /** Enqueues sphere state plus inherited LS state to the host. */
    void copySystemDeviceToHost(cudaStream_t stream) override;
    /** Writes sphere output plus inherited LS and interaction output. */
    void writeSystemVTU(int frameIndex) override;
    /** Returns combined sphere and inherited LS device memory. */
    double systemDeviceMemoryGB() const noexcept override;

    /** Adds user-defined sphere loads after all built-in force assembly. */
    virtual void addSphereExternalForceAndTorque(particleContainer::host_container_type&) {}
    /** Device counterpart of the external sphere load hook; work must use @p stream. */
    virtual void addSphereExternalForceAndTorque(particle::device_type, cudaStream_t) {}

private:
    /** Returns one mutable sphere after checked stable-index lookup. */
    particle& mutableSphere(int particleIndex);
    /** Assembles all sphere/LS forces in CPU mode. */
    void assembleForceAndTorque(Real historyTimeStep, cpu::mode);
    /** Assembles all sphere/LS forces in GPU mode. */
    void assembleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream);
    /** Assembles sphere device and LS host forces in Hybrid mode. */
    void assembleForceAndTorque(Real historyTimeStep, hybrid::mode, cudaStream_t stream);
    /** Searches and evaluates sphere interactions on the host. */
    void calculateSphereForceAndTorque(Real historyTimeStep, cpu::mode);
    /** Searches and evaluates sphere interactions on the device. */
    void calculateSphereForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream);
    /** Evaluates sphere contacts with hybrid synchronization. */
    void calculateSphereForceAndTorque(Real historyTimeStep, hybrid::mode, cudaStream_t stream);
    /** Uploads sphere, interaction, and spatial-grid state. */
    void initializeSphereDevice(cudaStream_t stream);
    /** Enqueues current sphere and sphere-interaction state to the host. */
    void copySphereDeviceToHost(cudaStream_t stream);
    /** Writes sphere particle/contact/bond groups for one frame. */
    void writeSphereVTU(int frameIndex);
    /** Returns sphere-specific device capacity in GiB. */
    double sphereDeviceMemoryGB() const noexcept;

    particleContainer spheres_;                   ///< Spherical rigid particles.
    interactionContainer sphereInteractions_;     ///< Sphere-sphere contacts, bonds, and histories.
    interactionContainer sphereLSInteractions_;   ///< Sphere-LS contacts and bonds.
    spatialGridContainer sphereSpatialGrid_;      ///< Device multi-level sphere grid.
    contactSearch sphereSearch_;                  ///< CPU sphere-sphere search.
    contactSearch sphereLSSearch_;                ///< CPU sphere-LS search.
    std::vector<sphereVTUField> sphereVTUFields_; ///< Optional sphere output fields.
};

} // namespace fundem
