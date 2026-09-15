/**
 * @file LSDEM.h
 * @brief Declares the level-set DEM solver and shared solid-system services.
 */
#pragma once

#include "data/vtuOutput.h"
#include "geometry/LSGeometry.h"
#include "interaction/contactSearch.h"
#include "interaction/interactionContainer.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "particle/spatialGrid.h"
#include "solver/solver.h"

#include <vector>

namespace fundem
{

namespace levelset
{
class LSInfo;
}

/**
 * Base level-set DEM solver owning materials, reusable geometries, rigid
 * LS particles, interactions, and CPU/GPU search infrastructure.
 */
class LSDEM : public solver
{
public:
    /** Exposes the common lazy initialization entry point. */
    void initialize() { solver::initialize(); }

    using Real = math::Real;
    using Vec3 = math::Vec3;

    /** Creates a CPU LS solver targeting @p device if GPU mode is selected later. */
    explicit LSDEM(int device = 0) : solver(device) {}
    ~LSDEM() override = default;
    LSDEM(const LSDEM&) = delete;
    LSDEM& operator=(const LSDEM&) = delete;
    LSDEM(LSDEM&&) = delete;
    LSDEM& operator=(LSDEM&&) = delete;

    const materialContainer& materials() const noexcept { return materials_; }
    const LSGeometryContainer& geometries() const noexcept { return geometries_; }
    const LSParticleContainer& LSParticles() const noexcept { return LSParticles_; }
    const interactionContainer& LSParticleInteractions() const noexcept { return LSParticleInteractions_; }

    /** Appends a valid level-set material and returns its stable index. */
    int addMaterial(const LSMaterial& value);
    /** Converts and appends one reusable level-set input geometry. */
    int addGeometry(const levelset::LSInfo& value, bool isFixed = false);
    /** Appends a valid LS particle whose material and geometry belong here. */
    int addLSParticle(const LSParticle& value);
    /** Appends a valid LS-LS bond referencing this solver's particle container. */
    virtual int addBond(const bond& value);

    /** Replaces one LS-particle position before or between solve calls. */
    void setLSParticlePosition(int particleIndex, const Vec3& value);
    /** Replaces one LS-particle linear velocity before or between solve calls. */
    void setLSParticleVelocity(int particleIndex, const Vec3& value);
    /** Replaces one LS-particle angular velocity before or between solve calls. */
    void setLSParticleAngularVelocity(int particleIndex, const Vec3& value);

    void setLSParticleVTUFields(std::vector<LSParticleVTUField> fields);
    void setContactVTUFields(std::vector<contactVTUField> fields);
    void setBondVTUFields(std::vector<bondVTUField> fields);

protected:
    /** Constructs an LS base directly in a derived solver's execution mode. */
    LSDEM(int device, executionMode mode) : solver(device, mode) { LSParticles_.setUseDevice(mode == executionMode::GPU); }
    /** Aggregates LS-particle, contact, and bond energy. */
    energyRecord systemEnergy() const override;
    /** Accepts only execution modes supported by a pure LS solver. */
    void validateExecutionMode(executionMode value) const override;
    /** Updates LS container storage selection after a mode change. */
    void executionModeChanged(executionMode value) override;
    /** Initializes LS search and interaction state on the host. */
    void initialize(cpu::mode) override;
    /** Initializes LS geometry, particles, search, and interactions on the device. */
    void initialize(gpu::mode, cudaStream_t stream) override;
    /** Rejects unsupported pure-LS hybrid initialization. */
    void initialize(hybrid::mode, cudaStream_t stream) override;
    /** Advances one host LS-DEM step. */
    void advanceCPU(Real timeStep) override;
    /** Advances one device LS-DEM step. */
    void advanceGPU(Real timeStep) override;
    /** Rejects unsupported pure-LS hybrid advancement. */
    void advanceHybrid(Real timeStep) override;
    /** Enqueues LS particle and interaction copies to host. */
    void copySystemDeviceToHost(cudaStream_t stream) override;
    /** Writes LS particle, contact, and bond VTU groups. */
    void writeSystemVTU(int frameIndex) override;
    /** Returns memory allocated by LS device containers. */
    double systemDeviceMemoryGB() const noexcept override;

    /** Adds user-defined LSParticle loads after all built-in force assembly. */
    virtual void addLSParticleExternalForceAndTorque(LSParticleContainer::host_container_type&) {}
    /** Device counterpart of the external LS load hook; work must use @p stream. */
    virtual void addLSParticleExternalForceAndTorque(LSParticle::device_type, cudaStream_t) {}
    /** Shared validated append path used by derived solvers. */
    int appendMaterial(const material& value);
    /** Verifies that a bond type and container ownership match this solver. */
    void validateBondForAddition(const bond& value) const;
    /** Appends a validated LS-LS bond and invalidates backend state. */
    int appendLSParticleBond(const bond& value);

    LSParticleContainer& mutableLSParticles() noexcept { return LSParticles_; }
    interactionContainer& mutableLSParticleInteractions() noexcept { return LSParticleInteractions_; }
    spatialGridContainer& mutableLSParticleSpatialGrid() noexcept { return LSParticleSpatialGrid_; }
    const std::vector<contactVTUField>& contactVTUFields() const noexcept { return contactVTUFields_; }
    const std::vector<bondVTUField>& bondVTUFields() const noexcept { return bondVTUFields_; }

    /** Assembles LS contact, bond, and user loads on the host; gravity is applied during velocity integration. */
    void assembleLSParticleForceAndTorque(Real historyTimeStep, cpu::mode);
    /** Assembles LS contact, bond, and user loads on the device; gravity is applied during velocity integration. */
    void assembleLSParticleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream);
    /** Searches and evaluates LS contacts and bonds on the host. */
    void calculateLSParticleForceAndTorque(Real historyTimeStep, cpu::mode);
    /** Searches and evaluates LS contacts and bonds on the device. */
    void calculateLSParticleForceAndTorque(Real historyTimeStep, gpu::mode, cudaStream_t stream);

    /** Uploads LS geometry, particle, mapping, interaction, and grid state. */
    void initializeLSParticleDevice(cudaStream_t stream);
    /** Copies only LS force and torque fields required by hybrid assembly. */
    void copyLSParticleForceAndTorqueToHost(cudaStream_t stream);

private:
    /** Returns one mutable LS particle after checked stable-index lookup. */
    LSParticle& mutableLSParticle(int particleIndex);

    materialContainer materials_;                         ///< Append-only shared material definitions.
    LSGeometryContainer geometries_;                      ///< Append-only reusable level-set geometries.
    LSParticleContainer LSParticles_;                     ///< Rigid level-set particles.
    interactionContainer LSParticleInteractions_;         ///< LS-LS contacts, bonds, ranges, and histories.
    spatialGridContainer LSParticleSpatialGrid_;          ///< Device multi-level LS search grid.
    contactSearch LSParticleSearch_;                      ///< CPU LS-LS contact search.
    std::vector<LSParticleVTUField> LSParticleVTUFields_; ///< Optional LS output fields.
    std::vector<contactVTUField> contactVTUFields_;       ///< Optional contact output fields.
    std::vector<bondVTUField> bondVTUFields_;             ///< Optional bond output fields.
};

} // namespace fundem
