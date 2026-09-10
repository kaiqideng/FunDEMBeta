/**
 * @file structureIceLoadSolver.h
 * @brief Declares level-ice generation and the structure–ice interaction solver.
 */
#pragma once

#include "solver/SphereDEM.h"

#include <memory>
#include <vector>

namespace fundem
{

/** One bond edge in the generated triangular level-ice lattice. */
struct levelIceConnection {
    int firstParticleIndex_{-1};  ///< First generated sphere index.
    int secondParticleIndex_{-1}; ///< Second generated sphere index.
};

/**
 * Generates a one-sphere-thick triangular lattice representing level ice and
 * evaluates the bond parameters used by the original IceDEM model.
 */
class levelIce
{
public:
    using Real = math::Real;
    using Vec3 = math::Vec3;

    explicit levelIce(Real thickness = 0.1, Real density = 910.0);

    Real thickness() const noexcept { return thickness_; }
    Real density() const noexcept { return density_; }
    const Vec3& minimumBoundary() const noexcept { return minimumBoundary_; }
    const Vec3& maximumBoundary() const noexcept { return maximumBoundary_; }
    const std::vector<Vec3>& particlePositions() const noexcept { return particlePositions_; }
    const std::vector<levelIceConnection>& connections() const noexcept { return connections_; }

    void setThickness(Real value);
    void setDensity(Real value);

    /**
     * Builds the planar hexagonal lattice and places its center elevation from
     * exact spherical buoyancy equilibrium.
     */
    void build(Real minimumX, Real maximumX, Real minimumY, Real maximumY, Real waterDensity, Real waterLevel);

    Real elementMass() const noexcept;
    Real elementRadius() const noexcept { return 0.5 * thickness_; }
    Real modeICriticalEnergyReleaseRate(Real flexuralStrength, Real youngsModulus) const noexcept;
    Real normalStiffness(Real youngsModulus) const noexcept;
    Real shearStiffness(Real youngsModulus) const noexcept;
    Real bendingStiffness(Real youngsModulus) const noexcept;
    Real torsionalStiffness(Real youngsModulus, Real poissonRatio) const noexcept;

private:
    /** Clears generated lattice data after a property change. */
    void clear() noexcept;
    Real equilibriumCenterZ(Real waterDensity, Real waterLevel) const noexcept;

    Real thickness_{0.1};                         ///< Ice-sheet thickness and element diameter.
    Real density_{910.0};                         ///< Ice density.
    Vec3 minimumBoundary_{Vec3::zero()};          ///< Bounds of generated spherical elements.
    Vec3 maximumBoundary_{Vec3::zero()};          ///< Bounds of generated spherical elements.
    std::vector<Vec3> particlePositions_;         ///< Generated sphere centers.
    std::vector<levelIceConnection> connections_; ///< Unique nearest-neighbor bond edges.
};

/**
 * SphereDEM specialization for level-ice interaction with LS structures.
 *
 * Every sphere receives buoyancy and quadratic water drag. Every LSParticle is
 * treated as one structure component, and its resultant force is averaged over
 * the DEM steps between output frames into structureIceLoad.dat.
 */
class structureIceLoadSolver : public SphereDEM
{
public:
    explicit structureIceLoadSolver(int device = 0);
    structureIceLoadSolver(executionMode mode, int device = 0);
    ~structureIceLoadSolver() override;

    structureIceLoadSolver(const structureIceLoadSolver&) = delete;
    structureIceLoadSolver& operator=(const structureIceLoadSolver&) = delete;
    structureIceLoadSolver(structureIceLoadSolver&&) = delete;
    structureIceLoadSolver& operator=(structureIceLoadSolver&&) = delete;

    levelIce& ice() noexcept { return ice_; }
    const levelIce& ice() const noexcept { return ice_; }
    const Vec3& waterVelocity() const noexcept { return waterVelocity_; }
    Real waterDensity() const noexcept { return waterDensity_; }
    Real waterLevel() const noexcept { return waterLevel_; }
    Real dragCoefficient() const noexcept { return dragCoefficient_; }

    /** Replaces the uniform water velocity, density, drag coefficient, and free-surface elevation. */
    void setWaterCondition(const Vec3& velocity, Real density = 1000.0, Real dragCoefficient = 0.1, Real waterLevel = 0.0);

protected:
    void initialize(cpu::mode mode) override;
    void initialize(gpu::mode mode, cudaStream_t stream) override;
    void initialize(hybrid::mode mode, cudaStream_t stream) override;
    void addSphereExternalForceAndTorque(particleContainer::host_container_type& particles) override;
    void addSphereExternalForceAndTorque(particle::device_type particles, cudaStream_t stream) override;
    void addLSParticleExternalForceAndTorque(LSParticleContainer::host_container_type& structures) override;
    void addLSParticleExternalForceAndTorque(LSParticle::device_type structures, cudaStream_t stream) override;
    void writeSystemVTU(int frameIndex) override;
    double systemDeviceMemoryGB() const noexcept override;

private:
    struct loadAccumulator;

    void ensureLoadAccumulator(int structureCount, bool useDevice, cudaStream_t stream);
    void writeStructureIceLoad(int frameIndex);
    void finishInitialization(bool useDevice, cudaStream_t stream);

    levelIce ice_;                                     ///< Level-ice lattice generator and parameter helper.
    Vec3 waterVelocity_{Vec3::zero()};                 ///< Uniform water velocity.
    Real waterDensity_{1000.0};                        ///< Water density.
    Real waterLevel_{0.0};                             ///< Horizontal free-surface elevation.
    Real dragCoefficient_{0.1};                        ///< Translational and rotational drag coefficient.
    std::unique_ptr<loadAccumulator> loadAccumulator_; ///< Host/device interval force sums.
    bool initializationAssembly_{false};               ///< Suppresses load sampling during initial force assembly.
};

} // namespace fundem
