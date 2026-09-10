#include "structureIceLoadSolver.h"

#include "data/HostAoSDeviceSoA.h"
#include "structureIceLoadFunctions.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <utility>

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "structureIceLoadKernel.cuh"
#endif

namespace fundem
{
namespace
{

constexpr math::Real forceToMegaNewtons = 1.0e-6;

bool positiveFinite(math::Real value) noexcept { return math::isFinite(value) && value > 0.0; }

} // namespace

struct structureIceLoadSolver::loadAccumulator {
    std::vector<Vec3> hostForces_;
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    host_device_detail::deviceArray<Vec3> deviceForces_;
#endif
    int structureCount_{-1};
    int deviceCount_{0};
    int sampleCount_{0};
    bool useDevice_{false};
};

levelIce::levelIce(Real thickness, Real density)
{
    setThickness(thickness);
    setDensity(density);
}

void levelIce::setThickness(Real value)
{
    if (!positiveFinite(value))
    {
        throw std::invalid_argument("Level-ice thickness must be finite and positive.");
    }
    if (value == thickness_)
    {
        return;
    }
    thickness_ = value;
    clear();
}

void levelIce::setDensity(Real value)
{
    if (!positiveFinite(value))
    {
        throw std::invalid_argument("Level-ice density must be finite and positive.");
    }
    if (value == density_)
    {
        return;
    }
    density_ = value;
    clear();
}

void levelIce::clear() noexcept
{
    minimumBoundary_ = Vec3::zero();
    maximumBoundary_ = Vec3::zero();
    particlePositions_.clear();
    connections_.clear();
}

void levelIce::build(Real minimumX, Real maximumX, Real minimumY, Real maximumY, Real waterDensity, Real waterLevel)
{
    if (!math::isFinite(minimumX) || !math::isFinite(maximumX) || !math::isFinite(minimumY) || !math::isFinite(maximumY) || maximumX <= minimumX || maximumY <= minimumY)
    {
        throw std::invalid_argument("Level-ice planar bounds must be finite and ordered.");
    }
    if (!positiveFinite(waterDensity) || !math::isFinite(waterLevel))
    {
        throw std::invalid_argument("Level-ice water density and level are invalid.");
    }

    const Real radius = elementRadius();
    const Real diameter = thickness_;
    const Real centerXMinimum = minimumX + radius;
    const Real centerXMaximum = maximumX - radius;
    const Real centerYMinimum = minimumY + radius;
    const Real centerYMaximum = maximumY - radius;
    if (centerXMaximum < centerXMinimum || centerYMaximum < centerYMinimum)
    {
        throw std::invalid_argument("The planar domain is too small for the level-ice thickness.");
    }

    clear();

    const Real centerZ = equilibriumCenterZ(waterDensity, waterLevel);
    const Real rowSpacing = std::sqrt(3.0) * radius;
    const Real tolerance = math::defaultTolerance * std::max({Real{1.0}, std::abs(maximumX), std::abs(maximumY)});
    const int rowCount = static_cast<int>(std::floor((centerYMaximum - centerYMinimum + tolerance) / rowSpacing)) + 1;
    std::vector<std::vector<int>> rowParticleIndices;
    rowParticleIndices.reserve(static_cast<std::size_t>(rowCount));

    for (int rowIndex = 0; rowIndex < rowCount; ++rowIndex)
    {
        const Real y = centerYMinimum + rowIndex * rowSpacing;
        const Real xOffset = rowIndex % 2 == 0 ? 0.0 : radius;
        const Real firstX = centerXMinimum + xOffset;
        std::vector<int> currentRow;
        if (firstX <= centerXMaximum + tolerance)
        {
            const int columnCount = static_cast<int>(std::floor((centerXMaximum - firstX + tolerance) / diameter)) + 1;
            currentRow.reserve(static_cast<std::size_t>(columnCount));
            for (int columnIndex = 0; columnIndex < columnCount; ++columnIndex)
            {
                currentRow.push_back(static_cast<int>(particlePositions_.size()));
                particlePositions_.push_back({firstX + columnIndex * diameter, y, centerZ});
            }
        }
        rowParticleIndices.push_back(std::move(currentRow));
    }

    if (particlePositions_.empty())
    {
        throw std::runtime_error("The level-ice lattice contains no particles.");
    }

    for (int rowIndex = 0; rowIndex < static_cast<int>(rowParticleIndices.size()); ++rowIndex)
    {
        const std::vector<int>& currentRow = rowParticleIndices[rowIndex];
        for (int columnIndex = 0; columnIndex + 1 < static_cast<int>(currentRow.size()); ++columnIndex)
        {
            connections_.push_back({currentRow[columnIndex], currentRow[columnIndex + 1]});
        }
        if (rowIndex == 0)
        {
            continue;
        }

        const std::vector<int>& previousRow = rowParticleIndices[rowIndex - 1];
        const bool currentRowIsOdd = rowIndex % 2 == 1;
        for (int columnIndex = 0; columnIndex < static_cast<int>(currentRow.size()); ++columnIndex)
        {
            if (columnIndex < static_cast<int>(previousRow.size()))
            {
                connections_.push_back({currentRow[columnIndex], previousRow[columnIndex]});
            }
            const int secondPreviousIndex = currentRowIsOdd ? columnIndex + 1 : columnIndex - 1;
            if (secondPreviousIndex >= 0 && secondPreviousIndex < static_cast<int>(previousRow.size()))
            {
                connections_.push_back({currentRow[columnIndex], previousRow[secondPreviousIndex]});
            }
        }
    }

    minimumBoundary_ = particlePositions_.front() - Vec3{radius};
    maximumBoundary_ = particlePositions_.front() + Vec3{radius};
    for (const Vec3& position : particlePositions_)
    {
        minimumBoundary_ = math::componentMin(minimumBoundary_, position - Vec3{radius});
        maximumBoundary_ = math::componentMax(maximumBoundary_, position + Vec3{radius});
    }
}

levelIce::Real levelIce::elementMass() const noexcept
{
    const Real radius = elementRadius();
    return (4.0 / 3.0) * math::pi * radius * radius * radius * density_;
}

levelIce::Real levelIce::modeICriticalEnergyReleaseRate(Real flexuralStrength, Real youngsModulus) const noexcept
{
    return positiveFinite(flexuralStrength) && positiveFinite(youngsModulus) ? flexuralStrength * flexuralStrength * thickness_ / youngsModulus : 0.0;
}

levelIce::Real levelIce::normalStiffness(Real youngsModulus) const noexcept
{
    const Real radius = elementRadius();
    return positiveFinite(youngsModulus) ? youngsModulus * math::pi * radius * radius / thickness_ : 0.0;
}

levelIce::Real levelIce::shearStiffness(Real youngsModulus) const noexcept
{
    const Real radius = elementRadius();
    const Real area = math::pi * radius * radius;
    const Real secondMoment = 0.25 * area * radius * radius;
    return positiveFinite(youngsModulus) ? 12.0 * youngsModulus * secondMoment / (thickness_ * thickness_ * thickness_) : 0.0;
}

levelIce::Real levelIce::bendingStiffness(Real youngsModulus) const noexcept
{
    const Real radius = elementRadius();
    const Real secondMoment = 0.25 * math::pi * radius * radius * radius * radius;
    return positiveFinite(youngsModulus) ? youngsModulus * secondMoment / thickness_ : 0.0;
}

levelIce::Real levelIce::torsionalStiffness(Real youngsModulus, Real poissonRatio) const noexcept
{
    if (!positiveFinite(youngsModulus) || !math::isFinite(poissonRatio) || poissonRatio <= -1.0)
    {
        return 0.0;
    }
    const Real radius = elementRadius();
    const Real polarMoment = 0.5 * math::pi * radius * radius * radius * radius;
    const Real shearModulus = youngsModulus / (2.0 * (1.0 + poissonRatio));
    return shearModulus * polarMoment / thickness_;
}

levelIce::Real levelIce::equilibriumCenterZ(Real waterDensity, Real waterLevel) const noexcept
{
    const Real radius = elementRadius();
    const Real densityRatio = density_ / waterDensity;
    if (densityRatio <= 0.0)
    {
        return waterLevel + radius;
    }
    if (densityRatio >= 1.0)
    {
        return waterLevel - radius;
    }

    const Real targetVolume = densityRatio * (4.0 / 3.0) * math::pi * radius * radius * radius;
    Real lowerHeight = 0.0;
    Real upperHeight = 2.0 * radius;
    for (int iteration = 0; iteration < 50; ++iteration)
    {
        const Real height = 0.5 * (lowerHeight + upperHeight);
        const Real submergedVolume = math::pi * height * height * (radius - height / 3.0);
        if (submergedVolume < targetVolume)
        {
            lowerHeight = height;
        }
        else
        {
            upperHeight = height;
        }
    }
    return waterLevel + radius - 0.5 * (lowerHeight + upperHeight);
}

structureIceLoadSolver::structureIceLoadSolver(int device) : SphereDEM(device), loadAccumulator_(std::make_unique<loadAccumulator>()) {}

structureIceLoadSolver::structureIceLoadSolver(executionMode mode, int device) : SphereDEM(mode, device), loadAccumulator_(std::make_unique<loadAccumulator>()) {}

structureIceLoadSolver::~structureIceLoadSolver() = default;

void structureIceLoadSolver::setWaterCondition(const Vec3& velocity, Real density, Real dragCoefficient, Real waterLevel)
{
    if (!math::isFinite(velocity) || !positiveFinite(density) || !math::isFinite(dragCoefficient) || dragCoefficient < 0.0 || !math::isFinite(waterLevel))
    {
        throw std::invalid_argument("The water condition must contain finite velocity/level, positive density, and non-negative drag.");
    }
    waterVelocity_ = velocity;
    waterDensity_ = density;
    waterLevel_ = waterLevel;
    dragCoefficient_ = dragCoefficient;
}

void structureIceLoadSolver::initialize(cpu::mode mode)
{
    initializationAssembly_ = true;
    try
    {
        SphereDEM::initialize(mode);
    }
    catch (...)
    {
        initializationAssembly_ = false;
        throw;
    }
    initializationAssembly_ = false;
    finishInitialization(false, nullptr);
}

void structureIceLoadSolver::initialize(gpu::mode mode, cudaStream_t stream)
{
    initializationAssembly_ = true;
    try
    {
        SphereDEM::initialize(mode, stream);
    }
    catch (...)
    {
        initializationAssembly_ = false;
        throw;
    }
    initializationAssembly_ = false;
    finishInitialization(true, stream);
}

void structureIceLoadSolver::initialize(hybrid::mode mode, cudaStream_t stream)
{
    initializationAssembly_ = true;
    try
    {
        SphereDEM::initialize(mode, stream);
    }
    catch (...)
    {
        initializationAssembly_ = false;
        throw;
    }
    initializationAssembly_ = false;
    finishInitialization(false, stream);
}

void structureIceLoadSolver::finishInitialization(bool useDevice, cudaStream_t stream) { ensureLoadAccumulator(static_cast<int>(LSParticles().hostSize()), useDevice, stream); }

void structureIceLoadSolver::addSphereExternalForceAndTorque(particleContainer::host_container_type& particles)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (particles.size() >= 1024)
#endif
    for (int particleIndex = 0; particleIndex < static_cast<int>(particles.size()); ++particleIndex)
    {
        particle& value = particles[particleIndex];
        Vec3 force = value.force();
        Vec3 torque = value.torque();
        execution::addSphereBuoyancyAndDrag(force,
                                            torque,
                                            value.position(),
                                            value.velocity(),
                                            value.angularVelocity(),
                                            value.radius(),
                                            gravity(),
                                            waterVelocity_,
                                            waterDensity_,
                                            waterLevel_,
                                            dragCoefficient_);
        value.setForce(force);
        value.setTorque(torque);
    }
}

void structureIceLoadSolver::addSphereExternalForceAndTorque(particle::device_type particles, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    cuda::launchSphereBuoyancyAndDrag(particles, gravity(), waterVelocity_, waterDensity_, waterLevel_, dragCoefficient_, stream);
#else
    (void)particles;
    (void)stream;
#endif
}

void structureIceLoadSolver::addLSParticleExternalForceAndTorque(LSParticleContainer::host_container_type& structures)
{
    if (initializationAssembly_)
    {
        return;
    }
    ensureLoadAccumulator(static_cast<int>(structures.size()), false, nullptr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (structures.size() >= 1024)
#endif
    for (int structureIndex = 0; structureIndex < static_cast<int>(structures.size()); ++structureIndex)
    {
        loadAccumulator_->hostForces_[structureIndex] += structures[structureIndex].force();
    }
    ++loadAccumulator_->sampleCount_;
}

void structureIceLoadSolver::addLSParticleExternalForceAndTorque(LSParticle::device_type structures, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (initializationAssembly_)
    {
        return;
    }
    ensureLoadAccumulator(structures.size_, true, stream);
    cuda::launchStructureForceAccumulation(loadAccumulator_->deviceForces_.data(), structures, stream);
    ++loadAccumulator_->sampleCount_;
#else
    (void)structures;
    (void)stream;
#endif
}

void structureIceLoadSolver::ensureLoadAccumulator(int structureCount, bool useDevice, cudaStream_t stream)
{
    loadAccumulator& accumulator = *loadAccumulator_;
    if (accumulator.structureCount_ == structureCount && accumulator.useDevice_ == useDevice)
    {
        return;
    }

    accumulator.hostForces_.assign(static_cast<std::size_t>(structureCount), Vec3::zero());
    accumulator.structureCount_ = structureCount;
    accumulator.sampleCount_ = 0;
    accumulator.useDevice_ = useDevice;

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (useDevice)
    {
        accumulator.deviceForces_.allocate(static_cast<std::size_t>(structureCount));
        accumulator.deviceCount_ = structureCount;
        cuda::launchStructureForceClear(accumulator.deviceForces_.data(), structureCount, stream);
    }
    else
    {
        accumulator.deviceForces_.reset();
        accumulator.deviceCount_ = 0;
    }
#else
    (void)stream;
    accumulator.deviceCount_ = 0;
#endif
}

void structureIceLoadSolver::writeSystemVTU(int frameIndex)
{
    SphereDEM::writeSystemVTU(frameIndex);
    writeStructureIceLoad(frameIndex);
}

void structureIceLoadSolver::writeStructureIceLoad(int frameIndex)
{
    const bool useDevice = mode() == executionMode::GPU;
    ensureLoadAccumulator(static_cast<int>(LSParticles().hostSize()), useDevice, gpuStream());
    loadAccumulator& accumulator = *loadAccumulator_;

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (useDevice && accumulator.structureCount_ > 0)
    {
        host_device_detail::checkCuda(cudaMemcpyAsync(accumulator.hostForces_.data(),
                                                      accumulator.deviceForces_.data(),
                                                      static_cast<std::size_t>(accumulator.structureCount_) * sizeof(Vec3),
                                                      cudaMemcpyDeviceToHost,
                                                      gpuStream()),
                                      "copy structure force accumulator to host");
        host_device_detail::synchronize(gpuStream());
    }
#endif

    const std::filesystem::path fileName = std::filesystem::path(outputDirectory()) / "structureIceLoad.dat";
    std::ofstream output(fileName, frameIndex == 0 ? std::ios::trunc : std::ios::app);
    if (!output)
    {
        throw std::runtime_error("Cannot open structure ice-load DAT file: " + fileName.string());
    }

    if (frameIndex == 0)
    {
        output << "# time";
        for (int structureIndex = 0; structureIndex < accumulator.structureCount_; ++structureIndex)
        {
            output << " Structure" << structureIndex << "_Fx(MN)" << " Structure" << structureIndex << "_Fy(MN)" << " Structure" << structureIndex << "_Fz(MN)";
        }
        output << '\n';
    }

    output << std::scientific << std::setprecision(16) << time();
    const Real inverseSampleCount = accumulator.sampleCount_ > 0 ? 1.0 / accumulator.sampleCount_ : 0.0;
    for (const Vec3& forceSum : accumulator.hostForces_)
    {
        const Vec3 averageForce = forceToMegaNewtons * inverseSampleCount * forceSum;
        output << ' ' << averageForce.x << ' ' << averageForce.y << ' ' << averageForce.z;
    }
    output << '\n';
    if (!output)
    {
        throw std::runtime_error("Failed while writing structure ice-load DAT file: " + fileName.string());
    }

    std::fill(accumulator.hostForces_.begin(), accumulator.hostForces_.end(), Vec3::zero());
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (useDevice)
    {
        cuda::launchStructureForceClear(accumulator.deviceForces_.data(), accumulator.structureCount_, gpuStream());
    }
#endif
    accumulator.sampleCount_ = 0;
}

double structureIceLoadSolver::systemDeviceMemoryGB() const noexcept
{
    const double accumulatorMemory = static_cast<double>(loadAccumulator_->deviceCount_) * sizeof(Vec3) / (1024.0 * 1024.0 * 1024.0);
    return SphereDEM::systemDeviceMemoryGB() + accumulatorMemory;
}

} // namespace fundem
