/**
 * @file tutorial3.cpp
 * @brief Reproduces the classical three-dimensional bore-in-a-box experiment.
 * @details Geometry follows Gomez-Gesteira and Dalrymple, J. Waterway, Port,
 * Coastal, Ocean Eng. 130(2), 63-69 (2004).
 */
#include "data/datWriter.h"
#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "solver/SPHDEM.h"
#include "tutorialOptions.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

/** Records the measured square-column reaction without adding any external load. */
class BoreInBoxSimulation final : public fundem::SPHDEM
{
public:
    BoreInBoxSimulation(fundem::executionMode mode, int device) : SPHDEM(mode, device) {}

    void setColumnParticleIndex(int value) noexcept { columnParticleIndex_ = value; }

    /** Allocates stable host storage for the initial force and every completed DEM step. */
    void prepareColumnForceHistory(int numberOfSteps)
    {
        if (numberOfSteps < 0)
        {
            throw std::invalid_argument("The force-history step count cannot be negative.");
        }
        forceSamples_.clear();
        forceSamples_.reserve(static_cast<std::size_t>(numberOfSteps) + 1);
        expectedForceSampleCount_ = numberOfSteps + 1;
    }

    /** Writes the raw history and the paper-compatible centered triangular moving average. */
    void writeColumnForceHistory(fundem::math::Real filterWindow) const
    {
        using fundem::math::Real;
        using fundem::math::Vec3;

        if (filterWindow <= fundem::math::defaultTolerance || forceSamples_.empty())
        {
            throw std::logic_error("Column-force filtering requires a positive window and at least one DEM-step sample.");
        }
        if (static_cast<int>(forceSamples_.size()) != expectedForceSampleCount_)
        {
            throw std::logic_error("The column-force history is incomplete.");
        }

        fundem::datWriter rawWriter((std::filesystem::path(outputDirectory()) / "columnForce.dat").string(), {"time", "forceX", "forceY", "forceZ", "horizontalForce"}, false);
        for (const ForceSample& sample : forceSamples_)
        {
            rawWriter.appendRow({sample.time_, sample.force_.x, sample.force_.y, sample.force_.z, std::hypot(sample.force_.x, sample.force_.y)});
        }
        rawWriter.flush();

        const int sampleCount = static_cast<int>(forceSamples_.size());
        std::vector<Real> sampleWidths(static_cast<std::size_t>(sampleCount), filterWindow);
        if (sampleCount > 1)
        {
            for (int sampleIndex = 1; sampleIndex < sampleCount; ++sampleIndex)
            {
                if (forceSamples_[sampleIndex].time_ <= forceSamples_[sampleIndex - 1].time_)
                {
                    throw std::logic_error("Column-force sample times must increase strictly.");
                }
            }
            sampleWidths.front() = forceSamples_[1].time_ - forceSamples_[0].time_;
            sampleWidths.back() = forceSamples_.back().time_ - forceSamples_[sampleCount - 2].time_;
            for (int sampleIndex = 1; sampleIndex + 1 < sampleCount; ++sampleIndex)
            {
                sampleWidths[sampleIndex] = 0.5 * (forceSamples_[sampleIndex + 1].time_ - forceSamples_[sampleIndex - 1].time_);
            }
        }

        fundem::datWriter filteredWriter((std::filesystem::path(outputDirectory()) / "columnForceFiltered.dat").string(), {"time", "forceX", "forceY", "forceZ", "horizontalForce"}, false);
        const Real halfWindow = 0.5 * filterWindow;
        int windowBegin = 0;
        int windowEnd = 0;
        for (int sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex)
        {
            const Real sampleTime = forceSamples_[sampleIndex].time_;
            while (windowBegin < sampleIndex && sampleTime - forceSamples_[windowBegin].time_ >= halfWindow)
            {
                ++windowBegin;
            }
            if (windowEnd < sampleIndex)
            {
                windowEnd = sampleIndex;
            }
            while (windowEnd + 1 < sampleCount && forceSamples_[windowEnd + 1].time_ - sampleTime < halfWindow)
            {
                ++windowEnd;
            }

            Vec3 filteredForce = Vec3::zero();
            Real totalWeight = 0.0;
            for (int neighborIndex = windowBegin; neighborIndex <= windowEnd; ++neighborIndex)
            {
                const Real normalizedDistance = std::abs(forceSamples_[neighborIndex].time_ - sampleTime) / halfWindow;
                const Real weight = (1.0 - normalizedDistance) * sampleWidths[neighborIndex];
                if (weight > 0.0)
                {
                    filteredForce += weight * forceSamples_[neighborIndex].force_;
                    totalWeight += weight;
                }
            }
            if (totalWeight <= fundem::math::defaultTolerance)
            {
                throw std::logic_error("The column-force filter produced an empty averaging window.");
            }
            filteredForce /= totalWeight;
            filteredWriter.appendRow({sampleTime, filteredForce.x, filteredForce.y, filteredForce.z, std::hypot(filteredForce.x, filteredForce.y)});
        }
        filteredWriter.flush();
    }

protected:
    void addLSParticleExternalForceAndTorque(fundem::LSParticleContainer::host_container_type& particles) override
    {
        validateColumnParticleIndex(static_cast<int>(particles.size()));
        recordColumnForce(nextColumnForceSampleTime(), particles[columnParticleIndex_].force());
    }

    void addLSParticleExternalForceAndTorque(fundem::LSParticle::device_type particles, cudaStream_t stream) override
    {
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        validateColumnParticleIndex(static_cast<int>(LSParticles().hostSize()));
        recordColumnForce(nextColumnForceSampleTime(), fundem::math::Vec3::zero());
        fundem::host_device_detail::checkCuda(cudaMemcpyAsync(&forceSamples_.back().force_, particles.force_ + columnParticleIndex_, sizeof(fundem::math::Vec3), cudaMemcpyDeviceToHost, stream),
                                              "cudaMemcpyAsync tutorial3 column force");
#else
        (void)particles;
        (void)stream;
#endif
    }

private:
    struct ForceSample {
        fundem::math::Real time_;
        fundem::math::Vec3 force_;
    };

    void validateColumnParticleIndex(int particleCount) const
    {
        if (columnParticleIndex_ < 0 || columnParticleIndex_ >= particleCount)
        {
            throw std::logic_error("The bore-in-a-box column particle has not been configured.");
        }
    }

    fundem::math::Real nextColumnForceSampleTime() const noexcept { return forceSamples_.empty() ? time() : time() + timeStep(); }

    void recordColumnForce(fundem::math::Real sampleTime, const fundem::math::Vec3& force)
    {
        if (expectedForceSampleCount_ <= 0)
        {
            throw std::logic_error("Prepare the column-force history before solving.");
        }
        if (static_cast<int>(forceSamples_.size()) >= expectedForceSampleCount_)
        {
            throw std::logic_error("The column-force history received more samples than reserved.");
        }
        forceSamples_.push_back({sampleTime, force});
    }

    int columnParticleIndex_{-1};           ///< Stable LSParticle index of the instrumented column.
    int expectedForceSampleCount_{0};       ///< Initial state plus one sample per requested DEM step.
    std::vector<ForceSample> forceSamples_; ///< Stable storage also used as asynchronous CUDA destinations.
};

} // namespace

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;

        constexpr math::Real gravityMagnitude = 9.81;
        constexpr math::Real referenceDensity = 1000.0;
        constexpr math::Real dynamicViscosity = 1.0e-3;
        constexpr math::Real spacing = 0.01;
        constexpr math::Real smoothingLength = 1.3 * spacing;
        constexpr math::Real timeStep = 2.0e-5;
        constexpr math::Real duration = 3.0;
        constexpr math::Real outputTimeInterval = 0.02;
        constexpr math::Real forceFilterWindow = 0.013;

        constexpr math::Vec3 tankSize{1.60, 0.61, 0.75};
        constexpr int3 reservoirParticleCount{40, 61, 30};
        constexpr math::Real initialWaterDepth = reservoirParticleCount.z * spacing;
        constexpr math::Vec3 columnSize{0.12, 0.12, 0.75};
        constexpr math::Vec3 columnCenter{0.96, 0.30, 0.375};

        const tutorial::TutorialOptions run = tutorial::parseTutorialOptions(argc, argv, "tutorial3_files");
        BoreInBoxSimulation simulation(run.mode_, run.device_);
        simulation.setSPHProperties(spacing, smoothingLength, referenceDensity, dynamicViscosity);
        simulation.setSPHMaximumVelocity(std::sqrt(2.0 * gravityMagnitude * initialWaterDepth));

        // The 0.40 m reservoir and 0.01 m downstream wet bed reproduce the experiment.
        simulation.addSPHBlock(math::Vec3::zero(), reservoirParticleCount);
        simulation.addSPHBlock({0.40, 0.0, 0.0}, {50, 61, 1});
        simulation.addSPHBlock({1.02, 0.0, 0.0}, {58, 61, 1});
        simulation.addSPHBlock({0.90, 0.0, 0.0}, {12, 24, 1});
        simulation.addSPHBlock({0.90, 0.36, 0.0}, {12, 25, 1});

        const int wallMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, 0.0, 1.0, 1.0});

        // Reversing a closed box places the boundary solid outside the fluid cavity.
        levelset::BoxWall tank{tankSize};
        tank.buildLSGrid(spacing, 3);
        tank.reverseSDFSign();
        const int tankGeometryIndex = simulation.addGeometry(tank, true);
        LSParticle fixedTank;
        fixedTank.setPosition(0.5 * tankSize);
        fixedTank.setMaterial(simulation.materials(), wallMaterialIndex);
        fixedTank.setGeometry(simulation.geometries(), tankGeometryIndex);
        simulation.addLSParticle(fixedTank);

        levelset::BoxWall squareColumn{columnSize};
        squareColumn.buildLSGrid(spacing, 3);
        const int columnGeometryIndex = simulation.addGeometry(squareColumn, true);
        LSParticle fixedColumn;
        fixedColumn.setPosition(columnCenter);
        fixedColumn.setMaterial(simulation.materials(), wallMaterialIndex);
        fixedColumn.setGeometry(simulation.geometries(), columnGeometryIndex);
        const int columnParticleIndex = simulation.addLSParticle(fixedColumn);
        simulation.setColumnParticleIndex(columnParticleIndex);

        // Extra vertical search space lets spray leave the open 0.75 m-high tank.
        simulation.setBoundary(math::Vec3::zero(), {tankSize.x, tankSize.y, 1.0});
        simulation.setGravity({0.0, 0.0, -gravityMagnitude});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory(run.outputDirectory_);
        simulation.setOutputStepInterval(std::max(1, static_cast<int>(std::lround(outputTimeInterval / timeStep))));
        simulation.setVTUOutputFormat(run.vtuFormat_);

        const int stepCount = run.stepCount_ >= 0 ? run.stepCount_ : static_cast<int>(std::ceil(duration / timeStep));
        simulation.prepareColumnForceHistory(stepCount);
        simulation.solve(stepCount);
        simulation.writeColumnForceHistory(forceFilterWindow);
    }
    catch (const std::exception& error)
    {
        std::cerr << "tutorial3 failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
