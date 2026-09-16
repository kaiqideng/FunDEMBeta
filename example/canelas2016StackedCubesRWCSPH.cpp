#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "solver/SPHDEM.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

using fundem::math::Real;

struct options {
    std::string outputDirectory_{"canelas2016StackedCubesRWCSPH_files"};
    int stepCount_{-1};
    int device_{0};
    fundem::vtuFormat vtuFormat_{fundem::vtuFormat::binary};
};

options parseOptions(int argc, char** argv)
{
    options result;
    for (int argumentIndex = 1; argumentIndex < argc; ++argumentIndex)
    {
        const std::string argument = argv[argumentIndex];
        if (argument == "--ascii")
        {
            result.vtuFormat_ = fundem::vtuFormat::ascii;
        }
        else if (argument == "--steps")
        {
            if (++argumentIndex >= argc)
            {
                throw std::invalid_argument("--steps requires a non-negative integer.");
            }
            result.stepCount_ = std::stoi(argv[argumentIndex]);
            if (result.stepCount_ < 0)
            {
                throw std::invalid_argument("--steps requires a non-negative integer.");
            }
        }
        else if (argument == "--device")
        {
            if (++argumentIndex >= argc)
            {
                throw std::invalid_argument("--device requires a non-negative integer.");
            }
            result.device_ = std::stoi(argv[argumentIndex]);
            if (result.device_ < 0)
            {
                throw std::invalid_argument("--device requires a non-negative integer.");
            }
        }
        else if (argument == "--output")
        {
            if (++argumentIndex >= argc)
            {
                throw std::invalid_argument("--output requires a directory.");
            }
            result.outputDirectory_ = argv[argumentIndex];
        }
        else
        {
            throw std::invalid_argument("Unknown option: " + argument);
        }
    }
    return result;
}

} // namespace

/**
 * Dam-break flow impacting three vertically stacked PVC cubes from Canelas et
 * al. (2016), Section 5.2.2.
 *
 * R. B. Canelas, R. M. L. Ferreira, A. J. C. Crespo and M. Dominguez,
 * "SPH-DCDEM model for arbitrary geometries in free surface solid-fluid flows,"
 * Computer Physics Communications 202 (2016), 131-140.
 */
int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;

        const options run = parseOptions(argc, argv);
        constexpr Real gravityMagnitude = 9.81;
        constexpr Real referenceDensity = 1000.0;
        constexpr Real dynamicViscosity = 1.0e-3;
        constexpr Real flumeLength = 8.00;
        constexpr Real downstreamLength = 3.50;
        constexpr Real upstreamLength = flumeLength - downstreamLength;
        constexpr Real flumeWidth = 0.70;
        constexpr Real flumeHeight = 0.70;
        constexpr Real initialWaterDepth = 0.40;
        constexpr Real cubeSide = 0.15;
        constexpr Real cubeDensity = 800.0;
        constexpr Real cubeYoungModulus = 3.0e9;
        constexpr Real cubeFrictionCoefficient = 0.0;  // 0.45
        constexpr Real flumeFrictionCoefficient = 0.0; // 0.35
        constexpr Real restitutionCoefficient = 1.0;
        constexpr Real cubeFrontX = 1.70;
        constexpr Real spacing = 0.0125;
        constexpr Real smoothingLength = 1.3 * spacing;
        constexpr Real endTime = 2.0;
        constexpr Real outputInterval = 0.05;

        const Real normalStiffnessPerUnitArea = cubeYoungModulus / cubeSide;
        const Real shearStiffnessPerUnitArea = (2.0 / 7.0) * normalStiffnessPerUnitArea;
        const Real cubeMass = cubeDensity * cubeSide * cubeSide * cubeSide;
        const Real restitutionLogarithm = std::log(restitutionCoefficient);
        const Real dampingRatio = -restitutionLogarithm / std::sqrt(restitutionLogarithm * restitutionLogarithm + math::pi * math::pi);
        const Real contactFrequency = std::sqrt(normalStiffnessPerUnitArea * cubeSide * cubeSide * (1.0 - dampingRatio * dampingRatio) / cubeMass);
        const Real timeStep = math::pi / (50.0 * contactFrequency);
        const int stepCount = run.stepCount_ >= 0 ? run.stepCount_ : static_cast<int>(std::ceil(endTime / timeStep));

        SPHDEM simulation(run.device_);
        simulation.setSPHProperties(spacing, smoothingLength, referenceDensity, dynamicViscosity);
        simulation.setSPHMaximumVelocity(std::sqrt(2.0 * gravityMagnitude * initialWaterDepth));

        const int3 waterParticleCount{static_cast<int>(std::ceil(upstreamLength / spacing)),
                                      static_cast<int>(std::ceil(flumeWidth / spacing)),
                                      static_cast<int>(std::ceil(initialWaterDepth / spacing))};
        simulation.addSPHBlock({-upstreamLength, 0.0, 0.0}, waterParticleCount);

        const int wallMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, flumeFrictionCoefficient, restitutionCoefficient, 1.0});
        const int cubeMaterialIndex = simulation.addMaterial(LSMaterial{normalStiffnessPerUnitArea, shearStiffnessPerUnitArea, cubeFrictionCoefficient, restitutionCoefficient, cubeDensity});

        levelset::BoxWall flumeWall{{flumeLength, flumeWidth, flumeHeight}};
        flumeWall.buildLSGrid(spacing, 3, true);
        flumeWall.reverseSDFSign();
        const int wallGeometryIndex = simulation.addGeometry(flumeWall);

        LSParticle wall;
        wall.setPosition({0.5 * (downstreamLength - upstreamLength), 0.5 * flumeWidth, 0.5 * flumeHeight});
        wall.setMaterial(simulation.materials(), wallMaterialIndex);
        wall.setGeometry(simulation.geometries(), wallGeometryIndex);
        simulation.addLSParticle(wall);

        levelset::Superellipsoid cube{0.5 * cubeSide, 0.5 * cubeSide, 0.5 * cubeSide, 0.1, 0.1};
        cube.buildSurfaceNode();
        cube.buildLSGrid(spacing, 3);
        const int cubeGeometryIndex = simulation.addGeometry(cube);
        const Real cubeCenterX = cubeFrontX + 0.5 * cubeSide;
        const Real cubeCenterY = 0.5 * flumeWidth;
        for (int cubeIndex = 0; cubeIndex < 3; ++cubeIndex)
        {
            LSParticle cubeParticle;
            cubeParticle.setPosition({cubeCenterX, cubeCenterY, (0.5 + Real(cubeIndex)) * cubeSide});
            cubeParticle.setMaterial(simulation.materials(), cubeMaterialIndex);
            cubeParticle.setGeometry(simulation.geometries(), cubeGeometryIndex);
            simulation.addLSParticle(cubeParticle);
        }

        simulation.setBoundary({-upstreamLength, 0.0, 0.0}, {downstreamLength, flumeWidth, flumeHeight});
        simulation.setGravity({0.0, 0.0, -gravityMagnitude});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory(run.outputDirectory_);
        simulation.setOutputStepInterval(std::max(1, static_cast<int>(std::round(outputInterval / timeStep))));
        simulation.setVTUOutputFormat(run.vtuFormat_);
        simulation.setExecutionMode(executionMode::Hybrid);
        simulation.solve(stepCount);
    }
    catch (const std::exception& error)
    {
        std::cerr << "[canelas2016StackedCubesRWCSPH] " << error.what() << '\n';
        return 1;
    }
    return 0;
}
