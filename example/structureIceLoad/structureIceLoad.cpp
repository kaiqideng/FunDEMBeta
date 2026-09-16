/**
 * @file structureIceLoad.cpp
 * @brief Demonstrates bonded level ice driven against a fixed LS structure.
 */
#include "structureIceLoadSolver.h"

#include "data/myLSObject.h"
#include "interaction/bond.h"
#include "material/LSMaterial.h"
#include "material/material.h"
#include "particle/LSParticle.h"
#include "particle/particle.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

using fundem::executionMode;
using fundem::vtuFormat;
using fundem::math::Real;
using fundem::math::Vec3;

struct options {
    std::string outputDirectory_{"structureIceLoad_files"};
    int stepCount_{-1};
    int device_{0};
    executionMode mode_{executionMode::CPU};
    vtuFormat vtuFormat_{vtuFormat::binary};
};

executionMode parseExecutionMode(const std::string& value)
{
    if (value == "cpu")
    {
        return executionMode::CPU;
    }
    if (value == "gpu")
    {
        return executionMode::GPU;
    }
    if (value == "hybrid")
    {
        return executionMode::Hybrid;
    }
    throw std::invalid_argument("--mode must be cpu, gpu, or hybrid.");
}

options parseOptions(int argc, char** argv)
{
    options result;
    for (int argumentIndex = 1; argumentIndex < argc; ++argumentIndex)
    {
        const std::string argument = argv[argumentIndex];
        if (argument == "--ascii")
        {
            result.vtuFormat_ = vtuFormat::ascii;
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
        else if (argument == "--mode")
        {
            if (++argumentIndex >= argc)
            {
                throw std::invalid_argument("--mode requires cpu, gpu, or hybrid.");
            }
            result.mode_ = parseExecutionMode(argv[argumentIndex]);
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

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;

        const options run = parseOptions(argc, argv);
        constexpr Real gravityMagnitude = 9.81;
        constexpr Real waterDensity = 1000.0;
        constexpr Real waterLevel = 0.0;
        constexpr Real iceThickness = 0.18;
        constexpr Real iceDensity = 910.0;
        constexpr Real iceYoungsModulus = 2.0e9;
        constexpr Real icePoissonRatio = 0.3;
        constexpr Real iceFlexuralStrength = 0.5e6;
        constexpr Real iceMinimumX = -1.8;
        constexpr Real iceMaximumX = 0.0;
        constexpr Real iceMinimumY = -0.6;
        constexpr Real iceMaximumY = 0.6;
        constexpr Real structureRadius = 0.3;
        constexpr Real duration = 2.0;
        constexpr Real outputInterval = 0.05;
        constexpr Vec3 waterVelocity{0.2, 0.0, 0.0};

        structureIceLoadSolver simulation{run.mode_, run.device_};
        simulation.setWaterCondition(waterVelocity, waterDensity, 0.1, waterLevel);
        simulation.ice().setThickness(iceThickness);
        simulation.ice().setDensity(iceDensity);
        simulation.ice().build(iceMinimumX, iceMaximumX, iceMinimumY, iceMaximumY, waterDensity, waterLevel);

        const levelIce& ice = simulation.ice();
        const Real normalStiffness = ice.normalStiffness(iceYoungsModulus);
        const Real shearStiffness = ice.shearStiffness(iceYoungsModulus);
        const Real bendingStiffness = ice.bendingStiffness(iceYoungsModulus);
        const Real torsionalStiffness = ice.torsionalStiffness(iceYoungsModulus, icePoissonRatio);
        const Real modeICriticalEnergy = ice.modeICriticalEnergyReleaseRate(iceFlexuralStrength, iceYoungsModulus);
        const int iceMaterialIndex = simulation.addMaterial(material{normalStiffness, shearStiffness, 0.0, 0.0, 0.1, 0.0, 0.0, 0.8, iceDensity});

        const Real fixedBoundaryWidth = 2.0 * ice.elementRadius();
        for (const Vec3& position : ice.particlePositions())
        {
            particle iceElement;
            iceElement.setPosition(position);
            iceElement.setVelocity(waterVelocity);
            iceElement.setRadius(ice.elementRadius());
            iceElement.setMaterial(simulation.materials(), iceMaterialIndex);
            if (position.x < ice.minimumBoundary().x + fixedBoundaryWidth || position.y < ice.minimumBoundary().y + fixedBoundaryWidth || position.y > ice.maximumBoundary().y - fixedBoundaryWidth)
            {
                iceElement.setInfiniteMass();
            }
            simulation.addSphere(iceElement);
        }

        const Real crossSectionArea = math::pi * ice.elementRadius() * ice.elementRadius();
        for (const levelIceConnection& indices : ice.connections())
        {
            const Vec3 direction = ice.particlePositions()[indices.firstParticleIndex_] - ice.particlePositions()[indices.secondParticleIndex_];
            bond connection{math::norm(direction)};
            if (!connection.setConnection(simulation.spheres(), indices.firstParticleIndex_, indices.secondParticleIndex_, math::normalizedOrZero(direction)) ||
                !connection.setStiffness(normalStiffness, shearStiffness, bendingStiffness, torsionalStiffness) || !connection.setCrossSectionArea(crossSectionArea) ||
                !connection.setModeICriticalEnergy(modeICriticalEnergy) || !connection.setModeIICriticalEnergy(2.0 * modeICriticalEnergy))
            {
                throw std::runtime_error("Failed to create a level-ice bond.");
            }
            simulation.addBond(connection);
        }

        levelset::ConeWall structureGeometry{{0.0, 0.0, -structureRadius}, {0.0, 0.0, structureRadius}, 1.5 * structureRadius, 0.5 * structureRadius};
        structureGeometry.buildLSGrid(0.25 * iceThickness, 3);
        const int structureGeometryIndex = simulation.addGeometry(structureGeometry, true);
        const int structureMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, 0.1, 1.0, 1.0});

        const Real structureCenterX = ice.maximumBoundary().x + structureRadius + ice.elementRadius();
        LSParticle structure;
        structure.setPosition({structureCenterX, 0.0, 0.0});
        structure.setMaterial(simulation.materials(), structureMaterialIndex);
        structure.setGeometry(simulation.geometries(), structureGeometryIndex);
        simulation.addLSParticle(structure);

        const Real timeStep = (math::pi / 50.0) * std::sqrt(ice.elementMass() / normalStiffness);
        const int stepCount = run.stepCount_ >= 0 ? run.stepCount_ : static_cast<int>(std::ceil(duration / timeStep));
        simulation.setBoundary({ice.minimumBoundary().x - iceThickness, ice.minimumBoundary().y - iceThickness, -2.0 * structureRadius},
                               {structureCenterX + 2.0 * structureRadius, ice.maximumBoundary().y + iceThickness, 2.0 * structureRadius});
        simulation.setGravity({0.0, 0.0, -gravityMagnitude});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory(run.outputDirectory_);
        simulation.setOutputStepInterval(std::max(1, static_cast<int>(std::lround(outputInterval / timeStep))));
        simulation.setVTUOutputFormat(run.vtuFormat_);
        simulation.solve(stepCount);
    }
    catch (const std::exception& error)
    {
        std::cerr << "[structureIceLoad] " << error.what() << '\n';
        return 1;
    }
    return 0;
}
