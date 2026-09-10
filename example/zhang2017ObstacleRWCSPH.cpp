#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "solver/SPHDEM.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

using fundem::SPHDEM;
using fundem::math::Real;
using fundem::math::Vec3;

struct options {
    std::string outputDirectory_{"zhang2017ObstacleRWCSPH_files"};
    int stepCount_{300000};
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

/** Three-dimensional dam break against a cuboidal obstacle from Zhang et al. (2017). */
int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;

        const options run = parseOptions(argc, argv);
        constexpr Real gravityMagnitude = 9.81;
        constexpr Real tankLength = 3.20;
        constexpr Real tankWidth = 1.00;
        constexpr Real tankHeight = 2.00;
        constexpr Real waterLength = 1.23;
        constexpr Real waterHeight = 0.55;
        constexpr Real spacing = waterHeight / 44.0;
        constexpr Real smoothingLength = 1.3 * spacing;
        constexpr Real referenceDensity = 1000.0;
        constexpr Real dynamicViscosity = 1.0e-3;
        constexpr Real obstacleLength = 0.1625;
        constexpr Real obstacleWidth = 0.40;
        constexpr Real obstacleHeight = 0.1625;
        constexpr Real obstacleCenterX = waterLength + 1.032;
        constexpr Real timeStep = 1.0e-5;

        SPHDEM simulation(run.device_);
        simulation.setSPHProperties(spacing, smoothingLength, referenceDensity, dynamicViscosity);
        simulation.setSPHMaximumVelocity(std::sqrt(2.0 * gravityMagnitude * waterHeight));

        const int3 waterParticleCount{static_cast<int>(waterLength / spacing), static_cast<int>(tankWidth / spacing), static_cast<int>(waterHeight / spacing)};
        simulation.addSPHBlock(Vec3::zero(), waterParticleCount);

        const int wallMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, 0.0, 1.0, 1.0});

        levelset::BoxWall tankWall{{tankLength, tankWidth, tankHeight}};
        tankWall.buildLSGrid(spacing, 3);
        tankWall.reverseSDFSign();
        const int tankGeometryIndex = simulation.addGeometry(tankWall, true);

        LSParticle tank;
        tank.setPosition({0.5 * tankLength, 0.5 * tankWidth, 0.5 * tankHeight});
        tank.setMaterial(simulation.materials(), wallMaterialIndex);
        tank.setGeometry(simulation.geometries(), tankGeometryIndex);
        simulation.addLSParticle(tank);

        levelset::BoxWall obstacle{{obstacleLength, obstacleWidth, obstacleHeight}};
        obstacle.buildLSGrid(spacing, 3);
        const int obstacleGeometryIndex = simulation.addGeometry(obstacle, true);

        LSParticle obstacleParticle;
        obstacleParticle.setPosition({obstacleCenterX, 0.5 * tankWidth, 0.5 * obstacleHeight});
        obstacleParticle.setMaterial(simulation.materials(), wallMaterialIndex);
        obstacleParticle.setGeometry(simulation.geometries(), obstacleGeometryIndex);
        simulation.addLSParticle(obstacleParticle);

        simulation.setBoundary(Vec3::zero(), {tankLength, tankWidth, tankHeight});
        simulation.setGravity({0.0, 0.0, -gravityMagnitude});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory(run.outputDirectory_);
        simulation.setOutputStepInterval(static_cast<int>(0.05 / timeStep));
        simulation.setVTUOutputFormat(run.vtuFormat_);
        simulation.setExecutionMode(executionMode::Hybrid);
        simulation.solve(run.stepCount_);
    }
    catch (const std::exception& error)
    {
        std::cerr << "[zhang2017ObstacleRWCSPH] " << error.what() << '\n';
        return 1;
    }
    return 0;
}
