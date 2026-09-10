/** CPU Gömböc self-righting tutorial. */
#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "solver/LSDEM.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

using Real = fundem::math::Real;

constexpr Real defaultTimeStep = 1.0e-4;
constexpr Real defaultDuration = 120.0;
constexpr Real outputTimeInterval = 0.05;
constexpr Real ceramicNormalStiffnessPerUnitArea = 2.4e10;
constexpr Real ceramicShearStiffnessPerUnitArea = 7.1e9;
constexpr Real ceramicFrictionCoefficient = 0.1;
constexpr Real ceramicRestitutionCoefficient = 0.4;
constexpr Real ceramicDensity = 2400.0;
constexpr Real initialRotationX = -87.1886819;
constexpr Real initialRotationY = -52.5462996;
constexpr Real initialRotationZ = -161.4623206;

struct options {
    std::filesystem::path outputDirectory_;
    Real timeStep_{defaultTimeStep};
    Real duration_{defaultDuration};
    int stepCount_{-1};
    fundem::vtuFormat vtuFormat_{fundem::vtuFormat::binary};
};

std::filesystem::path executableDirectory(const char* executable)
{
    std::error_code error;
    std::filesystem::path path = std::filesystem::weakly_canonical(executable, error);
    if (error)
    {
        error.clear();
        path = std::filesystem::absolute(executable, error);
    }
    return error ? std::filesystem::current_path() : path.parent_path();
}

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
        else if (argument == "--time-step")
        {
            if (++argumentIndex >= argc || (result.timeStep_ = std::stod(argv[argumentIndex])) <= 0.0)
            {
                throw std::invalid_argument("--time-step requires a positive value.");
            }
        }
        else if (argument == "--duration")
        {
            if (++argumentIndex >= argc || (result.duration_ = std::stod(argv[argumentIndex])) <= 0.0)
            {
                throw std::invalid_argument("--duration requires a positive value.");
            }
        }
        else if (argument == "--steps")
        {
            if (++argumentIndex >= argc || (result.stepCount_ = std::stoi(argv[argumentIndex])) < 0)
            {
                throw std::invalid_argument("--steps requires a non-negative integer.");
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

fundem::math::Quaternion initialOrientation() noexcept
{
    using fundem::math::Quaternion;
    using fundem::math::Vec3;
    constexpr Real degree = fundem::math::pi / 180.0;
    const Quaternion qx = Quaternion::fromUnitAxisAngle(Vec3::unitX(), initialRotationX * degree);
    const Quaternion qy = Quaternion::fromUnitAxisAngle(Vec3::unitY(), initialRotationY * degree);
    const Quaternion qz = Quaternion::fromUnitAxisAngle(Vec3::unitZ(), initialRotationZ * degree);
    return fundem::math::normalizedOrIdentity(qz * qy * qx);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;

        options run = parseOptions(argc, argv);
        const std::filesystem::path executablePath = executableDirectory(argv[0]);
        if (run.outputDirectory_.empty())
        {
            run.outputDirectory_ = executablePath / "tutorial0_files";
        }
        else if (run.outputDirectory_.is_relative())
        {
            run.outputDirectory_ = executablePath / run.outputDirectory_;
        }

        LSDEM simulation;
        const int gombocMaterialIndex =
            simulation.addMaterial(LSMaterial{ceramicNormalStiffnessPerUnitArea, ceramicShearStiffnessPerUnitArea, ceramicFrictionCoefficient, ceramicRestitutionCoefficient, ceramicDensity});
        const int wallMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, ceramicFrictionCoefficient, ceramicRestitutionCoefficient, 1.0});

        levelset::TriangleMesh gomboc;
        gomboc.loadOBJ((executablePath / "Gomboc.obj").string());
        gomboc.buildLSGrid(64);
        const int gombocGeometryIndex = simulation.addGeometry(gomboc);

        LSParticle fallingGomboc;
        const math::Quaternion orientation = initialOrientation();
        fallingGomboc.setOrientation(orientation);
        fallingGomboc.setMaterial(simulation.materials(), gombocMaterialIndex);
        fallingGomboc.setGeometry(simulation.geometries(), gombocGeometryIndex);
        fallingGomboc.setPosition({0.0, 0.0, 0.1});
        simulation.addLSParticle(fallingGomboc);

        levelset::PlaneWall plane{{0.0, 0.0, 1.0}, 1.0};
        plane.buildLSGrid(0.1, 2);
        const int planeGeometryIndex = simulation.addGeometry(plane, true);

        LSParticle fixedPlane;
        fixedPlane.setMaterial(simulation.materials(), wallMaterialIndex);
        fixedPlane.setGeometry(simulation.geometries(), planeGeometryIndex);
        simulation.addLSParticle(fixedPlane);

        simulation.setBoundary({-0.5, -0.5, 0.0}, {0.5, 0.5, 0.35});
        simulation.setGravity({0.0, 0.0, -9.81});
        simulation.setTimeStep(run.timeStep_);
        simulation.setOutputDirectory(run.outputDirectory_.string());
        simulation.setOutputStepInterval(std::max(1, static_cast<int>(std::lround(outputTimeInterval / run.timeStep_))));
        simulation.setVTUOutputFormat(run.vtuFormat_);

        const int stepCount = run.stepCount_ >= 0 ? run.stepCount_ : static_cast<int>(std::ceil(run.duration_ / run.timeStep_));
        simulation.solve(stepCount);
    }
    catch (const std::exception& error)
    {
        std::cerr << "tutorial0 failed: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
