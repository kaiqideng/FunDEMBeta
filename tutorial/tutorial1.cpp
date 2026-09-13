/** CPU level-set anchor-chain tutorial derived from the CollisionLab chain-link case. */
#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "solver/LSDEM.h"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <vector>

namespace
{

using fundem::math::Real;
using fundem::math::Vec3;

struct linkRing {
    Vec3 center_;
    Vec3 normal_;
};

fundem::levelset::TriangleMesh makeChainLink()
{
    constexpr Real straightHalfLength = 0.03;
    constexpr Real endRadius = 0.025;
    constexpr Real rodRadius = 0.01;
    constexpr int straightSegmentCount = 30;
    constexpr int arcSegmentCount = 40;
    constexpr int crossSectionSegmentCount = 32;

    std::vector<linkRing> rings;
    rings.reserve(2 * (straightSegmentCount + arcSegmentCount));
    for (int index = 0; index < straightSegmentCount; ++index)
    {
        const Real fraction = Real(index) / straightSegmentCount;
        rings.push_back({{-straightHalfLength + 2.0 * straightHalfLength * fraction, endRadius, 0.0}, Vec3::unitY()});
    }
    for (int index = 0; index < arcSegmentCount; ++index)
    {
        const Real angle = 0.5 * fundem::math::pi - fundem::math::pi * Real(index) / arcSegmentCount;
        rings.push_back({{straightHalfLength + endRadius * std::cos(angle), endRadius * std::sin(angle), 0.0}, {std::cos(angle), std::sin(angle), 0.0}});
    }
    for (int index = 0; index < straightSegmentCount; ++index)
    {
        const Real fraction = Real(index) / straightSegmentCount;
        rings.push_back({{straightHalfLength - 2.0 * straightHalfLength * fraction, -endRadius, 0.0}, -Vec3::unitY()});
    }
    for (int index = 0; index < arcSegmentCount; ++index)
    {
        const Real angle = -0.5 * fundem::math::pi - fundem::math::pi * Real(index) / arcSegmentCount;
        rings.push_back({{-straightHalfLength + endRadius * std::cos(angle), endRadius * std::sin(angle), 0.0}, {std::cos(angle), std::sin(angle), 0.0}});
    }

    const int ringCount = static_cast<int>(rings.size());
    std::vector<Vec3> vertices;
    vertices.reserve(ringCount * crossSectionSegmentCount);
    for (const linkRing& ring : rings)
    {
        for (int index = 0; index < crossSectionSegmentCount; ++index)
        {
            const Real angle = 2.0 * fundem::math::pi * Real(index) / crossSectionSegmentCount;
            vertices.push_back(ring.center_ + rodRadius * std::cos(angle) * ring.normal_ + Vec3{0.0, 0.0, rodRadius * std::sin(angle)});
        }
    }

    std::vector<int3> triangles;
    triangles.reserve(2 * ringCount * crossSectionSegmentCount);
    for (int ringIndex = 0; ringIndex < ringCount; ++ringIndex)
    {
        const int nextRingIndex = (ringIndex + 1) % ringCount;
        for (int sectionIndex = 0; sectionIndex < crossSectionSegmentCount; ++sectionIndex)
        {
            const int nextSectionIndex = (sectionIndex + 1) % crossSectionSegmentCount;
            const int first = ringIndex * crossSectionSegmentCount + sectionIndex;
            const int second = ringIndex * crossSectionSegmentCount + nextSectionIndex;
            const int third = nextRingIndex * crossSectionSegmentCount + nextSectionIndex;
            const int fourth = nextRingIndex * crossSectionSegmentCount + sectionIndex;
            triangles.push_back({first, second, third});
            triangles.push_back({first, third, fourth});
        }
    }
    return {vertices, triangles};
}

std::filesystem::path executableDirectory(const char* executable)
{
    std::error_code error;
    const std::filesystem::path absolutePath = std::filesystem::absolute(executable, error);
    return error ? std::filesystem::current_path() : absolutePath.parent_path();
}

} // namespace

int main(int, char** argv)
{
    try
    {
        using namespace fundem;

        constexpr int particleCount = 11;
        constexpr Real particleSpacing = 0.07;
        constexpr Real initialHeight = 0.3;
        constexpr Real gridSpacing = 0.002;
        constexpr Real timeStep = 2.0e-5;
        constexpr Real duration = 3.0;
        constexpr Real outputTimeInterval = 0.05;

        LSDEM simulation;
        const int materialIndex = simulation.addMaterial(LSMaterial{2.0e10, 8.0e9, 0.4, 0.2, 7850.0});

        levelset::TriangleMesh chainLink = makeChainLink();
        chainLink.buildLSGrid(gridSpacing, 2);
        const int geometryIndex = simulation.addGeometry(chainLink);
        const math::Quaternion perpendicularOrientation = math::Quaternion::fromUnitAxisAngle(Vec3::unitX(), 0.5 * math::pi);

        for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
        {
            LSParticle value;
            value.setPosition({(particleIndex - particleCount / 2) * particleSpacing, 0.0, initialHeight});
            value.setOrientation(particleIndex % 2 == 0 ? math::Quaternion::identity() : perpendicularOrientation);
            value.setMaterial(simulation.materials(), materialIndex);
            value.setGeometry(simulation.geometries(), geometryIndex);
            if (particleIndex == 0 || particleIndex == particleCount - 1)
            {
                value.setInfiniteMass();
            }
            simulation.addLSParticle(value);
        }

        simulation.setBoundary({-0.5, -0.2, -0.2}, {0.5, 0.2, 0.5});
        simulation.setGravity({0.0, 0.0, -9.81});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory((executableDirectory(argv[0]) / "tutorial1_files").string());
        simulation.setOutputStepInterval(static_cast<int>(std::lround(outputTimeInterval / timeStep)));
        simulation.solve(static_cast<int>(std::ceil(duration / timeStep)));
    }
    catch (const std::exception& error)
    {
        std::cerr << "tutorial1 failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
