/**
 * @file tutorial2.cpp
 * @brief Demonstrates a bonded spherical-particle cloth falling onto a level-set box.
 */
#include "data/myLSObject.h"
#include "interaction/bond.h"
#include "material/LSMaterial.h"
#include "material/material.h"
#include "particle/LSParticle.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"
#include "tutorialOptions.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace
{

using fundem::math::Real;
using fundem::math::Vec3;

/** Particle centers and neighbor edges of one triangular cloth lattice. */
struct ClothMesh {
    std::vector<Vec3> positions_;
    std::vector<std::array<int, 2>> edges_;
};

/** Builds a horizontal staggered triangular lattice and its nearest-neighbor edges. */
ClothMesh makeClothMesh(Real sideLength, Real particleDiameter, Real height)
{
    ClothMesh result;
    const Real radius = 0.5 * particleDiameter;
    const Real rowSpacing = std::sqrt(3.0) * radius;
    const Real maximumCenter = sideLength - radius;
    std::vector<std::vector<int>> rows;

    for (int rowIndex = 0;; ++rowIndex)
    {
        const Real y = radius + rowIndex * rowSpacing;
        if (y > maximumCenter + fundem::math::defaultTolerance)
        {
            break;
        }

        const Real xOffset = rowIndex % 2 == 0 ? 0.0 : radius;
        std::vector<int> row;
        for (Real x = radius + xOffset; x <= maximumCenter + fundem::math::defaultTolerance; x += particleDiameter)
        {
            row.push_back(static_cast<int>(result.positions_.size()));
            result.positions_.push_back({x, y, height});
        }
        if (!row.empty())
        {
            rows.push_back(std::move(row));
        }
    }

    for (int rowIndex = 0; rowIndex < static_cast<int>(rows.size()); ++rowIndex)
    {
        const std::vector<int>& row = rows[rowIndex];
        for (int columnIndex = 0; columnIndex + 1 < static_cast<int>(row.size()); ++columnIndex)
        {
            result.edges_.push_back({row[columnIndex], row[columnIndex + 1]});
        }
        if (rowIndex == 0)
        {
            continue;
        }

        const std::vector<int>& previousRow = rows[rowIndex - 1];
        for (int columnIndex = 0; columnIndex < static_cast<int>(row.size()); ++columnIndex)
        {
            if (columnIndex < static_cast<int>(previousRow.size()))
            {
                result.edges_.push_back({row[columnIndex], previousRow[columnIndex]});
            }
            const int secondNeighbor = rowIndex % 2 == 0 ? columnIndex - 1 : columnIndex + 1;
            if (secondNeighbor >= 0 && secondNeighbor < static_cast<int>(previousRow.size()))
            {
                result.edges_.push_back({row[columnIndex], previousRow[secondNeighbor]});
            }
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

        constexpr Real clothSideLength = 0.60;
        constexpr Real particleDiameter = 0.01;
        constexpr Real particleRadius = 0.5 * particleDiameter;
        constexpr Real clothDensity = 1500.0;
        constexpr Real youngsModulus = 1.0e6;
        constexpr Real poissonRatio = 0.3;
        constexpr Real frictionCoefficient = 0.4;
        constexpr Real restitutionCoefficient = 0.5;
        constexpr Real clothHeight = 0.38;
        constexpr Vec3 boxSize{0.24, 0.24, 0.24};
        constexpr Vec3 boxCenter{0.30, 0.30, 0.12};
        constexpr Real duration = 3.0;
        constexpr Real outputTimeInterval = 0.05;

        constexpr Real crossSectionArea = math::pi * particleRadius * particleRadius;
        constexpr Real areaMomentOfInertia = 0.25 * crossSectionArea * particleRadius * particleRadius;
        constexpr Real polarMomentOfInertia = 2.0 * areaMomentOfInertia;
        constexpr Real shearModulus = youngsModulus / (2.0 * (1.0 + poissonRatio));
        constexpr Real normalStiffness = youngsModulus * crossSectionArea / particleDiameter;
        constexpr Real shearStiffness = 12.0 * youngsModulus * areaMomentOfInertia / (particleDiameter * particleDiameter * particleDiameter);
        constexpr Real bendingStiffness = 0.01 * youngsModulus * areaMomentOfInertia / particleDiameter;
        constexpr Real torsionalStiffness = 0.01 * shearModulus * polarMomentOfInertia / particleDiameter;
        constexpr Real particleMass = (4.0 / 3.0) * math::pi * particleRadius * particleRadius * particleRadius * clothDensity;
        const Real timeStep = (math::pi / 50.0) * std::sqrt(particleMass / normalStiffness);

        const tutorial::TutorialOptions run = tutorial::parseTutorialOptions(argc, argv, "tutorial2_files");
        SphereDEM simulation(run.mode_, run.device_);
        const int clothMaterialIndex = simulation.addMaterial(material{normalStiffness, shearStiffness, 0.0, 0.0, frictionCoefficient, 0.0, 0.0, restitutionCoefficient, clothDensity});

        const ClothMesh cloth = makeClothMesh(clothSideLength, particleDiameter, clothHeight);
        for (const Vec3& position : cloth.positions_)
        {
            particle value;
            value.setPosition(position);
            value.setRadius(particleRadius);
            value.setMaterial(simulation.materials(), clothMaterialIndex);
            simulation.addSphere(value);
        }

        for (const std::array<int, 2>& edge : cloth.edges_)
        {
            const Vec3 direction = cloth.positions_[edge[0]] - cloth.positions_[edge[1]];
            const Real equivalentLength = math::norm(direction);
            bond connection{equivalentLength};
            if (!connection.setConnection(simulation.spheres(), edge[0], edge[1], math::normalizedOrZero(direction)) ||
                !connection.setStiffness(normalStiffness, shearStiffness, bendingStiffness, torsionalStiffness))
            {
                throw std::runtime_error("Failed to create the bonded cloth lattice.");
            }
            simulation.addBond(connection);
        }

        const int boxMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, frictionCoefficient, 1.0, 1.0});
        levelset::BoxWall box{boxSize};
        box.buildLSGrid(particleDiameter, 3);
        const int boxGeometryIndex = simulation.addGeometry(box, true);

        LSParticle fixedBox;
        fixedBox.setPosition(boxCenter);
        fixedBox.setMaterial(simulation.materials(), boxMaterialIndex);
        fixedBox.setGeometry(simulation.geometries(), boxGeometryIndex);
        simulation.addLSParticle(fixedBox);

        simulation.setBoundary({-0.20, -0.20, -0.60}, {0.80, 0.80, 0.70});
        simulation.setGravity({0.0, 0.0, -9.81});
        simulation.setTimeStep(timeStep);
        simulation.setOutputDirectory(run.outputDirectory_);
        simulation.setOutputStepInterval(std::max(1, static_cast<int>(std::lround(outputTimeInterval / timeStep))));
        simulation.setVTUOutputFormat(run.vtuFormat_);

        const int stepCount = run.stepCount_ >= 0 ? run.stepCount_ : static_cast<int>(std::ceil(duration / timeStep));
        simulation.solve(stepCount);
    }
    catch (const std::exception& error)
    {
        std::cerr << "tutorial2 failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
