/**
 * @file tutorial2.cpp
 * @brief Demonstrates a bonded spherical-particle cloth falling onto a level-set box.
 */
#include "data/myLSObject.h"
#include "execution/cpu/parallelPolicy.h"
#include "interaction/bond.h"
#include "material/LSMaterial.h"
#include "material/material.h"
#include "particle/LSParticle.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"
#include "tutorialOptions.h"

#if FUNDEM_HAS_CUDA
#include "tutorial2DampingKernel.cuh"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

using fundem::math::Real;
using fundem::math::Vec3;

/** Applies the Workbench sample's linear and angular particle damping after force assembly. */
class ClothSimulation : public fundem::SphereDEM
{
public:
    /** Configures the execution backend and the two independent drag coefficients. */
    ClothSimulation(fundem::executionMode mode, int device, Real velocityDampingRate, Real rotationalDragCoefficient)
        : SphereDEM(mode, device), velocityDampingRate_(velocityDampingRate), rotationalDragCoefficient_(rotationalDragCoefficient)
    {}

protected:
    /** Adds F = -gamma m v and T = -c omega to finite-mass cloth spheres. */
    void addSphereExternalForceAndTorque(fundem::particleContainer::host_container_type& particles) override
    {
        const int particleCount = static_cast<int>(particles.size());
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (particleCount >= fundem::cpu::parallelParticleThreshold)
#endif
        for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
        {
            fundem::particle& value = particles[particleIndex];
            if (value.inverseMass() > 0.0)
            {
                value.addForce(-(velocityDampingRate_ / value.inverseMass()) * value.velocity());
                value.addTorque(-rotationalDragCoefficient_ * value.angularVelocity());
            }
        }
    }

#if FUNDEM_HAS_CUDA
    /** Enqueues the same damping for both GPU and hybrid execution. */
    void addSphereExternalForceAndTorque(fundem::particle::device_type particles, cudaStream_t stream) override
    {
        fundem::tutorial::launchClothDamping(particles, velocityDampingRate_, rotationalDragCoefficient_, stream);
    }
#endif

private:
    Real velocityDampingRate_;        ///< Mass-proportional translational damping rate, in inverse seconds.
    Real rotationalDragCoefficient_; ///< Isotropic rotational drag coefficient, in N m s.
};

/** Particle centers and neighbor edges of one triangular cloth lattice. */
struct ClothMesh {
    std::vector<Vec3> positions_;          ///< Row-major sphere centers matching the software HCP packing.
    std::vector<std::array<int, 2>> edges_; ///< Unique nearest-neighbor pairs with the smaller index first.
};

/** Builds a horizontal staggered triangular lattice and its nearest-neighbor edges. */
ClothMesh makeClothMesh(const Vec3& origin, int columnCount, int rowCount, Real particleDiameter)
{
    ClothMesh result;
    const Real radius = 0.5 * particleDiameter;
    const Real rowSpacing = std::sqrt(3.0) * radius;
    result.positions_.reserve(columnCount * rowCount);
    result.edges_.reserve(rowCount * (columnCount - 1) + (rowCount - 1) * (2 * columnCount - 1));

    for (int rowIndex = 0; rowIndex < rowCount; ++rowIndex)
    {
        const Real xOffset = rowIndex % 2 == 0 ? 0.0 : radius;
        for (int columnIndex = 0; columnIndex < columnCount; ++columnIndex)
        {
            const int particleIndex = rowIndex * columnCount + columnIndex;
            result.positions_.push_back(origin + Vec3{columnIndex * particleDiameter + xOffset, rowIndex * rowSpacing, 0.0});
            if (columnIndex > 0)
            {
                result.edges_.push_back({particleIndex - 1, particleIndex});
            }
            if (rowIndex > 0)
            {
                result.edges_.push_back({particleIndex - columnCount, particleIndex});
                const int secondNeighbor = rowIndex % 2 == 0 ? columnIndex - 1 : columnIndex + 1;
                if (secondNeighbor >= 0 && secondNeighbor < columnCount)
                {
                    result.edges_.push_back({(rowIndex - 1) * columnCount + secondNeighbor, particleIndex});
                }
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
        constexpr int clothColumnCount = 93;
        constexpr int clothRowCount = 107;
        constexpr Real particleDiameter = clothSideLength / (clothColumnCount + 0.5);
        constexpr Real particleRadius = 0.5 * particleDiameter;
        constexpr Real clothDensity = 1500.0;
        constexpr Real youngsModulus = 1.0e6;
        constexpr Real poissonRatio = 0.3;
        constexpr Real frictionCoefficient = 0.6;
        constexpr Real restitutionCoefficient = 0.05;
        constexpr Real velocityDampingRate = 20.0;
        constexpr Real angularDampingRate = 40.0;
        constexpr Real clothHeight = 0.38;
        constexpr Vec3 boxSize{0.24, 0.24, 0.24};
        constexpr Vec3 boxCenter{0.30, 0.30, 0.12};
        constexpr Real duration = 3.0;
        constexpr Real timeStep = 1.0e-5;
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
        constexpr Real rotationalDragCoefficient = angularDampingRate * (2.0 / 5.0) * particleMass * particleRadius * particleRadius;

        const tutorial::TutorialOptions run = tutorial::parseTutorialOptions(argc, argv, "tutorial2_files");
        ClothSimulation simulation(run.mode_, run.device_, velocityDampingRate, rotationalDragCoefficient);
        const int clothMaterialIndex = simulation.addMaterial(material{normalStiffness, shearStiffness, 0.0, 0.0, frictionCoefficient, 0.0, 0.0, restitutionCoefficient, clothDensity});

        const Vec3 clothOrigin{particleRadius, 0.5 * (clothSideLength - (clothRowCount - 1) * std::sqrt(3.0) * particleRadius), clothHeight};
        const ClothMesh cloth = makeClothMesh(clothOrigin, clothColumnCount, clothRowCount, particleDiameter);
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
            connection.setCrossSectionArea(0.0);
            connection.setModeMixityExponent(1.0);
            connection.setDamageInitiationRatio(1.0);
            if (!connection.setConnection(simulation.spheres(), edge[0], edge[1], math::normalizedOrZero(direction)) ||
                !connection.setStiffness(normalStiffness, shearStiffness, bendingStiffness, torsionalStiffness))
            {
                throw std::runtime_error("Failed to create the bonded cloth lattice.");
            }
            simulation.addBond(connection);
        }

        const int boxMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, frictionCoefficient, restitutionCoefficient, 1.0});
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
