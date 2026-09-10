/** Checks level-set two-sphere central compression and surface-discretization convergence. */
#include "data/myLSObject.h"
#include "interaction/contactSearch.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "solver/LSDEM.h"
#include "validationCommon.h"

#include <algorithm>
#include <array>
#include <exception>

namespace
{

using namespace fundem;
using Real = math::Real;
using Vec3 = math::Vec3;

/** Summary of one surface-resolution sweep. */
struct resolutionResult {
    int surfaceNodeCount_{0};        ///< Number of surface quadrature nodes.
    Real maximumRelativeError_{0.0}; ///< Largest force error over all overlaps.
};

/** Evaluates all prescribed overlaps at one icosphere subdivision level. */
resolutionResult runResolution(int subdivisionLevel, std::ofstream& output)
{
    constexpr Real radius = 5.0e-3;
    constexpr Real normalStiffnessPerUnitArea = 1.0e12;
    constexpr Real density = 2500.0;
    constexpr std::array<Real, 3> overlaps{2.5e-5, 5.0e-5, 1.0e-4};

    LSDEM simulation;
    const int materialIndex = simulation.addMaterial(LSMaterial{normalStiffnessPerUnitArea, 0.0, 0.0, 1.0, density});
    levelset::Sphere geometry{radius};
    geometry.buildSurfaceNode(subdivisionLevel);
    geometry.buildLSGrid(60);
    const int geometryIndex = simulation.addGeometry(geometry);

    LSParticle movingSphere;
    movingSphere.setMaterial(simulation.materials(), materialIndex);
    movingSphere.setGeometry(simulation.geometries(), geometryIndex);
    simulation.addLSParticle(movingSphere);

    LSParticle fixedSphere;
    fixedSphere.setPosition({2.0 * radius, 0.0, 0.0});
    fixedSphere.setMaterial(simulation.materials(), materialIndex);
    fixedSphere.setGeometry(simulation.geometries(), geometryIndex);
    fixedSphere.setInfiniteMass();
    simulation.addLSParticle(fixedSphere);

    contactSearch search;
    contactContainer contacts;
    Real maximumRelativeError = 0.0;
    for (Real overlap : overlaps)
    {
        simulation.setLSParticlePosition(0, {overlap, 0.0, 0.0});
        search.findContacts(contacts, simulation.LSParticles());

        Vec3 resultantForce = Vec3::zero();
        for (contact& value : contacts.host())
        {
            validation::require(value.setMasterSlaveParticle(simulation.LSParticles(), value.masterParticleIndex(), value.slaveParticleIndex()), "Cannot refresh an LS contact.");
            validation::require(value.calculateForce(0.0), "Cannot evaluate an LS contact.");
            resultantForce += value.force();
        }

        const Real numericalForce = std::abs(resultantForce.x);
        const Real smallOverlapReference = 0.5 * math::pi * radius * normalStiffnessPerUnitArea * overlap * overlap;
        const Real relativeError = std::abs(numericalForce - smallOverlapReference) / smallOverlapReference;
        maximumRelativeError = std::max(maximumRelativeError, relativeError);
        output << subdivisionLevel << ' ' << simulation.geometries().hostSurfaceNodes().size() << ' ' << overlap << ' ' << numericalForce << ' ' << smallOverlapReference << ' ' << relativeError << ' '
               << contacts.hostSize() << '\n';
    }
    return {static_cast<int>(simulation.geometries().hostSurfaceNodes().size()), maximumRelativeError};
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        auto output = fundem::validation::openDataFile(fundem::validation::outputDirectory(argc, argv, "LSCentralCompression"), "LSCentralCompression.dat");
        output << "# subdivision surfaceNodes overlap numericalForce smallOverlapReference relativeError activeNodeContacts\n";

        const resolutionResult coarse = runResolution(2, output);
        const resolutionResult medium = runResolution(3, output);
        const resolutionResult refined = runResolution(4, output);
        const resolutionResult fine = runResolution(5, output);
        fundem::validation::require(coarse.surfaceNodeCount_ < medium.surfaceNodeCount_ && medium.surfaceNodeCount_ < refined.surfaceNodeCount_ && refined.surfaceNodeCount_ < fine.surfaceNodeCount_,
                                    "LS surface refinement did not increase node count.");
        fundem::validation::require(coarse.maximumRelativeError_ > medium.maximumRelativeError_ && medium.maximumRelativeError_ > refined.maximumRelativeError_ &&
                                        refined.maximumRelativeError_ > fine.maximumRelativeError_,
                                    "LS central-compression error did not decrease under surface refinement.");

        constexpr fundem::math::Real tolerance = 0.10;
        fundem::validation::require(fine.maximumRelativeError_ <= tolerance, "The finest LS central-compression response exceeded the small-overlap reference tolerance.");
        fundem::validation::printPass("LSCentralCompression", fine.maximumRelativeError_, tolerance);
    }
    catch (const std::exception& error)
    {
        std::cerr << "LSCentralCompression: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
