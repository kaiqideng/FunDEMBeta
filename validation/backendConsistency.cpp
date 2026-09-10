/** Compares an identical sphere-level-set contact trajectory across supported backends. */
#include "data/CudaTypes.h"
#include "data/myLSObject.h"
#include "material/LSMaterial.h"
#include "material/material.h"
#include "particle/LSParticle.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"
#include "validationCommon.h"

#include <algorithm>
#include <exception>

namespace
{

using namespace fundem;
using Real = math::Real;
using Vec3 = math::Vec3;

/** Final sphere kinematics retained for backend comparison. */
struct state {
    Vec3 position_{Vec3::zero()};        ///< Final center position.
    Vec3 velocity_{Vec3::zero()};        ///< Final linear velocity.
    Vec3 angularVelocity_{Vec3::zero()}; ///< Final angular velocity.
};

/** Runs the same short sphere-plane trajectory in one execution mode. */
state run(executionMode mode)
{
    SphereDEM simulation;
    simulation.setExecutionMode(mode);
    const int sphereMaterialIndex = simulation.addMaterial(material{2.0e5, 8.0e4, 0.0, 0.0, 0.25, 0.0, 0.0, 0.7, 1000.0});
    const int wallMaterialIndex = simulation.addMaterial(LSMaterial{0.0, 0.0, 0.25, 0.7, 1.0});

    particle sphere;
    sphere.setPosition({0.0, 0.0, 0.0495});
    sphere.setVelocity({0.1, 0.0, -0.05});
    sphere.setRadius(0.05);
    sphere.setMaterial(simulation.materials(), sphereMaterialIndex);
    simulation.addSphere(sphere);

    levelset::PlaneWall plane{Vec3::unitZ(), 1.0};
    plane.buildLSGrid(0.05, 2);
    const int geometryIndex = simulation.addGeometry(plane, true);
    LSParticle wall;
    wall.setMaterial(simulation.materials(), wallMaterialIndex);
    wall.setGeometry(simulation.geometries(), geometryIndex);
    simulation.addLSParticle(wall);

    simulation.setBoundary({-0.5, -0.5, -0.2}, {0.5, 0.5, 0.5});
    simulation.setGravity({0.0, 0.0, -9.81});
    simulation.setTimeStep(1.0e-5);
    simulation.solve(200);

    const particle& result = simulation.spheres().host().front();
    return {result.position(), result.velocity(), result.angularVelocity()};
}

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
/** Returns the largest scale-aware kinematic difference from a reference. */
Real stateError(const state& value, const state& reference) noexcept
{
    return std::max({validation::relativeError(value.position_, reference.position_),
                     validation::relativeError(value.velocity_, reference.velocity_),
                     validation::relativeError(value.angularVelocity_, reference.angularVelocity_)});
}
#endif

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const state CPUState = run(fundem::executionMode::CPU);
        auto output = fundem::validation::openDataFile(fundem::validation::outputDirectory(argc, argv, "backendConsistency"), "backendConsistency.dat");
        output << "# backend x y z vx vy vz omegaX omegaY omegaZ relativeErrorFromCPU\n";
        output << "CPU " << CPUState.position_.x << ' ' << CPUState.position_.y << ' ' << CPUState.position_.z << ' ' << CPUState.velocity_.x << ' ' << CPUState.velocity_.y << ' '
               << CPUState.velocity_.z << ' ' << CPUState.angularVelocity_.x << ' ' << CPUState.angularVelocity_.y << ' ' << CPUState.angularVelocity_.z << " 0\n";

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
        int deviceCount = 0;
        if (cudaGetDeviceCount(&deviceCount) != cudaSuccess || deviceCount <= 0)
        {
            std::cout << "backendConsistency: SKIP GPU and Hybrid (no CUDA device available)\n";
            return 0;
        }

        const state GPUState = run(fundem::executionMode::GPU);
        const state hybridState = run(fundem::executionMode::Hybrid);
        const Real GPUError = stateError(GPUState, CPUState);
        const Real hybridError = stateError(hybridState, CPUState);
        const Real maximumError = std::max(GPUError, hybridError);
        constexpr Real tolerance = 1.0e-9;
        output << "GPU " << GPUState.position_.x << ' ' << GPUState.position_.y << ' ' << GPUState.position_.z << ' ' << GPUState.velocity_.x << ' ' << GPUState.velocity_.y << ' '
               << GPUState.velocity_.z << ' ' << GPUState.angularVelocity_.x << ' ' << GPUState.angularVelocity_.y << ' ' << GPUState.angularVelocity_.z << ' ' << GPUError << '\n';
        output << "Hybrid " << hybridState.position_.x << ' ' << hybridState.position_.y << ' ' << hybridState.position_.z << ' ' << hybridState.velocity_.x << ' ' << hybridState.velocity_.y << ' '
               << hybridState.velocity_.z << ' ' << hybridState.angularVelocity_.x << ' ' << hybridState.angularVelocity_.y << ' ' << hybridState.angularVelocity_.z << ' ' << hybridError << '\n';
        fundem::validation::require(maximumError <= tolerance, "CPU, GPU, and Hybrid trajectories are inconsistent.");
        fundem::validation::printPass("backendConsistency", maximumError, tolerance);
#else
        std::cout << "backendConsistency: SKIP GPU and Hybrid (CUDA backend not built)\n";
#endif
    }
    catch (const std::exception& error)
    {
        std::cerr << "backendConsistency: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
