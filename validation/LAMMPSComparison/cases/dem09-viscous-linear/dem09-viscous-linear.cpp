/** @file dem09-viscous-linear.cpp
 * @brief Direct native-core LAMMPS comparison: dem09-viscous-linear.
 * Reference: https://github.com/lammps/lammps/blob/29dd3fd7642a6fe9e8b4dd99120df4068b221250/unittest/granular/tests/dem09-viscous-linear-3d-si.yaml
 * Historical PASS is provenance; current acceptance uses only the declared checks.
 */
#include "comparisonCommon.h"

#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
/** One-particle Stokes drag, F = -rate * mass * velocity, in addition to gravity. */
class ViscousDragDEM final : public fundem::SphereDEM
{
public:
    using SphereDEM::SphereDEM;

protected:
    void addSphereExternalForceAndTorque(fundem::particleContainer::host_container_type& particles) override
    {
        for (fundem::particle& body : particles)
            if (body.inverseMass() > 0.0)
                body.addForce(-body.velocity() * (dragRate / body.inverseMass()));
    }

    void addSphereExternalForceAndTorque(fundem::particle::device_type particles, cudaStream_t stream) override
    {
#if FUNDEM_HAS_CUDA
        // Stage this one-particle fixture on the solver's stream; no alternate integrator.
        fundem::math::Vec3 force, velocity;
        fundem::math::Real inverseMass{};
        const auto check = [](cudaError_t value)
        {
            if (value != cudaSuccess)
                throw std::runtime_error(cudaGetErrorString(value));
        };
        check(cudaMemcpyAsync(&force, particles.force_, sizeof(force), cudaMemcpyDeviceToHost, stream));
        check(cudaMemcpyAsync(&velocity, particles.velocity_, sizeof(velocity), cudaMemcpyDeviceToHost, stream));
        check(cudaMemcpyAsync(&inverseMass, particles.inverseMass_, sizeof(inverseMass), cudaMemcpyDeviceToHost, stream));
        check(cudaStreamSynchronize(stream));
        if (inverseMass > 0.0)
            force -= velocity * (dragRate / inverseMass);
        check(cudaMemcpyAsync(particles.force_, &force, sizeof(force), cudaMemcpyHostToDevice, stream));
        check(cudaStreamSynchronize(stream));
#else
        (void)particles;
        (void)stream;
        throw std::runtime_error("This build has no CUDA backend.");
#endif
    }

private:
    static constexpr fundem::math::Real dragRate = 95.4929658551372; // 1/s; gamma / sphere mass.
};
} // namespace

int main(int argc, char** argv)
{
    try
    {
        using namespace fundem;
        const auto options = comparison::parseOptions(argc, argv, "dem09-viscous-linear");
        if (options.help)
        {
            comparison::printUsage(argv[0]);
            return 0;
        }

        ViscousDragDEM simulation(options.mode);
        const comparison::CaseSettings settings{
            "dem09-viscous-linear", "LAMMPS", "PASS",
            0.0001, 800, {0.0, 0.0, -9.8100000000000005},
            {
            }
        };
        std::vector<comparison::Body> bodies;

        // Bulk materials: Young's modulus, Poisson ratio, restitution, density.
        const int sphereMaterial = simulation.addMaterial(material{10000.0, 0.33000000000000002, 1.0, 2000.0});
        // Initial body states, in the original stable insertion order.
        particle sphereBody;
        sphereBody.setRadius(5.0000000000000002e-05);
        sphereBody.setMaterial(simulation.materials(), sphereMaterial);
        sphereBody.setPosition({0.0, 0.0, 0.0080000000000000002});
        sphereBody.setVelocity({0.0, 0.0, 0.0});
        sphereBody.setAngularVelocity({0.0, 0.0, 0.0});
        sphereBody.setOrientation(math::Quaternion{1.0, 0.0, 0.0, 0.0});
        const int sphereIndex = simulation.addSphere(sphereBody);
        bodies.push_back({"sphere", true, sphereIndex});

        // Pair switches apply to the complete four-direction stiffness/damping group.
        // Search bounds include the full original duration; these are not physical walls.
        simulation.setBoundary({-0.0001, -0.0001, -0.023492000000000006}, {0.0001, 0.0001, 0.039492000000000006});
        simulation.setGravity(settings.gravity);
        simulation.setTimeStep(settings.timeStep / (options.halfStep ? 2.0 : 1.0));
        return comparison::execute(simulation, settings, bodies, options);
    }
    catch (const std::exception& error)
    {
        std::cerr << "dem09-viscous-linear: " << error.what() << '\n';
        return 1;
    }
}
