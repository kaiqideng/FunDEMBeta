#include "SPHNeighborhood.h"
#include "execution/cpu/parallelPolicy.h"

#include <algorithm>
#include <limits>

namespace fundem
{

void SPHNeighborhood::captureHost(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries)
{
    auto& reference = positions_.host();
    const int fluidCount = static_cast<int>(particles.hostSize());
    const int count = fluidCount + static_cast<int>(boundaries.hostSize());
    reference.resize(count);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (count >= cpu::parallelParticleThreshold)
#endif
    for (int index = 0; index < count; ++index)
        reference[index].position_ = index < fluidCount ? particles.host()[index].position() : boundaries.host()[index - fluidCount].position();
    valid_ = true;
}

void SPHNeighborhood::captureDevice(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries, cudaStream_t stream)
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    const auto fluidCount = particles.deviceSize();
    const auto boundaryCount = boundaries.deviceSize();
    positions_.resizeDevice(fluidCount + boundaryCount);
    auto* reference = positions_.device<SPHNeighborPosition::positionField>();
    if (fluidCount > 0)
        host_device_detail::checkCuda(cudaMemcpyAsync(reference, particles.device<pointMass::positionField>(), fluidCount * sizeof(math::Vec3), cudaMemcpyDeviceToDevice, stream), "SPH neighbor reference copy");
    if (boundaryCount > 0)
        host_device_detail::checkCuda(cudaMemcpyAsync(reference + fluidCount, boundaries.device<pointMass::positionField>(), boundaryCount * sizeof(math::Vec3), cudaMemcpyDeviceToDevice, stream), "SPH boundary reference copy");
    valid_ = true;
#else
    (void)particles;
    (void)boundaries;
    (void)stream;
    throw std::runtime_error("Device neighborhood tracking requires CUDA.");
#endif
}

bool SPHNeighborhood::needsRebuild(const SPHParticleContainer& particles, const virtualParticleContainer& boundaries, math::Real skin) const
{
    const int fluidCount = static_cast<int>(particles.hostSize());
    const int count = fluidCount + static_cast<int>(boundaries.hostSize());
    if (!valid_ || positions_.hostSize() != static_cast<std::size_t>(count))
        return true;
    math::Real maximumSquaredDisplacement = 0.0;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (count >= cpu::parallelParticleThreshold) reduction(max : maximumSquaredDisplacement)
#endif
    for (int index = 0; index < count; ++index)
    {
        const math::Vec3 position = index < fluidCount ? particles.host()[index].position() : boundaries.host()[index - fluidCount].position();
        const math::Vec3 displacement = position - positions_.host()[index].position_;
        const math::Real squared = math::dot(displacement, displacement);
        maximumSquaredDisplacement = std::max(maximumSquaredDisplacement, math::isFinite(squared) ? squared : std::numeric_limits<math::Real>::infinity());
    }
    if (!math::isFinite(maximumSquaredDisplacement))
        throw std::runtime_error("Non-finite SPH position while checking the neighborhood displacement.");
    return exceedsSkin(maximumSquaredDisplacement, skin);
}

} // namespace fundem
