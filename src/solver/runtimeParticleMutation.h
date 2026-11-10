#pragma once

#include "execution/motionIntegration.h"
#include "particle/rigidBody.h"

#include <stdexcept>
#include <vector>

namespace fundem::runtime_detail
{

inline std::vector<int> survivorMap(std::size_t count, const std::vector<int>& removed)
{
    std::vector<int> result(count, 0);
    for (const int index : removed)
    {
        if (index < 0 || index >= static_cast<int>(count) || result[index] < 0)
            throw std::invalid_argument("A runtime particle removal index is invalid or duplicated.");
        result[index] = -1;
    }
    int next = 0;
    for (int& index : result)
        if (index >= 0)
            index = next++;
    return result;
}

template <class Storage> void compact(Storage& storage, const std::vector<int>& map)
{
    typename Storage::host_container_type result;
    result.reserve(storage.hostSize());
    for (std::size_t index = 0; index < map.size(); ++index)
        if (map[index] >= 0)
            result.push_back(storage.host()[index]);
    storage.host() = std::move(result);
}

template <class Storage> math::Real resetVelocities(Storage& storage, const std::vector<int>& indices, cudaStream_t stream, bool resetLinear, bool resetAngular)
{
    math::Real removedEnergy = 0.0;
    for (const int index : indices)
    {
        auto& particle = storage.host()[index];
        removedEnergy += execution::kineticEnergy(resetLinear ? particle.velocity() : math::Vec3::zero(),
                                                 resetAngular ? particle.angularVelocity() : math::Vec3::zero(),
                                                 particle.orientation(), particle.inverseMass(), particle.inertiaTensor());
        if (resetLinear)
            particle.setVelocity(math::Vec3::zero());
        if (resetAngular)
            particle.setAngularVelocity(math::Vec3::zero());
    }
    if (storage.usesDevice())
    {
        if (resetLinear)
            storage.template copyHostFieldToDeviceAsync<pointMass::velocityField>(stream);
        if (resetAngular)
            storage.template copyHostFieldToDeviceAsync<rigidBody::angularVelocityField>(stream);
        host_device_detail::synchronize(stream);
    }
    return removedEnergy;
}

} // namespace fundem::runtime_detail
