#include "data/HostAoSDeviceSoA.h"
#include "interaction/contact.h"
#include "interaction/surfaceNodeMapping.h"
#include "material/LSMaterial.h"
#include "particle/LSParticle.h"
#include "TestSupport.h"

#include <cstddef>
#include <utility>

namespace
{

struct baseRecord {
    float position{0.0F};

    struct positionField {
        inline static constexpr auto member = &baseRecord::position;
    };

    using DeviceLayout = fundem::deviceLayout<positionField>;
};

struct record : baseRecord {
    float velocity{0.0F};
    int identifier{0};
    double hostOnly{0.0};

    struct velocityField {
        inline static constexpr auto member = &record::velocity;
    };

    struct identifierField {
        inline static constexpr auto member = &record::identifier;
    };

    using DeviceLayout = fundem::concatDeviceLayoutsT<baseRecord::DeviceLayout, fundem::deviceLayout<velocityField, identifierField>>;
};

__global__ void advance(float* positions, float* velocities, int* identifiers, std::size_t count)
{
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < count)
    {
        positions[index] += velocities[index];
        identifiers[index] += 10;
    }
}

} // namespace

int main()
{
    using Storage = fundem::hostAoSDeviceSoA<record>;

    Storage storage;
    FUNDEM_TEST_REQUIRE(!storage.usesDevice());
    storage.setUseDevice(true);
    FUNDEM_TEST_REQUIRE(storage.usesDevice());
    storage.setUseDevice(false);
    FUNDEM_TEST_REQUIRE(!storage.usesDevice());
    storage.setUseDevice(true);
    storage.host().resize(3);

    for (std::size_t index = 0; index < storage.hostSize(); ++index)
    {
        storage.host()[index].position = static_cast<float>(index);
        storage.host()[index].velocity = 0.5F;
        storage.host()[index].identifier = static_cast<int>(index);
        storage.host()[index].hostOnly = 100.0 + static_cast<double>(index);
    }

    storage.copyHostToDevice();

    FUNDEM_TEST_REQUIRE(storage.deviceSize() == 3);
    FUNDEM_TEST_REQUIRE(storage.deviceCapacity() == 3);

    advance<<<1, 32>>>(storage.device<baseRecord::positionField>(), storage.device<record::velocityField>(), storage.device<record::identifierField>(), storage.deviceSize());

    FUNDEM_TEST_REQUIRE_CUDA(cudaGetLastError());
    storage.copyDeviceToHost();

    for (std::size_t index = 0; index < storage.hostSize(); ++index)
    {
        FUNDEM_TEST_REQUIRE(storage.host()[index].position == static_cast<float>(index) + 0.5F);
        FUNDEM_TEST_REQUIRE(storage.host()[index].identifier == static_cast<int>(index) + 10);
        FUNDEM_TEST_REQUIRE(storage.host()[index].hostOnly == 100.0 + static_cast<double>(index));
    }

    storage.host()[0].velocity = 4.0F;
    storage.copyHostFieldToDevice<record::velocityField>();
    storage.host()[0].velocity = -1.0F;
    storage.copyDeviceFieldToHost<record::velocityField>();
    FUNDEM_TEST_REQUIRE(storage.host()[0].velocity == 4.0F);

    storage.setDeviceSize(2);
    storage.copyDeviceToHost();
    FUNDEM_TEST_REQUIRE(storage.hostSize() == 2);

    Storage moved = std::move(storage);
    FUNDEM_TEST_REQUIRE(!storage.usesDevice());
    FUNDEM_TEST_REQUIRE(moved.usesDevice());
    FUNDEM_TEST_REQUIRE(storage.deviceSize() == 0);
    FUNDEM_TEST_REQUIRE(storage.deviceCapacity() == 0);
    FUNDEM_TEST_REQUIRE(moved.deviceSize() == 2);
    FUNDEM_TEST_REQUIRE(moved.deviceCapacity() == 3);

    moved.resetDevice();
    FUNDEM_TEST_REQUIRE(moved.deviceSize() == 0);
    FUNDEM_TEST_REQUIRE(moved.deviceCapacity() == 0);

    fundem::hostAoSDeviceSoA<fundem::LSParticle> particles;
    particles.host().emplace_back();
    particles.copyHostToDevice();
    particles.copyDeviceToHost();

    FUNDEM_TEST_REQUIRE(particles.hostSize() == 1);
    FUNDEM_TEST_REQUIRE(particles.deviceSize() == 1);

    fundem::contactContainer contacts;
    contacts.host().emplace_back();
    contacts.host()[0].setParticleSurfaceNodeIndex(17);
    contacts.copyHostToDevice();
    contacts.host()[0].setParticleSurfaceNodeIndex(-1);
    contacts.copyDeviceToHost();
    FUNDEM_TEST_REQUIRE(contacts.host()[0].particleSurfaceNodeIndex() == 17);

    fundem::surfaceNodeMappingContainer mappings;
    mappings.host().push_back({11, 3, 7});
    mappings.copyHostToDevice();
    mappings.host()[0].particleSurfaceNodeIndex_ = -1;
    mappings.copyDeviceToHost();
    FUNDEM_TEST_REQUIRE(mappings.host()[0].particleSurfaceNodeIndex_ == 7);

    fundem::materialContainer materials;
    materials.host().push_back(fundem::material{1.0e5, 5.0e4, 0.0, 0.0, 0.4, 0.0, 0.0, 0.5, 1000.0});
    materials.host().push_back(fundem::LSMaterial{1.0e5, 5.0e4, 0.4, 0.5, 1000.0});
    materials.copyHostToDevice();
    fundem::materialType deviceMaterialTypes[2]{};
    FUNDEM_TEST_REQUIRE_CUDA(cudaMemcpy(deviceMaterialTypes, materials.device<fundem::material::typeField>(), sizeof(deviceMaterialTypes), cudaMemcpyDeviceToHost));
    FUNDEM_TEST_REQUIRE(deviceMaterialTypes[0] == fundem::materialType::standard);
    FUNDEM_TEST_REQUIRE(deviceMaterialTypes[1] == fundem::materialType::levelSet);
}
