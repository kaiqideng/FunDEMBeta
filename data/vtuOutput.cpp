#include "vtuOutput.h"

#include "execution/motionIntegration.h"

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

namespace fundem
{
namespace
{

using Real = math::Real;
using Vec3 = math::Vec3;

constexpr int parallelParticleThreshold = 1024;
constexpr int parallelSurfacePointThreshold = 4096;

template <class Storage, class Getter> auto collectValues(const Storage& storage, Getter getter)
{
    using Value = std::decay_t<decltype(getter(std::declval<const typename Storage::host_type&>()))>;
    const auto& host = storage.host();
    const int valueCount = static_cast<int>(host.size());
    std::vector<Value> values(valueCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (valueCount >= parallelParticleThreshold)
#endif
    for (int index = 0; index < valueCount; ++index)
    {
        values[index] = getter(host[index]);
    }
    return values;
}

template <class MasterStorage, class SlaveStorage, class Getter>
auto collectBondParticleValues(const bondContainer& bonds, const MasterStorage& masterParticles, const SlaveStorage& slaveParticles, Getter getter)
{
    using Value = std::decay_t<decltype(getter(std::declval<const bond&>(), std::declval<const typename MasterStorage::host_type&>(), std::declval<const typename SlaveStorage::host_type&>()))>;
    const auto& bondHost = bonds.host();
    const int bondCount = static_cast<int>(bondHost.size());
    std::vector<Value> values(bondCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (bondCount >= parallelParticleThreshold)
#endif
    for (int bondIndex = 0; bondIndex < bondCount; ++bondIndex)
    {
        const bond& value = bondHost[bondIndex];
        values[bondIndex] = getter(value, masterParticles.host()[value.masterParticleIndex()], slaveParticles.host()[value.slaveParticleIndex()]);
    }
    return values;
}

template <class Storage, class Getter> std::vector<Real> collectQuaternions(const Storage& storage, Getter getter)
{
    const auto& host = storage.host();
    const int valueCount = static_cast<int>(host.size());
    std::vector<Real> values(4 * valueCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (valueCount >= parallelParticleThreshold)
#endif
    for (int index = 0; index < valueCount; ++index)
    {
        const math::Quaternion& orientation = getter(host[index]);
        const int offset = 4 * index;
        values[offset] = orientation.w;
        values[offset + 1] = orientation.x;
        values[offset + 2] = orientation.y;
        values[offset + 3] = orientation.z;
    }
    return values;
}

template <class Storage, class Getter> std::vector<Real> collectMatrices(const Storage& storage, Getter getter)
{
    const auto& host = storage.host();
    const int valueCount = static_cast<int>(host.size());
    std::vector<Real> values(math::Mat3::size() * valueCount);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (valueCount >= parallelParticleThreshold)
#endif
    for (int index = 0; index < valueCount; ++index)
    {
        const math::Mat3& matrix = getter(host[index]);
        std::copy(matrix.data(), matrix.data() + math::Mat3::size(), values.begin() + math::Mat3::size() * index);
    }
    return values;
}

bool selectLSParticle(const LSParticle& particle, bool fixed) noexcept { return (particle.inverseMass() == 0.0) == fixed; }

std::vector<int> surfacePointOffsets(const LSParticleContainer& particles, bool fixed)
{
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    std::vector<int> offsets(particleCount + 1, 0);
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const LSParticle& particle = particleHost[particleIndex];
        if (selectLSParticle(particle, fixed))
        {
            offsets[particleIndex + 1] = static_cast<int>(particle.surfaceNodes().size());
        }
    }
    std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
    return offsets;
}

template <class Getter> auto collectSurfacePointValues(const LSParticleContainer& particles, bool fixed, Getter getter)
{
    using Value = std::decay_t<decltype(getter(std::declval<const LSParticle&>(), 0, 0))>;
    const auto offsets = surfacePointOffsets(particles, fixed);
    const int particleCount = static_cast<int>(particles.hostSize());
    std::vector<Value> values(offsets.back());
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 4) if (offsets.back() >= parallelSurfacePointThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const LSParticle& particle = particles.host()[particleIndex];
        if (!selectLSParticle(particle, fixed))
        {
            continue;
        }
        for (int surfaceNodeIndex = 0; surfaceNodeIndex < static_cast<int>(particle.surfaceNodes().size()); ++surfaceNodeIndex)
        {
            values[offsets[particleIndex] + surfaceNodeIndex] = getter(particle, particleIndex, surfaceNodeIndex);
        }
    }
    return values;
}

template <class Getter> std::vector<Real> collectSurfacePointQuaternions(const LSParticleContainer& particles, bool fixed, Getter getter)
{
    const auto offsets = surfacePointOffsets(particles, fixed);
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    std::vector<Real> values(4 * offsets.back());
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 4) if (offsets.back() >= parallelSurfacePointThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const LSParticle& particle = particleHost[particleIndex];
        if (!selectLSParticle(particle, fixed))
        {
            continue;
        }
        const math::Quaternion& orientation = getter(particle);
        for (int surfaceNodeIndex = 0; surfaceNodeIndex < static_cast<int>(particle.surfaceNodes().size()); ++surfaceNodeIndex)
        {
            const int offset = 4 * (offsets[particleIndex] + surfaceNodeIndex);
            values[offset] = orientation.w;
            values[offset + 1] = orientation.x;
            values[offset + 2] = orientation.y;
            values[offset + 3] = orientation.z;
        }
    }
    return values;
}

template <class Getter> std::vector<Real> collectSurfacePointMatrices(const LSParticleContainer& particles, bool fixed, Getter getter)
{
    const auto offsets = surfacePointOffsets(particles, fixed);
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    std::vector<Real> values(math::Mat3::size() * offsets.back());
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 4) if (offsets.back() >= parallelSurfacePointThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const LSParticle& particle = particleHost[particleIndex];
        if (!selectLSParticle(particle, fixed))
        {
            continue;
        }
        const math::Mat3& matrix = getter(particle);
        for (int surfaceNodeIndex = 0; surfaceNodeIndex < static_cast<int>(particle.surfaceNodes().size()); ++surfaceNodeIndex)
        {
            std::copy(matrix.data(), matrix.data() + math::Mat3::size(), values.begin() + math::Mat3::size() * (offsets[particleIndex] + surfaceNodeIndex));
        }
    }
    return values;
}

void addSphereField(vtuWriter& writer, const particleContainer& spheres, sphereVTUField field, const Vec3& gravity)
{
    switch (field)
    {
    case sphereVTUField::orientation:
        writer.addPointData("orientation", collectQuaternions(spheres, [](const particle& value) -> const math::Quaternion& { return value.orientation(); }), 4);
        break;
    case sphereVTUField::angularVelocity:
        writer.addPointData("angularVelocity", collectValues(spheres, [](const particle& value) { return value.angularVelocity(); }));
        break;
    case sphereVTUField::force:
        writer.addPointData("force", collectValues(spheres, [](const particle& value) { return value.force(); }));
        break;
    case sphereVTUField::torque:
        writer.addPointData("torque", collectValues(spheres, [](const particle& value) { return value.torque(); }));
        break;
    case sphereVTUField::inverseMass:
        writer.addPointData("inverseMass", collectValues(spheres, [](const particle& value) { return value.inverseMass(); }));
        break;
    case sphereVTUField::inertiaTensor:
        writer.addPointData("inertiaTensor", collectMatrices(spheres, [](const particle& value) -> const math::Mat3& { return value.inertiaTensor(); }), 9);
        break;
    case sphereVTUField::materialIndex:
        writer.addPointData("materialIndex", collectValues(spheres, [](const particle& value) { return value.materialIndex(); }));
        break;
    case sphereVTUField::kineticEnergy:
        writer.addPointData("kineticEnergy",
                            collectValues(spheres,
                                          [](const particle& value)
                                          { return execution::kineticEnergy(value.velocity(), value.angularVelocity(), value.orientation(), value.inverseMass(), value.inertiaTensor()); }));
        break;
    case sphereVTUField::gravitationalPotentialEnergy:
        writer.addPointData("gravitationalPotentialEnergy",
                            collectValues(spheres, [&gravity](const particle& value) { return execution::gravitationalPotentialEnergy(value.position(), value.inverseMass(), gravity); }));
        break;
    }
}

void addLSParticleField(vtuWriter& writer, const LSParticleContainer& particles, bool fixed, LSParticleVTUField field, const Vec3& gravity)
{
    switch (field)
    {
    case LSParticleVTUField::particleIndex:
        writer.addPointData("particleIndex", collectSurfacePointValues(particles, fixed, [](const LSParticle&, int particleIndex, int) { return particleIndex; }));
        break;
    case LSParticleVTUField::distanceToCenter:
        writer.addPointData(
            "distanceToCenter",
            collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int surfaceNodeIndex) { return math::norm(value.surfaceNodes()[surfaceNodeIndex].position_); }));
        break;
    case LSParticleVTUField::orientation:
        writer.addPointData("orientation", collectSurfacePointQuaternions(particles, fixed, [](const LSParticle& value) -> const math::Quaternion& { return value.orientation(); }), 4);
        break;
    case LSParticleVTUField::angularVelocity:
        writer.addPointData("angularVelocity", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.angularVelocity(); }));
        break;
    case LSParticleVTUField::force:
        writer.addPointData("force", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.force(); }));
        break;
    case LSParticleVTUField::torque:
        writer.addPointData("torque", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.torque(); }));
        break;
    case LSParticleVTUField::inertiaTensor:
        writer.addPointData("inertiaTensor", collectSurfacePointMatrices(particles, fixed, [](const LSParticle& value) -> const math::Mat3& { return value.inertiaTensor(); }), 9);
        break;
    case LSParticleVTUField::boundingRadius:
        writer.addPointData("boundingRadius", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.boundingRadius(); }));
        break;
    case LSParticleVTUField::materialIndex:
        writer.addPointData("materialIndex", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.materialIndex(); }));
        break;
    case LSParticleVTUField::volume:
        writer.addPointData("volume", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.volume(); }));
        break;
    case LSParticleVTUField::geometryIndex:
        writer.addPointData("geometryIndex", collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int) { return value.geometryIndex(); }));
        break;
    case LSParticleVTUField::surfaceNodeArea:
        writer.addPointData("surfaceNodeArea",
                            collectSurfacePointValues(particles, fixed, [](const LSParticle& value, int, int surfaceNodeIndex) { return value.surfaceNodes()[surfaceNodeIndex].area_; }));
        break;
    case LSParticleVTUField::kineticEnergy:
        writer.addPointData("kineticEnergy",
                            collectSurfacePointValues(particles,
                                                      fixed,
                                                      [](const LSParticle& value, int, int) {
                                                          return execution::kineticEnergy(value.velocity(), value.angularVelocity(), value.orientation(), value.inverseMass(), value.inertiaTensor());
                                                      }));
        break;
    case LSParticleVTUField::gravitationalPotentialEnergy:
        writer.addPointData("gravitationalPotentialEnergy",
                            collectSurfacePointValues(particles,
                                                      fixed,
                                                      [&gravity](const LSParticle& value, int, int)
                                                      { return execution::gravitationalPotentialEnergy(value.position(), value.inverseMass(), gravity); }));
        break;
    }
}

void addSPHParticleField(vtuWriter& writer, const SPHParticleContainer& particles, SPHParticleVTUField field, const Vec3& gravity)
{
    switch (field)
    {
    case SPHParticleVTUField::force:
        writer.addPointData("force", collectValues(particles, [](const SPHParticle& value) { return value.force(); }));
        break;
    case SPHParticleVTUField::pressure:
        writer.addPointData("pressure", collectValues(particles, [](const SPHParticle& value) { return value.pressure(); }));
        break;
    case SPHParticleVTUField::densityRate:
        writer.addPointData("densityRate", collectValues(particles, [](const SPHParticle& value) { return value.densityRate(); }));
        break;
    case SPHParticleVTUField::freeSurface:
        writer.addPointData("freeSurface", collectValues(particles, [](const SPHParticle& value) { return value.isFreeSurface() ? 1 : 0; }));
        break;
    case SPHParticleVTUField::kineticEnergy:
        writer.addPointData(
            "kineticEnergy",
            collectValues(particles, [](const SPHParticle& value) { return value.inverseMass() > 0.0 ? 0.5 * math::dot(value.velocity(), value.velocity()) / value.inverseMass() : 0.0; }));
        break;
    case SPHParticleVTUField::gravitationalPotentialEnergy:
        writer.addPointData("gravitationalPotentialEnergy",
                            collectValues(particles, [&gravity](const SPHParticle& value) { return execution::gravitationalPotentialEnergy(value.position(), value.inverseMass(), gravity); }));
        break;
    }
}

void addContactField(vtuWriter& writer, const contactContainer& contacts, contactVTUField field)
{
    switch (field)
    {
    case contactVTUField::force:
        writer.addPointData("force", collectValues(contacts, [](const contact& value) { return value.force(); }));
        break;
    case contactVTUField::torque:
        writer.addPointData("torque", collectValues(contacts, [](const contact& value) { return value.torque(); }));
        break;
    case contactVTUField::masterParticleIndex:
        writer.addPointData("masterParticleIndex", collectValues(contacts, [](const contact& value) { return value.masterParticleIndex(); }));
        break;
    case contactVTUField::slaveParticleIndex:
        writer.addPointData("slaveParticleIndex", collectValues(contacts, [](const contact& value) { return value.slaveParticleIndex(); }));
        break;
    case contactVTUField::overlap:
        writer.addPointData("overlap", collectValues(contacts, [](const contact& value) { return value.overlap(); }));
        break;
    case contactVTUField::area:
        writer.addPointData("area", collectValues(contacts, [](const contact& value) { return value.area(); }));
        break;
    case contactVTUField::effectiveMass:
        writer.addPointData("effectiveMass", collectValues(contacts, [](const contact& value) { return value.effectiveMass(); }));
        break;
    case contactVTUField::effectiveRadius:
        writer.addPointData("effectiveRadius", collectValues(contacts, [](const contact& value) { return value.effectiveRadius(); }));
        break;
    case contactVTUField::normalForceMagnitude:
        writer.addPointData("normalForceMagnitude", collectValues(contacts, [](const contact& value) { return value.normalForceMagnitude(); }));
        break;
    case contactVTUField::normalElasticEnergy:
        writer.addPointData("normalElasticEnergy", collectValues(contacts, [](const contact& value) { return value.normalElasticEnergy(); }));
        break;
    case contactVTUField::slidingElasticEnergy:
        writer.addPointData("slidingElasticEnergy", collectValues(contacts, [](const contact& value) { return value.slidingElasticEnergy(); }));
        break;
    case contactVTUField::rollingElasticEnergy:
        writer.addPointData("rollingElasticEnergy", collectValues(contacts, [](const contact& value) { return value.rollingElasticEnergy(); }));
        break;
    case contactVTUField::torsionalElasticEnergy:
        writer.addPointData("torsionalElasticEnergy", collectValues(contacts, [](const contact& value) { return value.torsionalElasticEnergy(); }));
        break;
    case contactVTUField::slidingSpringDeformation:
        writer.addPointData("slidingSpringDeformation", collectValues(contacts, [](const contact& value) { return value.slidingSpringDeformation(); }));
        break;
    case contactVTUField::rollingSpringDeformation:
        writer.addPointData("rollingSpringDeformation", collectValues(contacts, [](const contact& value) { return value.rollingSpringDeformation(); }));
        break;
    case contactVTUField::torsionalSpringDeformation:
        writer.addPointData("torsionalSpringDeformation", collectValues(contacts, [](const contact& value) { return value.torsionalSpringDeformation(); }));
        break;
    }
}

template <class MasterStorage, class SlaveStorage>
void addBondField(vtuWriter& writer, const bondContainer& bonds, const MasterStorage& masterParticles, const SlaveStorage& slaveParticles, bondVTUField field)
{
    switch (field)
    {
    case bondVTUField::force:
        writer.addPointData("force", collectValues(bonds, [](const bond& value) { return value.force(); }));
        break;
    case bondVTUField::masterTorque:
        writer.addPointData("masterTorque", collectValues(bonds, [](const bond& value) { return value.masterTorque(); }));
        break;
    case bondVTUField::slaveTorque:
        writer.addPointData("slaveTorque", collectValues(bonds, [](const bond& value) { return value.slaveTorque(); }));
        break;
    case bondVTUField::masterParticleIndex:
        writer.addPointData("masterParticleIndex", collectValues(bonds, [](const bond& value) { return value.masterParticleIndex(); }));
        break;
    case bondVTUField::slaveParticleIndex:
        writer.addPointData("slaveParticleIndex", collectValues(bonds, [](const bond& value) { return value.slaveParticleIndex(); }));
        break;
    case bondVTUField::equivalentLength:
        writer.addPointData("equivalentLength", collectValues(bonds, [](const bond& value) { return value.equivalentLength(); }));
        break;
    case bondVTUField::normalElasticEnergy:
        writer.addPointData("normalElasticEnergy", collectValues(bonds, [](const bond& value) { return value.normalElasticEnergy(); }));
        break;
    case bondVTUField::shearElasticEnergy:
        writer.addPointData("shearElasticEnergy", collectValues(bonds, [](const bond& value) { return value.shearElasticEnergy(); }));
        break;
    case bondVTUField::bendingElasticEnergy:
        writer.addPointData("bendingElasticEnergy", collectValues(bonds, [](const bond& value) { return value.bendingElasticEnergy(); }));
        break;
    case bondVTUField::torsionalElasticEnergy:
        writer.addPointData("torsionalElasticEnergy", collectValues(bonds, [](const bond& value) { return value.torsionalElasticEnergy(); }));
        break;
    case bondVTUField::masterEndpointNormal:
        writer.addPointData("masterEndpointNormal",
                            collectBondParticleValues(bonds,
                                                      masterParticles,
                                                      slaveParticles,
                                                      [](const bond& value, const auto& master, const auto&) { return math::rotateUnit(master.orientation(), value.masterEndpointLocalNormal()); }));
        break;
    case bondVTUField::masterEndpointTangent1:
        writer.addPointData("masterEndpointTangent1",
                            collectBondParticleValues(bonds,
                                                      masterParticles,
                                                      slaveParticles,
                                                      [](const bond& value, const auto& master, const auto&) { return math::rotateUnit(master.orientation(), value.masterEndpointLocalTangent1()); }));
        break;
    case bondVTUField::masterEndpointTangent2:
        writer.addPointData("masterEndpointTangent2",
                            collectBondParticleValues(bonds,
                                                      masterParticles,
                                                      slaveParticles,
                                                      [](const bond& value, const auto& master, const auto&) { return math::rotateUnit(master.orientation(), value.masterEndpointLocalTangent2()); }));
        break;
    case bondVTUField::slaveEndpointNormal:
        writer.addPointData("slaveEndpointNormal",
                            collectBondParticleValues(bonds,
                                                      masterParticles,
                                                      slaveParticles,
                                                      [](const bond& value, const auto&, const auto& slave) { return math::rotateUnit(slave.orientation(), value.slaveEndpointLocalNormal()); }));
        break;
    case bondVTUField::slaveEndpointTangent1:
        writer.addPointData("slaveEndpointTangent1",
                            collectBondParticleValues(bonds,
                                                      masterParticles,
                                                      slaveParticles,
                                                      [](const bond& value, const auto&, const auto& slave) { return math::rotateUnit(slave.orientation(), value.slaveEndpointLocalTangent1()); }));
        break;
    case bondVTUField::slaveEndpointTangent2:
        writer.addPointData("slaveEndpointTangent2",
                            collectBondParticleValues(bonds,
                                                      masterParticles,
                                                      slaveParticles,
                                                      [](const bond& value, const auto&, const auto& slave) { return math::rotateUnit(slave.orientation(), value.slaveEndpointLocalTangent2()); }));
        break;
    case bondVTUField::fractureArea:
        writer.addPointData("fractureArea", collectValues(bonds, [](const bond& value) { return value.fractureArea(); }));
        break;
    case bondVTUField::damageFactor:
        writer.addPointData("damageFactor", collectValues(bonds, [](const bond& value) { return value.damageFactor(); }));
        break;
    case bondVTUField::modeICriticalEnergy:
        writer.addPointData("modeICriticalEnergy", collectValues(bonds, [](const bond& value) { return value.modeICriticalEnergy(); }));
        break;
    case bondVTUField::modeIICriticalEnergy:
        writer.addPointData("modeIICriticalEnergy", collectValues(bonds, [](const bond& value) { return value.modeIICriticalEnergy(); }));
        break;
    case bondVTUField::modeMixityExponent:
        writer.addPointData("modeMixityExponent", collectValues(bonds, [](const bond& value) { return value.modeMixityExponent(); }));
        break;
    case bondVTUField::damageInitiationRatio:
        writer.addPointData("damageInitiationRatio", collectValues(bonds, [](const bond& value) { return value.damageInitiationRatio(); }));
        break;
    case bondVTUField::maximumEnergyReleaseRatio:
        writer.addPointData("maximumEnergyReleaseRatio", collectValues(bonds, [](const bond& value) { return value.maximumEnergyReleaseRatio(); }));
        break;
    }
}

template <class Storage> void setInteractionDefaults(vtuWriter& writer, const Storage& interactions)
{
    writer.setPoints(collectValues(interactions, [](const typename Storage::host_type& value) { return value.point(); }));
    writer.setVertexCells();
    writer.addPointData("normal", collectValues(interactions, [](const typename Storage::host_type& value) { return value.normal(); }));
}

template <class MasterStorage, class SlaveStorage>
vtuWriter makeBondVTUWriterImpl(std::string fileName, const bondContainer& bonds, const MasterStorage& masterParticles, const SlaveStorage& slaveParticles, const std::vector<bondVTUField>& fields)
{
    vtuWriter writer(std::move(fileName));
    setInteractionDefaults(writer, bonds);
    for (bondVTUField field : fields)
    {
        addBondField(writer, bonds, masterParticles, slaveParticles, field);
    }
    return writer;
}

} // namespace

vtuWriter makeSphereVTUWriter(std::string fileName, const particleContainer& spheres, const Vec3& gravity, const std::vector<sphereVTUField>& fields)
{
    vtuWriter writer(std::move(fileName));
    writer.setPoints(collectValues(spheres, [](const particle& value) { return value.position(); }));
    writer.setVertexCells();
    writer.addPointData("radius", collectValues(spheres, [](const particle& value) { return value.radius(); }));
    writer.addPointData("velocity", collectValues(spheres, [](const particle& value) { return value.velocity(); }));
    for (sphereVTUField field : fields)
    {
        addSphereField(writer, spheres, field, gravity);
    }
    return writer;
}

namespace
{

vtuWriter makeLSParticleVTUWriterImpl(std::string fileName, const LSParticleContainer& particles, const Vec3& gravity, bool fixed, const std::vector<LSParticleVTUField>& fields)
{
    const auto& particleHost = particles.host();
    const int particleCount = static_cast<int>(particleHost.size());
    const std::vector<int> pointOffsets = surfacePointOffsets(particles, fixed);
    std::vector<int> cellOffsets(particleCount + 1, 0);
    std::vector<int> connectivityOffsets(particleCount + 1, 0);
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const LSParticle& particle = particleHost[particleIndex];
        if (!selectLSParticle(particle, fixed))
        {
            continue;
        }
        const auto& surfaceTriangles = particle.surfaceTriangles();
        const int surfaceNodeCount = static_cast<int>(particle.surfaceNodes().size());
        const int surfaceTriangleCount = static_cast<int>(surfaceTriangles.size());
        cellOffsets[particleIndex + 1] = surfaceTriangles.empty() ? surfaceNodeCount : surfaceTriangleCount;
        connectivityOffsets[particleIndex + 1] = surfaceTriangles.empty() ? surfaceNodeCount : 3 * surfaceTriangleCount;
    }
    std::partial_sum(cellOffsets.begin(), cellOffsets.end(), cellOffsets.begin());
    std::partial_sum(connectivityOffsets.begin(), connectivityOffsets.end(), connectivityOffsets.begin());

    std::vector<Vec3> points(pointOffsets.back());
    std::vector<Vec3> velocities(pointOffsets.back());
    std::vector<int> connectivity(connectivityOffsets.back());
    std::vector<int> offsets(cellOffsets.back());
    std::vector<unsigned char> types(cellOffsets.back());
#if defined(_OPENMP)
#pragma omp parallel for schedule(guided, 4) if (pointOffsets.back() >= parallelSurfacePointThreshold)
#endif
    for (int particleIndex = 0; particleIndex < particleCount; ++particleIndex)
    {
        const LSParticle& particle = particleHost[particleIndex];
        if (!selectLSParticle(particle, fixed))
        {
            continue;
        }
        const auto& surfaceNodes = particle.surfaceNodes();
        const auto& surfaceTriangles = particle.surfaceTriangles();
        const int pointOffset = pointOffsets[particleIndex];
        const int cellOffset = cellOffsets[particleIndex];
        const int connectivityOffset = connectivityOffsets[particleIndex];
        for (int surfaceNodeIndex = 0; surfaceNodeIndex < static_cast<int>(surfaceNodes.size()); ++surfaceNodeIndex)
        {
            const Vec3 worldOffset = math::rotateUnit(particle.orientation(), surfaceNodes[surfaceNodeIndex].position_);
            points[pointOffset + surfaceNodeIndex] = particle.position() + worldOffset;
            velocities[pointOffset + surfaceNodeIndex] = particle.velocity() + math::cross(particle.angularVelocity(), worldOffset);
        }

        if (surfaceTriangles.empty())
        {
            for (int surfaceNodeIndex = 0; surfaceNodeIndex < static_cast<int>(surfaceNodes.size()); ++surfaceNodeIndex)
            {
                connectivity[connectivityOffset + surfaceNodeIndex] = pointOffset + surfaceNodeIndex;
                offsets[cellOffset + surfaceNodeIndex] = connectivityOffset + surfaceNodeIndex + 1;
                types[cellOffset + surfaceNodeIndex] = 1;
            }
        }
        else
        {
            for (int triangleIndex = 0; triangleIndex < static_cast<int>(surfaceTriangles.size()); ++triangleIndex)
            {
                const int3& triangle = surfaceTriangles[triangleIndex].connectivity_;
                const int triangleConnectivityOffset = connectivityOffset + 3 * triangleIndex;
                connectivity[triangleConnectivityOffset] = pointOffset + triangle.x;
                connectivity[triangleConnectivityOffset + 1] = pointOffset + triangle.y;
                connectivity[triangleConnectivityOffset + 2] = pointOffset + triangle.z;
                offsets[cellOffset + triangleIndex] = triangleConnectivityOffset + 3;
                types[cellOffset + triangleIndex] = 5;
            }
        }
    }

    vtuWriter writer(std::move(fileName));
    writer.setPoints(points);
    writer.setCells(connectivity, offsets, types);
    writer.addPointData("velocity", velocities);
    for (LSParticleVTUField field : fields)
    {
        addLSParticleField(writer, particles, fixed, field, gravity);
    }
    return writer;
}

} // namespace

vtuWriter makeLSParticleVTUWriter(std::string fileName, const LSParticleContainer& particles, const Vec3& gravity, const std::vector<LSParticleVTUField>& fields)
{
    return makeLSParticleVTUWriterImpl(std::move(fileName), particles, gravity, false, fields);
}

vtuWriter makeFixedLSParticleVTUWriter(std::string fileName, const LSParticleContainer& particles, const Vec3& gravity, const std::vector<LSParticleVTUField>& fields)
{
    return makeLSParticleVTUWriterImpl(std::move(fileName), particles, gravity, true, fields);
}

vtuWriter makeSPHParticleVTUWriter(std::string fileName, const SPHParticleContainer& particles, const Vec3& gravity, const std::vector<SPHParticleVTUField>& fields)
{
    vtuWriter writer(std::move(fileName));
    writer.setPoints(collectValues(particles, [](const SPHParticle& value) { return value.position(); }));
    writer.setVertexCells();
    writer.addPointData("velocity", collectValues(particles, [](const SPHParticle& value) { return value.velocity(); }));
    writer.addPointData("mass", collectValues(particles, [](const SPHParticle& value) { return value.mass(); }));
    writer.addPointData("density", collectValues(particles, [](const SPHParticle& value) { return value.density(); }));
    for (SPHParticleVTUField field : fields)
    {
        addSPHParticleField(writer, particles, field, gravity);
    }
    return writer;
}

vtuWriter makeContactVTUWriter(std::string fileName, const contactContainer& contacts, const std::vector<contactVTUField>& fields)
{
    vtuWriter writer(std::move(fileName));
    setInteractionDefaults(writer, contacts);
    for (contactVTUField field : fields)
    {
        addContactField(writer, contacts, field);
    }
    return writer;
}

vtuWriter makeBondVTUWriter(std::string fileName, const bondContainer& bonds, const particleContainer& spheres, const std::vector<bondVTUField>& fields)
{
    return makeBondVTUWriterImpl(std::move(fileName), bonds, spheres, spheres, fields);
}

vtuWriter makeBondVTUWriter(std::string fileName, const bondContainer& bonds, const LSParticleContainer& LSParticles, const std::vector<bondVTUField>& fields)
{
    return makeBondVTUWriterImpl(std::move(fileName), bonds, LSParticles, LSParticles, fields);
}

vtuWriter makeBondVTUWriter(std::string fileName,
                            const bondContainer& bonds,
                            const particleContainer& masterSpheres,
                            const LSParticleContainer& slaveLSParticles,
                            const std::vector<bondVTUField>& fields)
{
    return makeBondVTUWriterImpl(std::move(fileName), bonds, masterSpheres, slaveLSParticles, fields);
}

} // namespace fundem
