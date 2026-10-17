#include "periodicShear.h"

#include "data/energyOutput.h"
#include "execution/cpu/particleFunctions.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace fundem::validation
{

namespace
{

/** Rejects an unsupported input size before converting it to a container index. */
int sourceCount(const std::vector<initialSphere>& initial)
{
    if (initial.empty() || initial.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Periodic shear requires a non-empty, integer-indexable source configuration.");
    }
    return static_cast<int>(initial.size());
}

} // namespace

periodicShear::periodicShear(const std::vector<initialSphere>& initial, const Vec3& minimum, const Vec3& maximum, shearParameters parameters)
    : SphereDEM(executionMode::CPU), parameters_(parameters), contactForce_(minimum, maximum, initial.empty() ? 0.0 : initial.front().radius_, sourceCount(initial), parameters)
{
    if (!math::isFinite(parameters_.wallSpeed_) || parameters_.wallSpeed_ < 0.0)
    {
        throw std::invalid_argument("Periodic shear requires a finite non-negative wall speed.");
    }

    const Real radius = initial.front().radius_;
    sourceIds_.reserve(initial.size());
    wallKinds_.reserve(initial.size());
    wallZ_.reserve(initial.size());
    mutableSpheres().host().reserve(initial.size());

    // Material supplies rho_p = 1 and the sphere mass/inertia only. The
    // external-force component owns every spring, dashpot and friction value.
    const material grainMaterial(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
    const int materialIndex = addMaterial(grainMaterial);
    std::unordered_set<int> sourceIds;
    sourceIds.reserve(initial.size());
    for (const initialSphere& source : initial)
    {
        if (!math::isFinite(source.position_) || !math::isFinite(source.radius_) || !math::nearlyEqual(source.radius_, radius) || (source.wall_ != -1 && source.wall_ != 0 && source.wall_ != 1) ||
            !sourceIds.insert(source.sourceId_).second)
        {
            throw std::invalid_argument("Periodic shear requires finite, monodisperse source spheres with distinct IDs and valid wall roles.");
        }
        particle sphere;
        sphere.setMaterial(materials(), materialIndex);
        sphere.setRadius(radius);
        sphere.setPosition(wrapPosition(source.position_));
        sphere.setVelocity({source.wall_ * parameters_.wallSpeed_, 0.0, 0.0});
        if (source.wall_ != 0)
        {
            sphere.setInfiniteMass();
        }
        addSphere(sphere);
        sourceIds_.push_back(source.sourceId_);
        wallKinds_.push_back(source.wall_);
        wallZ_.push_back(source.position_.z);
    }

    setBoundary(minimum, maximum);
    setGravity(Vec3::zero());
    setTimeStep(7.2320e-5);
}

Vec3 periodicShear::wrapPosition(Vec3 position) const noexcept { return contactForce_.wrapPosition(position); }

Vec3 periodicShear::minimumImage(Vec3 displacement) const noexcept { return contactForce_.minimumImage(displacement); }

void periodicShear::validateExecutionMode(executionMode value) const
{
    if (value != executionMode::CPU)
    {
        throw std::invalid_argument("The periodic-shear validation implements only the CPU backend.");
    }
}

void periodicShear::validateSystemConfiguration() const
{
    if (spheres().hostSize() != sourceIds_.size() || LSParticles().hostSize() != 0 || sphereInteractions().bonds().hostSize() != 0 || sphereLSInteractions().bonds().hostSize() != 0 ||
        LSParticleInteractions().bonds().hostSize() != 0)
    {
        throw std::invalid_argument("The periodic-shear validation population cannot be extended with other particles or bonds.");
    }
    if (gravity() != Vec3::zero() || minimumBoundary() != contactForce_.minimumBoundary() || maximumBoundary() != contactForce_.maximumBoundary())
    {
        throw std::invalid_argument("The periodic-shear validation requires zero gravity and its original periodic domain.");
    }
}

void periodicShear::prepareKinematics()
{
    auto& host = mutableSpheres().host();
    for (int index = 0; index < static_cast<int>(host.size()); ++index)
    {
        particle& sphere = host[index];
        Vec3 position = wrapPosition(sphere.position());
        if (wallKinds_[index] != 0)
        {
            position.z = wallZ_[index];
            sphere.setVelocity({wallKinds_[index] * parameters_.wallSpeed_, 0.0, 0.0});
            sphere.setAngularVelocity(Vec3::zero());
        }
        sphere.setPosition(position);
    }
}

void periodicShear::assembleContacts(Real historyTimeStep)
{
    cpu::clearForceAndTorque(mutableSpheres());
    historyTimeStep_ = historyTimeStep;
    addSphereExternalForceAndTorque(mutableSpheres().host());
}

void periodicShear::addSphereExternalForceAndTorque(particleContainer::host_container_type& particles) { contactForce_.apply(particles, historyTimeStep_); }

void periodicShear::initialize(cpu::mode)
{
    prepareKinematics();
    assembleContacts(0.0);
}

void periodicShear::initialize(gpu::mode, cudaStream_t) { throw std::logic_error("Periodic shear has no GPU implementation."); }

void periodicShear::initialize(hybrid::mode, cudaStream_t) { throw std::logic_error("Periodic shear has no Hybrid implementation."); }

void periodicShear::advanceCPU(Real timeStep)
{
    const Real halfTimeStep = 0.5 * timeStep;
    cpu::integrateVelocity(mutableSpheres(), Vec3::zero(), halfTimeStep);
    prepareKinematics();
    cpu::integratePosition(mutableSpheres(), timeStep);
    prepareKinematics();
    assembleContacts(timeStep);
    cpu::integrateVelocity(mutableSpheres(), Vec3::zero(), halfTimeStep);
}

void periodicShear::advanceGPU(Real) { throw std::logic_error("Periodic shear has no GPU implementation."); }

void periodicShear::advanceHybrid(Real) { throw std::logic_error("Periodic shear has no Hybrid implementation."); }

energyRecord periodicShear::systemEnergy() const
{
    energyRecord result;
    addParticleEnergy(result, spheres(), Vec3::zero());
    for (const shearContact& contact : contacts())
    {
        result.contactNormalElasticEnergy_ += 0.5 * parameters_.normalStiffness_ * contact.overlap_ * contact.overlap_;
        result.contactSlidingElasticEnergy_ += 0.5 * parameters_.slidingStiffness_ * math::normSquared(contact.tangentialSpring_);
    }
    return result;
}

} // namespace fundem::validation
