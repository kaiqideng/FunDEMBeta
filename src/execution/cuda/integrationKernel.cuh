/**
 * @file integrationKernel.cuh
 * @brief Declares rigid-body device views and streamed integration launch interfaces.
 */
#pragma once

#include "particle/rigidBody.h"

#include <cuda_runtime_api.h>

namespace fundem::cuda
{
namespace detail
{

/** Raw rigid-body arrays consumed by integration kernels. */
struct integrationDeviceView {
    math::Vec3* position_{nullptr};                   ///< World positions.
    math::Quaternion* orientation_{nullptr};          ///< Unit orientations.
    math::Vec3* velocity_{nullptr};                   ///< Linear velocities.
    math::Vec3* angularVelocity_{nullptr};            ///< Angular velocities.
    const math::Vec3* force_{nullptr};                ///< Accumulated forces.
    const math::Vec3* torque_{nullptr};               ///< Accumulated torques.
    const math::Real* inverseMass_{nullptr};          ///< Reciprocal masses.
    const math::Mat3* inertiaTensor_{nullptr};        ///< Body-frame inertias.
    const math::Mat3* inverseInertiaTensor_{nullptr}; ///< Inverse body-frame inertias.
};

/** Launches velocity integration on @p stream without synchronizing. */
void launchVelocityAndAngularVelocityIntegration(integrationDeviceView particles, const math::Vec3& gravity, math::Real timeStep, int particleCount, cudaStream_t stream);
/** Launches position/orientation integration on @p stream without synchronizing. */
void launchPositionAndOrientationIntegration(integrationDeviceView particles, math::Real timeStep, int particleCount, cudaStream_t stream);

template <class ParticleStorage> integrationDeviceView makeIntegrationDeviceView(ParticleStorage& particles) noexcept
{
    integrationDeviceView view;
    view.position_ = particles.template device<rigidBody::positionField>();
    view.orientation_ = particles.template device<rigidBody::orientationField>();
    view.velocity_ = particles.template device<rigidBody::velocityField>();
    view.angularVelocity_ = particles.template device<rigidBody::angularVelocityField>();
    view.force_ = particles.template device<rigidBody::forceField>();
    view.torque_ = particles.template device<rigidBody::torqueField>();
    view.inverseMass_ = particles.template device<rigidBody::inverseMassField>();
    view.inertiaTensor_ = particles.template device<rigidBody::inertiaTensorField>();
    view.inverseInertiaTensor_ = particles.template device<rigidBody::inverseInertiaTensorField>();
    return view;
}

} // namespace detail

/** Launches velocity integration for an entire device container. */
template <class ParticleStorage> void launchVelocityAndAngularVelocityIntegration(ParticleStorage& particles, const math::Vec3& gravity, math::Real timeStep, cudaStream_t stream = nullptr)
{
    detail::launchVelocityAndAngularVelocityIntegration(detail::makeIntegrationDeviceView(particles), gravity, timeStep, static_cast<int>(particles.deviceSize()), stream);
}

/** Launches position and orientation integration for an entire device container. */
template <class ParticleStorage> void launchPositionAndOrientationIntegration(ParticleStorage& particles, math::Real timeStep, cudaStream_t stream = nullptr)
{
    detail::launchPositionAndOrientationIntegration(detail::makeIntegrationDeviceView(particles), timeStep, static_cast<int>(particles.deviceSize()), stream);
}

} // namespace fundem::cuda
