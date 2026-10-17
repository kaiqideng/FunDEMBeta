/**
 * @file periodicContactForceTest.cpp
 * @brief Output-free checks of additive periodic contact loads and the solver hook.
 */
#include "periodicContactForce.h"
#include "periodicShear.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using namespace fundem;
using namespace fundem::validation;
using Real = math::Real;
using Vec3 = math::Vec3;

constexpr Vec3 lower{0.0, 0.0, 0.0};
constexpr Vec3 upper{10.0, 8.0, 6.0};
constexpr Real radius = 0.5;
constexpr Real stepSize = 1.0e-5;

/** Throws independently of build-time assertion settings. */
void require(bool condition, const std::string& message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

/** Checks one finite scalar against a known analytical result. */
void near(Real actual, Real expected, Real tolerance, const std::string& message) { require(math::isFinite(actual) && std::abs(actual - expected) <= tolerance, message); }

/** Checks one finite vector against a known analytical result. */
void near(const Vec3& actual, const Vec3& expected, Real tolerance, const std::string& message) { require(math::isFinite(actual) && math::norm(actual - expected) <= tolerance, message); }

/** Returns direct, independently specified spring and dashpot coefficients. */
contactForceParameters forceParameters()
{
    contactForceParameters result;
    result.normalStiffness_ = 10000.0;
    result.slidingStiffness_ = 4000.0;
    result.normalDamping_ = 7.0;
    result.slidingDamping_ = 5.0;
    result.friction_ = 0.3;
    return result;
}

/** Owns two seam-contacting particles with a material used for mass only. */
struct pairFixture {
    materialContainer materials_; ///< Material lifetime exceeds all particle snapshots.
    particleContainer particles_; ///< Stable source order; no solver initialization is needed.

    pairFixture()
    {
        materials_.host().emplace_back(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0);
        for (const Vec3& position : {Vec3{0.2, 3.0, 3.0}, Vec3{9.3, 3.0, 3.0}})
        {
            particle sphere;
            sphere.setMaterial(materials_, 0);
            sphere.setRadius(radius);
            sphere.setPosition(position);
            particles_.host().push_back(sphere);
        }
    }

    /** Clears only test loads so repeated external evaluations remain additive. */
    void clearLoads()
    {
        for (auto& sphere : particles_.host())
        {
            sphere.setForce(Vec3::zero());
            sphere.setTorque(Vec3::zero());
        }
    }
};

/** Checks exact dashpots, additive loads, periodic lever arms and action/reaction. */
void additiveContactLoads()
{
    pairFixture fixture;
    auto& host = fixture.particles_.host();
    const Vec3 firstLoad{1.0, 2.0, 3.0};
    const Vec3 secondLoad{-4.0, 5.0, -6.0};
    const Vec3 firstMoment{0.1, -0.2, 0.3};
    const Vec3 secondMoment{-0.4, 0.5, -0.6};
    host[0].setVelocity({-2.0, 0.3, 0.0});
    host[0].setForce(firstLoad);
    host[1].setForce(secondLoad);
    host[0].setTorque(firstMoment);
    host[1].setTorque(secondMoment);
    periodicContactForce contactForce(lower, upper, radius, 2, forceParameters());
    contactForce.apply(host, 0.0);

    require(contactForce.contacts().size() == 1, "External seam contact was omitted or duplicated.");
    const auto& contact = contactForce.contacts().front();
    require(contact.first_ == 0 && contact.second_ == 1, "External contact changed stable endpoint indices.");
    near(contact.overlap_, 0.1, 1.0e-12, "External seam overlap is incorrect.");
    near(contact.normalForce_, 1014.0, 1.0e-8, "The external normal dashpot was mass-scaled or combined with material stiffness.");
    near(contact.tangentialForce_, {0.0, -1.5, 0.0}, 1.0e-12, "The external tangential dashpot did not use ct times relative velocity.");
    const Vec3 contribution{1014.0, -1.5, 0.0};
    near(host[0].force(), firstLoad + contribution, 1.0e-8, "External contact erased the first preload or added its force twice.");
    near(host[1].force(), secondLoad - contribution, 1.0e-8, "External contact erased the second preload or added its force twice.");
    near(host[0].torque(), firstMoment + Vec3{0.0, 0.0, 0.675}, 1.0e-12, "External contact erased the first torque preload or used the wrong arm.");
    near(host[1].torque(), secondMoment + Vec3{0.0, 0.0, 0.675}, 1.0e-12, "External contact erased the second torque preload or used the wrong arm.");
    near(host[0].force() - firstLoad + host[1].force() - secondLoad, Vec3::zero(), 1.0e-10, "External contact violates Newton's third law.");
    const Vec3 separation = contactForce.minimumImage(host[1].position() - host[0].position());
    const Vec3 momentBalance = host[0].torque() - firstMoment + host[1].torque() - secondMoment + math::cross(separation, host[1].force() - secondLoad);
    near(momentBalance, Vec3::zero(), 1.0e-10, "External periodic contact violates angular momentum balance.");
    near(host[0].getMaterial().normalStiffness(), 0.0, 0.0, "The external law unexpectedly changed material properties.");
}

/** Checks that a separating dashpot never creates attraction. */
void nonAttractiveNormalForce()
{
    pairFixture fixture;
    auto& host = fixture.particles_.host();
    host[0].setVelocity({1000.0, 0.0, 0.0});
    periodicContactForce contactForce(lower, upper, radius, 2, forceParameters());
    contactForce.apply(host, 0.0);
    near(contactForce.contacts().front().normalForce_, 0.0, 0.0, "The external normal dashpot produced an attractive force.");
    near(host[0].force(), Vec3::zero(), 1.0e-12, "A clamped normal contact left a nonzero resultant.");
}

/** Checks the explicit Mercury release variant without changing the default law. */
void signedNormalRelease()
{
    pairFixture fixture;
    auto& host = fixture.particles_.host();
    host[0].setVelocity({1000.0, 500.0, 0.0});
    auto coefficients = forceParameters();
    coefficients.normalRelease_ = normalReleasePolicy::retainUntilSeparation;
    periodicContactForce contactForce(lower, upper, radius, 2, coefficients);
    contactForce.apply(host, 0.0);
    const auto& contact = contactForce.contacts().front();
    near(contact.normalForce_, -6000.0, 1.0e-8, "The signed-release option still truncated the unloading dashpot force.");
    near(contact.tangentialForce_, {0.0, -1800.0, 0.0}, 1.0e-8, "Signed-release friction did not use the absolute normal force.");
    near(contact.tangentialSpring_, {0.0, -0.175, 0.0}, 1.0e-12, "The Coulomb return map is inconsistent with the complete tangential force.");
    near(host[0].force(), {-6000.0, -1800.0, 0.0}, 1.0e-8, "The signed force was not applied to the particle.");
    near(host[0].force() + host[1].force(), Vec3::zero(), 1.0e-10, "The signed release variant violates action/reaction.");
}

/** Checks history reuse at zero elapsed time, across seams and after separation. */
void historyLifecycle()
{
    pairFixture fixture;
    auto& host = fixture.particles_.host();
    host[0].setVelocity({0.0, 0.3, 0.0});
    periodicContactForce contactForce(lower, upper, radius, 2, forceParameters());
    contactForce.apply(host, 1.0e-4);
    const Vec3 spring = contactForce.contacts().front().tangentialSpring_;
    const Vec3 force = contactForce.contacts().front().tangentialForce_;
    near(spring, {0.0, 3.0e-5, 0.0}, 1.0e-14, "External history advanced by an incorrect interval.");

    fixture.clearLoads();
    contactForce.apply(host, 0.0);
    near(contactForce.contacts().front().tangentialSpring_, spring, 1.0e-14, "A zero-time external evaluation advanced contact history.");
    near(contactForce.contacts().front().tangentialForce_, force, 1.0e-12, "A zero-time external evaluation changed contact force.");

    for (auto& sphere : host)
    {
        sphere.setPosition(contactForce.wrapPosition(sphere.position() + Vec3{9.7, 7.8, 0.0}));
    }
    fixture.clearLoads();
    contactForce.apply(host, 0.0);
    near(contactForce.contacts().front().tangentialSpring_, spring, 1.0e-14, "Crossing periodic images discarded external contact history.");
    near(contactForce.contacts().front().tangentialForce_, force, 1.0e-12, "Crossing periodic images changed external tangential force.");

    host[0].setPosition({2.0, 2.0, 3.0});
    host[1].setPosition({5.0, 5.0, 3.0});
    fixture.clearLoads();
    contactForce.apply(host, 0.0);
    require(contactForce.contacts().empty(), "Separated particles retained an external contact.");
    host[0].setPosition({0.2, 3.0, 3.0});
    host[1].setPosition({9.3, 3.0, 3.0});
    host[0].setVelocity(Vec3::zero());
    contactForce.apply(host, 0.0);
    require(contactForce.contacts().size() == 1, "External contact did not reopen after separation.");
    near(contactForce.contacts().front().tangentialSpring_, Vec3::zero(), 1.0e-14, "External recontact reused a terminated tangential history.");
}

/** Checks Coulomb return and the zero-stiffness, pure tangential dashpot case. */
void tangentialLimitAndPureDashpot()
{
    pairFixture fixture;
    auto& host = fixture.particles_.host();
    auto properties = forceParameters();
    properties.slidingDamping_ = 100.0;
    host[0].setVelocity({0.0, 1000.0, 0.0});
    periodicContactForce contactForce(lower, upper, radius, 2, properties);
    contactForce.apply(host, 0.0);
    const auto first = contactForce.contacts().front();
    near(math::norm(first.tangentialForce_), properties.friction_ * first.normalForce_, 1.0e-8, "External spring-plus-dashpot force escaped the Coulomb disk.");
    fixture.clearLoads();
    contactForce.apply(host, 0.0);
    near(contactForce.contacts().front().tangentialSpring_, first.tangentialSpring_, 1.0e-10, "A returned sliding spring changed during a zero-time evaluation.");
    near(contactForce.contacts().front().tangentialForce_, first.tangentialForce_, 1.0e-8, "The Coulomb return was inconsistent with its saved spring.");

    fixture.clearLoads();
    properties.slidingStiffness_ = 0.0;
    properties.slidingDamping_ = 5.0;
    host[0].setVelocity({0.0, 0.3, 0.0});
    periodicContactForce dashpot(lower, upper, radius, 2, properties);
    dashpot.apply(host, 1.0e-4);
    near(dashpot.contacts().front().tangentialForce_, {0.0, -1.5, 0.0}, 1.0e-12, "The zero-stiffness tangential dashpot was not evaluated directly.");
    near(dashpot.contacts().front().tangentialSpring_, Vec3::zero(), 0.0, "The pure dashpot retained a fictitious elastic spring.");
}

/** Confirms that fixed roles follow inverse mass, not a separate wall array. */
void infiniteMassEndpoint()
{
    pairFixture fixture;
    auto& host = fixture.particles_.host();
    host[0].setInfiniteMass();
    periodicContactForce contactForce(lower, upper, radius, 2, forceParameters());
    contactForce.apply(host, 0.0);
    require(contactForce.contacts().size() == 1, "External finite-fixed contact was lost.");
    const auto& contact = contactForce.contacts().front();
    require(contact.first_ == 1 && contact.second_ == 0, "The infinite-mass endpoint was not placed second.");
    near(contact.normalForce_, 1000.0, 1.0e-8, "External finite-fixed stiffness was rescaled.");
    near(host[0].force() + host[1].force(), Vec3::zero(), 1.0e-10, "External fixed reaction is not equal and opposite.");
    host[1].setInfiniteMass();
    fixture.clearLoads();
    contactForce.apply(host, 0.0);
    require(contactForce.contacts().empty(), "External law created an infinite-mass/infinite-mass contact.");
}

/** Proves an external boundary provider is used without changing the standalone law. */
void independentBoundaryCallbacks()
{
    pairFixture standaloneFixture;
    pairFixture callbackFixture;
    standaloneFixture.particles_.host()[0].setVelocity({-2.0, 0.3, 0.0});
    callbackFixture.particles_.host()[0].setVelocity({-2.0, 0.3, 0.0});
    periodicContactForce standalone(lower, upper, radius, 2, forceParameters());
    struct probe {
        const periodicContactForce* geometry_;
        mutable int wrapCalls_{0};
        mutable int imageCalls_{0};
    } boundary{&standalone};
    const periodicBoundaryFunctions functions{
        &boundary,
        [](const void* context, Vec3 position) noexcept
        {
            const auto& service = *static_cast<const probe*>(context);
            ++service.wrapCalls_;
            return service.geometry_->wrapPosition(position);
        },
        [](const void* context, Vec3 displacement) noexcept
        {
            const auto& service = *static_cast<const probe*>(context);
            ++service.imageCalls_;
            return service.geometry_->minimumImage(displacement);
        }};
    periodicContactForce external(lower, upper, radius, 2, forceParameters(), functions);
    for (Real timeStep : {Real{1.0e-4}, Real{0.0}})
    {
        const int previousWrapCalls = boundary.wrapCalls_;
        const int previousImageCalls = boundary.imageCalls_;
        standaloneFixture.clearLoads();
        callbackFixture.clearLoads();
        standalone.apply(standaloneFixture.particles_.host(), timeStep);
        external.apply(callbackFixture.particles_.host(), timeStep);
        require(boundary.wrapCalls_ > previousWrapCalls && boundary.imageCalls_ > previousImageCalls, "Pair evaluation bypassed the external wrapping or minimum-image provider.");
        require(standalone.contacts().size() == 1 && external.contacts().size() == 1, "The boundary callback changed the physical seam-contact count.");
        const auto& expected = standalone.contacts().front();
        const auto& actual = external.contacts().front();
        near(actual.point_, expected.point_, 0.0, "The external boundary changed the standalone contact point.");
        near(actual.normal_, expected.normal_, 0.0, "The external boundary changed the standalone contact normal.");
        near(actual.normalForce_, expected.normalForce_, 0.0, "The external boundary changed the normal-force law.");
        near(actual.tangentialForce_, expected.tangentialForce_, 0.0, "The external boundary changed tangential force.");
        near(actual.tangentialSpring_, expected.tangentialSpring_, 0.0, "The external boundary changed tangential-history advancement.");
        for (int index = 0; index < 2; ++index)
        {
            near(callbackFixture.particles_.host()[index].force(), standaloneFixture.particles_.host()[index].force(), 0.0, "External boundary callbacks changed the additive load.");
            near(callbackFixture.particles_.host()[index].torque(), standaloneFixture.particles_.host()[index].torque(), 0.0, "External boundary callbacks changed the periodic lever arm.");
        }
    }
    near(external.wrapPosition({20.2, -5.0, 3.0}), standalone.wrapPosition({20.2, -5.0, 3.0}), 0.0, "Explicit wrapping did not use the provider's geometry.");
    bool rejected = false;
    try
    {
        periodicContactForce incomplete(lower, upper, radius, 2, forceParameters(), {&boundary, functions.wrapPosition_, nullptr});
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "A partial boundary provider silently mixed external and built-in periodic geometry.");
}

/** Checks that linked-cell search follows an external direction selector, including Z. */
void selectableBoundaryCallbacks()
{
    struct selectedBoundary {
        int mask_;
    } boundary{1};
    const periodicBoundaryFunctions functions{
        &boundary,
        [](const void* context, Vec3 position) noexcept
        {
            const auto& service = *static_cast<const selectedBoundary*>(context);
            for (int axis = 0; axis < 3; ++axis)
                if (service.mask_ & (1 << axis))
                {
                    position[axis] = std::fmod(position[axis] - lower[axis], upper[axis] - lower[axis]);
                    if (position[axis] < 0.0) position[axis] += upper[axis] - lower[axis];
                    position[axis] += lower[axis];
                }
            return position;
        },
        [](const void* context, Vec3 displacement) noexcept
        {
            const auto& service = *static_cast<const selectedBoundary*>(context);
            for (int axis = 0; axis < 3; ++axis)
                if (service.mask_ & (1 << axis)) displacement[axis] -= (upper[axis] - lower[axis]) * std::round(displacement[axis] / (upper[axis] - lower[axis]));
            return displacement;
        },
        [](const void* context, int axis) noexcept { return (static_cast<const selectedBoundary*>(context)->mask_ & (1 << axis)) != 0; }};
    for (int axis = 0; axis < 3; ++axis)
    {
        pairFixture fixture;
        auto& host = fixture.particles_.host();
        Vec3 first{3.0, 3.0, 3.0};
        Vec3 second = first;
        first[axis] = 0.2;
        second[axis] = upper[axis] - 0.7;
        host[0].setPosition(first);
        host[1].setPosition(second);
        boundary.mask_ = 1 << axis;
        periodicContactForce enabled(lower, upper, radius, 2, forceParameters(), functions);
        enabled.apply(host, 0.0);
        require(enabled.contacts().size() == 1 && enabled.contacts().front().first_ == 0 && enabled.contacts().front().second_ == 1,
                "An enabled single-axis seam was lost, duplicated or replaced by a self-contact.");
        near(enabled.contacts().front().overlap_, 0.1, 1.0e-12, "Selected-axis periodic overlap is incorrect.");
        fixture.clearLoads();
        boundary.mask_ = 1 << ((axis + 1) % 3);
        periodicContactForce disabled(lower, upper, radius, 2, forceParameters(), functions);
        disabled.apply(host, 0.0);
        require(disabled.contacts().empty(), "Linked-cell search wrapped a nonperiodic direction.");
    }
    boundary.mask_ = 0;
    bool rejected = false;
    try { periodicContactForce allOff(lower, upper, radius, 2, forceParameters(), functions); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "The external direction selector accepted an all-disabled periodic boundary.");
}

/** Exposes kinematic setup and counts real calls through the existing external-force hook. */
class countingShear : public periodicShear
{
public:
    using periodicShear::periodicShear;

    int externalCalls() const noexcept { return externalCalls_; }

    void velocity(int index, const Vec3& value) { mutableSpheres().host().at(index).setVelocity(value); }

protected:
    /** Preserves the adapter's external law while proving virtual dispatch is used. */
    void addSphereExternalForceAndTorque(particleContainer::host_container_type& host) override
    {
        ++externalCalls_;
        periodicShear::addSphereExternalForceAndTorque(host);
    }

private:
    int externalCalls_{0}; ///< Number of actual force assembly passes, not solver API calls.
};

/** Ensures the adapter disables ordinary forces and advances external history once per step. */
void solverExternalHook()
{
    shearParameters properties;
    const auto direct = forceParameters();
    properties.normalStiffness_ = direct.normalStiffness_;
    properties.slidingStiffness_ = direct.slidingStiffness_;
    properties.normalDamping_ = direct.normalDamping_;
    properties.slidingDamping_ = direct.slidingDamping_;
    properties.friction_ = direct.friction_;
    properties.wallSpeed_ = 0.0;
    const std::vector<initialSphere> initial{{11, {0.2, 3.0, 3.0}, radius, 0}, {29, {9.3, 3.0, 3.0}, radius, 0}};
    countingShear model(initial, lower, upper, properties);
    model.setTimeStep(stepSize);
    model.velocity(0, {-2.0, 0.3, 0.0});
    model.initialize();
    require(model.externalCalls() == 1, "Initialization did not call the external-force hook exactly once.");
    const auto& materials = model.materials().host();
    require(materials.size() == 1, "The adapter unexpectedly created multiple mass-only materials.");
    near(materials.front().normalStiffness(), 0.0, 0.0, "Ordinary normal stiffness was not disabled.");
    near(materials.front().slidingStiffness(), 0.0, 0.0, "Ordinary tangential stiffness was not disabled.");
    near(materials.front().slidingFrictionCoefficient(), 0.0, 0.0, "Ordinary contact friction was not disabled.");
    near(model.contacts().front().normalForce_, 1014.0, 1.0e-8, "The adapter did not use the direct external normal law.");
    near(model.spheres().host()[0].force(), {1014.0, -1.5, 0.0}, 1.0e-8, "Initialization duplicated or omitted the external contact force.");
    near(model.contacts().front().tangentialSpring_, Vec3::zero(), 1.0e-14, "Initialization incorrectly advanced tangential history.");
    model.initialize();
    require(model.externalCalls() == 1, "Repeated initialization duplicated the external assembly.");
    model.step();
    require(model.externalCalls() == 2, "One DEM step did not assemble the external law exactly once.");
    const auto contact = model.contacts().front();
    require(math::norm(contact.tangentialSpring_) > 0.0, "The first DEM step did not advance external contact history.");
    const Vec3 resultant = contact.normalForce_ * contact.normal_ + contact.tangentialForce_;
    near(model.spheres().host()[0].force(), resultant, 1.0e-8, "The adapter added a duplicate ordinary force after external assembly.");

    model.solve(0);
    model.solve(0);
    require(model.externalCalls() == 2, "A zero-step solve continuation unexpectedly reassembled the external law.");
    near(model.contacts().front().tangentialSpring_, contact.tangentialSpring_, 1.0e-14, "A zero-step continuation advanced external history.");
    model.step();
    require(model.externalCalls() == 3 && model.stepCount() == 2, "Resumed stepping bypassed or duplicated the external-force hook.");
    near(model.time(), 2.0 * stepSize, 1.0e-14, "External adapter advancement changed the represented simulation time.");
}

} // namespace

int main()
{
    try
    {
        additiveContactLoads();
        nonAttractiveNormalForce();
        signedNormalRelease();
        historyLifecycle();
        tangentialLimitAndPureDashpot();
        infiniteMassEndpoint();
        independentBoundaryCallbacks();
        selectableBoundaryCallbacks();
        solverExternalHook();
        std::cout << "periodicContactForceTest: PASS (additive loads, direct dashpots, periodic history, external solver hook).\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "periodicContactForceTest: FAIL: " << error.what() << '\n';
        return 1;
    }
}
