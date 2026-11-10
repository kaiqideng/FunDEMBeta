/** Small, output-free checks of the periodic shear benchmark implementation. */
#include "periodicShear.h"

#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <map>
#include <set>
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
constexpr Real stepSize = 1.0e-5;

void require(bool condition, const std::string& message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

void near(Real actual, Real expected, Real tolerance, const std::string& message) { require(math::isFinite(actual) && std::abs(actual - expected) <= tolerance, message); }

void near(const Vec3& actual, const Vec3& expected, Real tolerance, const std::string& message) { require(math::isFinite(actual) && math::norm(actual - expected) <= tolerance, message); }

shearParameters parameters()
{
    shearParameters result;
    result.normalStiffness_ = 10000.0;
    result.slidingStiffness_ = 4000.0;
    result.normalDamping_ = 0.0;
    result.slidingDamping_ = 0.0;
    result.friction_ = 0.3;
    result.wallSpeed_ = 0.5;
    return result;
}

/** Allows controlled coordinate changes without resetting the solver's contact history. */
class testShear : public periodicShear
{
public:
    using periodicShear::periodicShear;

    void translate(const Vec3& displacement)
    {
        for (auto& sphere : mutableSpheres().host())
        {
            sphere.setPosition(wrapPosition(sphere.position() + displacement));
        }
    }

    void place(int index, const Vec3& position) { mutableSpheres().host().at(index).setPosition(position); }

    void velocity(int index, const Vec3& value)
    {
        auto& sphere = mutableSpheres().host().at(index);
        sphere.setVelocity(value);
        sphere.setAngularVelocity(Vec3::zero());
    }
};

std::vector<initialSphere> pair(const Vec3& first, const Vec3& second) { return {{11, first, 0.5, 0}, {29, second, 0.5, 0}}; }

void periodicCoordinates()
{
    testShear model(pair({0.2, 3.0, 3.0}, {9.3, 3.0, 3.0}), lower, upper, parameters());
    near(model.minimumImage({9.1, 0.0, 0.2}), {-0.9, 0.0, 0.2}, 1.0e-12, "X minimum image is incorrect.");
    near(model.minimumImage({0.0, -7.4, -0.2}), {0.0, 0.6, -0.2}, 1.0e-12, "Y minimum image is incorrect.");
    near(model.minimumImage({9.4, 7.4, 6.5}), {-0.6, -0.6, 6.5}, 1.0e-12, "Corner minimum image changed the nonperiodic Z coordinate.");
    near(model.wrapPosition({-0.1, 8.1, -0.3}), {9.9, 0.1, -0.3}, 1.0e-12, "Position wrapping is incorrect.");
    near(model.wrapPosition({10.0, 8.0, 6.0}), {0.0, 0.0, 6.0}, 1.0e-12, "Upper periodic boundaries are not half open.");
    near(model.wrapPosition({32.0, -17.0, 3.0}), {2.0, 7.0, 3.0}, 1.0e-12, "Wrapping multiple box lengths failed.");
}

void normalContact(const Vec3& first, const Vec3& second, const Vec3& maximum = upper)
{
    const auto properties = parameters();
    testShear model(pair(first, second), lower, maximum, properties);
    model.setTimeStep(stepSize);
    model.initialize();
    require(model.sourceIds() == std::vector<int>({11, 29}), "Source-file IDs were reordered or replaced by container indices.");
    require(model.contacts().size() == 1, "A seam pair was omitted or duplicated.");
    const auto& contact = model.contacts().front();
    const auto& spheres = model.spheres().host();
    for (const auto& sphere : spheres)
    {
        near(1.0 / sphere.inverseMass(), math::pi / 6.0, 1.0e-14, "Unit-density sphere mass is not pi/6.");
        for (int axis = 0; axis < 3; ++axis)
        {
            near(sphere.inertiaTensor()(axis, axis), math::pi / 60.0, 1.0e-14, "Unit-density sphere inertia is not pi/60.");
        }
    }
    const Vec3 separation = model.minimumImage(first - second);
    const Real distance = math::norm(separation);
    const Real overlap = 1.0 - distance;
    const Real normalForce = properties.normalStiffness_ * overlap;
    require(contact.first_ == 0 && contact.second_ == 1, "Contact endpoints did not preserve stable particle indices.");
    near(contact.overlap_, overlap, 1.0e-12, "Periodic overlap is incorrect.");
    near(contact.normal_, separation / distance, 1.0e-12, "Periodic contact normal is incorrect.");
    near(contact.normalForce_, normalForce, 1.0e-8, "Normal contact force differs from stiffness times overlap.");
    near(spheres[0].force(), normalForce * contact.normal_, 1.0e-8, "Force has the wrong direction or magnitude.");
    near(spheres[0].force() + spheres[1].force(), Vec3::zero(), 1.0e-10, "Contact violates Newton's third law.");
    near(spheres[0].torque(), Vec3::zero(), 1.0e-10, "A normal seam contact produced a spurious first-particle torque.");
    near(spheres[1].torque(), Vec3::zero(), 1.0e-10, "A normal seam contact produced a spurious second-particle torque.");
}

void normalDashpot()
{
    auto properties = parameters();
    properties.normalDamping_ = 7.0;
    const auto check = [&](Real velocity, Real expectedForce)
    {
        testShear model(pair({0.2, 3.0, 3.0}, {9.3, 3.0, 3.0}), lower, upper, properties);
        model.velocity(0, {velocity, 0.0, 0.0});
        model.initialize();
        require(model.contacts().size() == 1, "Normal dashpot check lost its overlapping pair.");
        near(model.contacts().front().normalForce_, expectedForce, 1.0e-8, "Normal dashpot coefficient or nonattractive force clamp is incorrect.");
    };
    check(-2.0, 1014.0);
    check(1000.0, 0.0);
}

void wallContactCoefficients()
{
    auto properties = parameters();
    properties.wallSpeed_ = 0.0;
    properties.normalDamping_ = 7.0;
    const std::vector<initialSphere> initial{{47, {9.3, 3.0, 3.0}, 0.5, 1}, {71, {0.2, 3.0, 3.0}, 0.5, 0}};
    const auto check = [&](Real velocity, Real expectedForce)
    {
        testShear model(initial, lower, upper, properties);
        model.velocity(1, {velocity, 0.0, 0.0});
        model.initialize();
        require(model.contacts().size() == 1, "A dynamic-fixed seam pair was omitted or duplicated.");
        const auto& contact = model.contacts().front();
        const auto& spheres = model.spheres().host();
        require(contact.first_ == 1 && contact.second_ == 0, "A wall-first input was not normalized to a finite first endpoint.");
        near(spheres[0].inverseMass(), 0.0, 0.0, "A wall grain has finite mass.");
        near(contact.normalForce_, expectedForce, 1.0e-8, "Dynamic-fixed contact rescaled the direct spring/dashpot coefficients.");
        near(spheres[1].force(), {expectedForce, 0.0, 0.0}, 1.0e-8, "Dynamic-fixed seam force has the wrong direction.");
        near(spheres[0].force() + spheres[1].force(), Vec3::zero(), 1.0e-10, "The prescribed wall reaction is not equal and opposite.");
    };
    check(0.0, 1000.0);
    check(-2.0, 1014.0);
}

void tangentialBalance()
{
    auto properties = parameters();
    properties.slidingDamping_ = 5.0;
    testShear model(pair({0.2, 3.0, 3.0}, {9.3, 3.0, 3.0}), lower, upper, properties);
    model.velocity(0, {0.0, 0.3, 0.0});
    model.setTimeStep(stepSize);
    model.initialize();
    const auto& spheres = model.spheres().host();
    const Vec3 separation = model.minimumImage(spheres[1].position() - spheres[0].position());
    near(spheres[0].force() + spheres[1].force(), Vec3::zero(), 1.0e-10, "Tangential force violates Newton's third law.");
    const Vec3 moment = spheres[0].torque() + spheres[1].torque() + math::cross(separation, spheres[1].force());
    near(moment, Vec3::zero(), 1.0e-10, "Tangential seam contact used a nonperiodic lever arm.");
    near(math::norm(model.contacts().front().tangentialForce_), 1.5, 1.0e-12, "Tangential dashpot coefficient was mass-scaled or otherwise rescaled.");
}

void historyAcrossSeams()
{
    const auto initial = pair({0.2, 0.2, 3.0}, {9.6, 7.6, 3.0});
    testShear reference(initial, lower, upper, parameters());
    testShear translated(initial, lower, upper, parameters());
    reference.setTimeStep(stepSize);
    translated.setTimeStep(stepSize);
    reference.velocity(0, {0.1, -0.1, 0.0});
    translated.velocity(0, {0.1, -0.1, 0.0});
    reference.initialize();
    translated.initialize();
    for (int step = 0; step < 10; ++step)
    {
        reference.step();
        translated.step();
    }
    require(math::norm(reference.contacts().front().tangentialSpring_) > 0.0, "History check did not build tangential deformation.");
    translated.translate({9.6, 7.7, 0.0});
    require(translated.initialized(), "Test translation unexpectedly invalidated the contact state.");
    reference.step();
    translated.step();
    require(reference.contacts().size() == 1 && translated.contacts().size() == 1, "Changing periodic images lost the contact.");
    const auto& first = reference.contacts().front();
    const auto& second = translated.contacts().front();
    near(second.tangentialSpring_, first.tangentialSpring_, 1.0e-11, "Tangential history was reset while crossing a periodic seam.");
    near(second.tangentialForce_, first.tangentialForce_, 1.0e-7, "Friction changed under a common periodic translation.");
    for (int index = 0; index < 2; ++index)
    {
        near(translated.spheres().host()[index].force(), reference.spheres().host()[index].force(), 1.0e-7, "Periodic translation changed particle force.");
        near(translated.spheres().host()[index].torque(), reference.spheres().host()[index].torque(), 1.0e-7, "Periodic translation changed particle torque.");
    }

    translated.place(0, {2.0, 2.0, 3.0});
    translated.place(1, {5.0, 5.0, 3.0});
    translated.step();
    require(translated.contacts().empty(), "Separated particles retained a contact.");
    translated.place(0, {0.2, 3.0, 3.0});
    translated.place(1, {9.3, 3.0, 3.0});
    translated.velocity(0, Vec3::zero());
    translated.velocity(1, Vec3::zero());
    translated.step();
    require(translated.contacts().size() == 1, "Recontact was not detected.");
    near(translated.contacts().front().tangentialSpring_, Vec3::zero(), 1.0e-12, "A reopened contact inherited stale tangential history.");
}

void coulombLimit()
{
    auto properties = parameters();
    properties.slidingDamping_ = 100.0;
    testShear model(pair({0.2, 3.0, 3.0}, {9.3, 3.0, 3.0}), lower, upper, properties);
    model.setTimeStep(stepSize);
    model.velocity(0, {0.0, 1000.0, 0.0});
    model.initialize();
    for (int step = 0; step < 3; ++step)
    {
        require(model.contacts().size() == 1, "Sliding check lost its overlapping pair.");
        const auto& contact = model.contacts().front();
        const Real limit = properties.friction_ * contact.normalForce_;
        near(math::norm(contact.tangentialForce_), limit, 1.0e-8, "Total elastic-plus-damping tangential force was not returned to the Coulomb disk.");
        near(math::dot(contact.tangentialForce_, contact.normal_), 0.0, 1.0e-9, "Tangential force has a normal component.");
        model.step();
    }
}

void wallMotion()
{
    const auto properties = parameters();
    const std::vector<initialSphere> initial{{2, {2.0, 2.0, 1.0}, 0.5, -1}, {3, {2.1, 2.0, 1.0}, 0.5, -1}, {7, {6.0, 2.0, 5.0}, 0.5, 1}, {8, {6.1, 2.0, 5.0}, 0.5, 1}, {11, {4.0, 5.0, 3.0}, 0.5, 0}};
    testShear model(initial, lower, upper, properties);
    model.setTimeStep(stepSize);
    model.initialize();
    require(model.contacts().empty(), "Fixed-fixed contacts should be skipped.");
    for (int step = 0; step < 50; ++step)
    {
        model.step();
    }
    require(model.contacts().empty(), "Fixed-fixed pairs were added while moving the walls.");
    const auto& spheres = model.spheres().host();
    for (int index = 0; index < static_cast<int>(initial.size()); ++index)
    {
        const auto& sphere = spheres[index];
        require(math::isFinite(sphere.position()) && math::isFinite(sphere.velocity()) && math::isFinite(sphere.angularVelocity()) && math::isFinite(sphere.orientation()) &&
                    math::isFinite(sphere.force()) && math::isFinite(sphere.torque()),
                "Short integration produced a nonfinite value.");
        const Vec3 velocity{initial[index].wall_ * properties.wallSpeed_, 0.0, 0.0};
        near(sphere.velocity(), velocity, 1.0e-12, "Wall velocity did not remain prescribed.");
        near(sphere.position(), model.wrapPosition(initial[index].position_ + model.time() * velocity), 1.0e-10, "A wall drifted away from its prescribed trajectory.");
        near(sphere.position().z, initial[index].position_.z, 1.0e-12, "Prescribed wall motion drifted in Z.");
        near(sphere.angularVelocity(), Vec3::zero(), 1.0e-12, "An unloaded wall developed spin.");
    }
}

void shortIntegration()
{
    auto properties = parameters();
    properties.normalDamping_ = 10.0;
    properties.slidingDamping_ = 3.0;
    testShear model(pair({0.2, 0.2, 3.0}, {9.6, 7.6, 3.0}), lower, upper, properties);
    model.setTimeStep(stepSize);
    model.velocity(0, {0.1, -0.1, 0.0});
    model.initialize();
    for (int step = 0; step < 500; ++step)
    {
        model.step();
        for (const auto& sphere : model.spheres().host())
        {
            require(math::isFinite(sphere.position()) && math::isFinite(sphere.velocity()) && math::isFinite(sphere.angularVelocity()) && math::isFinite(sphere.orientation()) &&
                        math::isFinite(sphere.force()) && math::isFinite(sphere.torque()),
                    "Short interacting integration produced a nonfinite value.");
            require(sphere.position().x >= lower.x && sphere.position().x < upper.x && sphere.position().y >= lower.y && sphere.position().y < upper.y,
                    "An integrated particle escaped its periodic cell.");
        }
    }
    near(model.time(), 500 * stepSize, 1.0e-14, "Short integration advanced the wrong represented time.");
}

void cellSearchAgainstBruteForce()
{
    const auto image = [](Vec3 displacement)
    {
        displacement.x -= 10.0 * std::round(displacement.x / 10.0);
        displacement.y -= 8.0 * std::round(displacement.y / 8.0);
        return displacement;
    };
    std::uint32_t state = 0x71a09f3du;
    const auto uniform = [&]()
    {
        state = 1664525u * state + 1013904223u;
        return (static_cast<Real>(state) + 0.5) / 4294967296.0;
    };
    std::vector<initialSphere> initial{{101, {0.1, 0.1, 2.0}, 0.5, 0}, {104, {9.85, 7.85, 2.0}, 0.5, 0}, {107, {1.01, 2.99, 3.0}, 0.5, 0}, {110, {1.99, 3.01, 3.0}, 0.5, 0}};
    for (int attempt = 0; initial.size() < 120 && attempt < 10000; ++attempt)
    {
        const Vec3 position{10.0 * uniform(), 8.0 * uniform(), 6.0 * uniform()};
        bool separated = true;
        for (const auto& sphere : initial)
        {
            if (math::normSquared(image(position - sphere.position_)) < 0.04)
            {
                separated = false;
                break;
            }
        }
        if (separated)
        {
            initial.push_back({101 + 3 * static_cast<int>(initial.size()), position, 0.5, 0});
        }
    }
    require(initial.size() == 120, "Could not construct the deterministic neighbor-search fixture.");

    using pairIndex = std::pair<int, int>;
    std::map<pairIndex, Real> expected;
    for (int first = 0; first < static_cast<int>(initial.size()); ++first)
    {
        for (int second = first + 1; second < static_cast<int>(initial.size()); ++second)
        {
            const Real distance = math::norm(image(initial[first].position_ - initial[second].position_));
            if (distance < 1.0)
            {
                expected.emplace(pairIndex{first, second}, 1.0 - distance);
            }
        }
    }
    require(expected.size() > 20, "Brute-force fixture did not exercise enough simultaneous contacts.");
    testShear model(initial, lower, upper, parameters());
    model.initialize();
    require(model.contacts().size() == expected.size(), "Periodic cell search and brute force returned different contact counts.");
    std::set<pairIndex> seen;
    for (const auto& contact : model.contacts())
    {
        const pairIndex key{contact.first_, contact.second_};
        require(contact.first_ < contact.second_ && seen.insert(key).second, "Periodic cell search returned a duplicate or reversed finite pair.");
        const auto found = expected.find(key);
        require(found != expected.end(), "Periodic cell search returned a pair absent from the brute-force reference.");
        near(contact.overlap_, found->second, 1.0e-12, "Cell-search overlap differs from the brute-force reference.");
        near(contact.normalForce_, parameters().normalStiffness_ * found->second, 1.0e-8, "Cell-search force differs from direct stiffness times overlap.");
        near(contact.tangentialForce_, Vec3::zero(), 1.0e-12, "A stationary brute-force fixture developed tangential force.");
    }
}

} // namespace

int main()
{
    try
    {
        periodicCoordinates();
        normalContact({0.2, 3.0, 3.0}, {9.3, 3.0, 3.0});
        normalContact({3.0, 0.2, 3.0}, {3.0, 7.3, 3.0});
        normalContact({0.2, 0.2, 3.0}, {9.6, 7.6, 3.0});
        normalContact({0.2, 0.2, 3.0}, {2.3, 2.3, 3.0}, {2.5, 2.5, 6.0});
        normalDashpot();
        wallContactCoefficients();
        tangentialBalance();
        historyAcrossSeams();
        coulombLimit();
        wallMotion();
        shortIntegration();
        cellSearchAgainstBruteForce();
        std::cout << "periodicShearTest: PASS (periodic contacts, balance, history, friction, wall motion).\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "periodicShearTest: FAIL: " << error.what() << '\n';
        return 1;
    }
}
