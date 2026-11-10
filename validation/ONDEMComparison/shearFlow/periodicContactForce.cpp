#include "periodicContactForce.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace fundem::validation
{

namespace
{

/** Returns a representable cell count without reducing the interaction horizon. */
int dimensionCellCount(Real length, Real diameter)
{
    const Real count = std::max(1.0, std::floor(length / diameter));
    if (count > static_cast<Real>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("The periodic-contact cell count exceeds the supported integer index range.");
    }
    return static_cast<int>(count);
}

/** Wraps an integer coordinate into a periodic axis, including negative offsets. */
int wrappedCell(int coordinate, int count) noexcept
{
    const int remainder = coordinate % count;
    return remainder < 0 ? remainder + count : remainder;
}

} // namespace

periodicContactForce::periodicContactForce(const Vec3& minimum, const Vec3& maximum, Real radius, int particleCount, contactForceParameters parameters, periodicBoundaryFunctions boundary)
    : parameters_(parameters), domainMinimum_(minimum), domainLength_(maximum - minimum), boundary_(boundary)
{
    if ((boundary_.wrapPosition_ == nullptr) != (boundary_.minimumImage_ == nullptr))
        throw std::invalid_argument("The external periodic boundary must provide both wrapping and minimum-image callbacks.");
    for (int axis = 0; axis < 3; ++axis)
        periodicAxes_[axis] = boundary_.isPeriodic_ ? boundary_.isPeriodic_(boundary_.context_, axis) : axis < 2;
    if (std::none_of(periodicAxes_.begin(), periodicAxes_.end(), [](bool enabled) { return enabled; }))
        throw std::invalid_argument("At least one periodic axis must be enabled.");
    if (!math::isFinite(parameters.normalStiffness_) || parameters.normalStiffness_ <= 0.0 || !math::isFinite(parameters.slidingStiffness_) || parameters.slidingStiffness_ < 0.0 ||
        !math::isFinite(parameters.normalDamping_) || parameters.normalDamping_ < 0.0 || !math::isFinite(parameters.slidingDamping_) || parameters.slidingDamping_ < 0.0 ||
        !math::isFinite(parameters.friction_) || parameters.friction_ < 0.0)
    {
        throw std::invalid_argument("Periodic contact force requires positive normal stiffness and finite non-negative coefficients.");
    }
    if (particleCount <= 0)
    {
        throw std::invalid_argument("Periodic contact force requires a positive fixed particle count.");
    }
    if (!math::isFinite(minimum) || !math::isFinite(maximum) || !math::isFinite(domainLength_) || domainLength_.x <= 0.0 || domainLength_.y <= 0.0 || domainLength_.z <= 0.0)
    {
        throw std::invalid_argument("Periodic contact force requires finite, ordered domain bounds.");
    }

    const Real diameter = 2.0 * radius;
    if (!math::isFinite(radius) || radius <= 0.0 || !math::isFinite(diameter))
    {
        throw std::invalid_argument("Periodic contact force requires a finite positive common radius.");
    }
    for (int axis = 0; axis < 3; ++axis)
        if (periodicAxes_[axis] && domainLength_[axis] <= 2.0 * diameter)
            throw std::invalid_argument("Each enabled periodic dimension must be longer than two particle diameters.");
    cellCountX_ = dimensionCellCount(domainLength_.x, diameter);
    cellCountY_ = dimensionCellCount(domainLength_.y, diameter);
    cellCountZ_ = dimensionCellCount(domainLength_.z, diameter);
    const std::size_t cellCountXY = static_cast<std::size_t>(cellCountX_) * static_cast<std::size_t>(cellCountY_);
    if (cellCountXY > static_cast<std::size_t>(std::numeric_limits<int>::max()) / static_cast<std::size_t>(cellCountZ_))
    {
        throw std::invalid_argument("The periodic-contact cell list exceeds the supported integer index range.");
    }
    cellWidth_ = {domainLength_.x / cellCountX_, domainLength_.y / cellCountY_, domainLength_.z / cellCountZ_};
    cellHeads_.resize(cellCountXY * static_cast<std::size_t>(cellCountZ_), -1);
    next_.resize(static_cast<std::size_t>(particleCount), -1);
}

Vec3 periodicContactForce::wrapPosition(Vec3 position) const noexcept
{
    if (boundary_.wrapPosition_)
        return boundary_.wrapPosition_(boundary_.context_, position);
    for (int axis = 0; axis < 3; ++axis)
    {
        if (!periodicAxes_[axis])
            continue;
        Real coordinate = std::fmod(position[axis] - domainMinimum_[axis], domainLength_[axis]);
        if (coordinate < 0.0)
        {
            coordinate += domainLength_[axis];
        }
        // Rounding after addition can hit the upper endpoint of the half-open cell.
        position[axis] = domainMinimum_[axis] + coordinate;
        if (position[axis] >= domainMinimum_[axis] + domainLength_[axis])
        {
            position[axis] = domainMinimum_[axis];
        }
    }
    return position;
}

Vec3 periodicContactForce::minimumImage(Vec3 displacement) const noexcept
{
    if (boundary_.minimumImage_)
        return boundary_.minimumImage_(boundary_.context_, displacement);
    for (int axis = 0; axis < 3; ++axis)
        if (periodicAxes_[axis])
            displacement[axis] -= domainLength_[axis] * std::round(displacement[axis] / domainLength_[axis]);
    return displacement;
}

int periodicContactForce::cellIndex(const Vec3& position) const noexcept
{
    const Real xCoordinate = (position.x - domainMinimum_.x) / cellWidth_.x;
    const int x = xCoordinate <= 0.0 ? 0 : (xCoordinate >= cellCountX_ ? cellCountX_ - 1 : static_cast<int>(xCoordinate));
    const Real yCoordinate = (position.y - domainMinimum_.y) / cellWidth_.y;
    const int y = yCoordinate <= 0.0 ? 0 : (yCoordinate >= cellCountY_ ? cellCountY_ - 1 : static_cast<int>(yCoordinate));
    const Real zCoordinate = (position.z - domainMinimum_.z) / cellWidth_.z;
    const int z = zCoordinate <= 0.0 ? 0 : (zCoordinate >= cellCountZ_ ? cellCountZ_ - 1 : static_cast<int>(zCoordinate));
    return x + cellCountX_ * (y + cellCountY_ * z);
}

void periodicContactForce::buildCells(const particleContainer::host_container_type& particles)
{
    std::fill(cellHeads_.begin(), cellHeads_.end(), -1);
    for (int index = static_cast<int>(particles.size()) - 1; index >= 0; --index)
    {
        const int cell = cellIndex(particles[index].position());
        next_[index] = cellHeads_[cell];
        cellHeads_[cell] = index;
    }
}

std::uint64_t periodicContactForce::pairKey(int first, int second) noexcept
{
    const auto lower = static_cast<std::uint32_t>(std::min(first, second));
    const auto upper = static_cast<std::uint32_t>(std::max(first, second));
    return (static_cast<std::uint64_t>(lower) << 32) | upper;
}

void periodicContactForce::evaluatePair(particleContainer::host_container_type& particles, int first, int second, Real historyTimeStep)
{
    if (particles[first].inverseMass() == 0.0 && particles[second].inverseMass() == 0.0)
    {
        return;
    }
    if (particles[first].inverseMass() == 0.0)
    {
        std::swap(first, second);
    }
    particle& firstSphere = particles[first];
    particle& secondSphere = particles[second];
    const Vec3 separation = minimumImage(firstSphere.position() - secondSphere.position());
    const Real radiusSum = firstSphere.radius() + secondSphere.radius();
    const Real distanceSquared = math::normSquared(separation);
    if (distanceSquared >= radiusSum * radiusSum)
    {
        return;
    }
    const Real distance = std::sqrt(distanceSquared);
    if (distance <= math::defaultTolerance)
    {
        throw std::runtime_error("Coincident sphere centers have no defined shear contact normal.");
    }
    const Vec3 normal = separation / distance;
    const Real overlap = radiusSum - distance;
    const Vec3 firstArm = -(firstSphere.radius() - 0.5 * overlap) * normal;
    const Vec3 secondArm = (secondSphere.radius() - 0.5 * overlap) * normal;
    const Vec3 relativeVelocity = firstSphere.velocity() + math::cross(firstSphere.angularVelocity(), firstArm) - secondSphere.velocity() - math::cross(secondSphere.angularVelocity(), secondArm);
    const Real normalVelocity = math::dot(relativeVelocity, normal);
    const Vec3 tangentialVelocity = relativeVelocity - normalVelocity * normal;
    const Real signedNormalForce = parameters_.normalStiffness_ * overlap - parameters_.normalDamping_ * normalVelocity;
    const Real normalForce = parameters_.normalRelease_ == normalReleasePolicy::retainUntilSeparation ? signedNormalForce : std::max(0.0, signedNormalForce);

    contactHistory& history = histories_[pairKey(first, second)];
    history.generation_ = generation_;
    const Real previousLength = math::norm(history.spring_);
    Vec3 spring = history.spring_ - math::dot(history.spring_, normal) * normal;
    const Real projectedLength = math::norm(spring);
    // Use a relative degeneracy test so small but valid histories still accumulate.
    if (projectedLength > math::defaultTolerance * previousLength)
    {
        spring *= previousLength / projectedLength;
    }
    else
    {
        spring = Vec3::zero();
    }
    spring += historyTimeStep * tangentialVelocity;
    Vec3 tangentialForce = -parameters_.slidingStiffness_ * spring - parameters_.slidingDamping_ * tangentialVelocity;
    const Real tangentialMagnitude = math::norm(tangentialForce);
    const Real frictionLimit = parameters_.friction_ * std::abs(normalForce);
    if (tangentialMagnitude > frictionLimit)
    {
        tangentialForce *= frictionLimit / tangentialMagnitude;
        spring = parameters_.slidingStiffness_ > 0.0 ? -(tangentialForce + parameters_.slidingDamping_ * tangentialVelocity) / parameters_.slidingStiffness_ : Vec3::zero();
    }
    if (parameters_.slidingStiffness_ == 0.0)
    {
        spring = Vec3::zero();
    }
    history.spring_ = spring;

    const Vec3 force = normalForce * normal + tangentialForce;
    firstSphere.addForce(force);
    secondSphere.addForce(-force);
    firstSphere.addTorque(math::cross(firstArm, tangentialForce));
    secondSphere.addTorque(math::cross(secondArm, -tangentialForce));
    contacts_.push_back({first, second, wrapPosition(firstSphere.position() + firstArm), normal, tangentialForce, spring, overlap, normalForce});
}

void periodicContactForce::apply(particleContainer::host_container_type& particles, Real historyTimeStep)
{
    if (particles.size() != next_.size())
    {
        throw std::invalid_argument("The periodic contact-force population must retain its original count and index ordering.");
    }
    contacts_.clear();
    buildCells(particles);
    ++generation_;
    for (int first = 0; first < static_cast<int>(particles.size()); ++first)
    {
        const int cell = cellIndex(particles[first].position());
        const int cellX = cell % cellCountX_;
        const int cellY = (cell / cellCountX_) % cellCountY_;
        const int cellZ = cell / (cellCountX_ * cellCountY_);
        std::array<int, 27> visited{};
        int visitedCount = 0;
        for (int dz = -1; dz <= 1; ++dz)
        {
            const int z = periodicAxes_[2] ? wrappedCell(cellZ + dz, cellCountZ_) : cellZ + dz;
            if (z < 0 || z >= cellCountZ_)
            {
                continue;
            }
            for (int dy = -1; dy <= 1; ++dy)
            {
                const int y = periodicAxes_[1] ? wrappedCell(cellY + dy, cellCountY_) : cellY + dy;
                if (y < 0 || y >= cellCountY_)
                    continue;
                for (int dx = -1; dx <= 1; ++dx)
                {
                    const int x = periodicAxes_[0] ? wrappedCell(cellX + dx, cellCountX_) : cellX + dx;
                    if (x < 0 || x >= cellCountX_)
                        continue;
                    const int neighborCell = x + cellCountX_ * (y + cellCountY_ * z);
                    if (std::find(visited.begin(), visited.begin() + visitedCount, neighborCell) != visited.begin() + visitedCount)
                    {
                        continue;
                    }
                    visited[visitedCount++] = neighborCell;
                    for (int second = cellHeads_[neighborCell]; second >= 0; second = next_[second])
                    {
                        if (second > first)
                        {
                            evaluatePair(particles, first, second, historyTimeStep);
                        }
                    }
                }
            }
        }
    }
    for (auto iterator = histories_.begin(); iterator != histories_.end();)
    {
        if (iterator->second.generation_ != generation_)
        {
            iterator = histories_.erase(iterator);
        }
        else
        {
            ++iterator;
        }
    }
}

} // namespace fundem::validation
