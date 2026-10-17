#include "benchmarkIO.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace fundem::validation
{
namespace
{

/** Reads one nonempty record, permitting blank lines and whole-line comments. */
std::string nextRecord(std::istream& stream)
{
    std::string line;
    while (std::getline(stream, line))
    {
        const auto first = line.find_first_not_of(" \t\r");
        if (first != std::string::npos && line[first] != '#')
            return line;
    }
    return {};
}

/** Uses a locale-independent numeric format and rejects trailing record fields. */
void requireEnd(std::istringstream& stream, const char* record)
{
    stream >> std::ws;
    if (!stream.eof())
        throw std::runtime_error(std::string("Unexpected fields in the initial ") + record + " record.");
}

/** Compares source decimal values with the library's scale-aware tolerance. */
bool close(Real first, Real second) { return std::abs(first - second) <= math::defaultTolerance * std::max({Real{1.0}, std::abs(first), std::abs(second)}); }

/** Writes a vector as three space-separated values. */
void vectorFields(std::ostream& stream, const Vec3& value) { stream << value.x << ' ' << value.y << ' ' << value.z; }

/** Exports a normalized Hamilton orientation as XYZ roll/pitch/yaw angles. */
Vec3 eulerAngles(math::Quaternion value)
{
    if (!math::tryNormalize(value))
        throw std::runtime_error("A particle has an invalid orientation during Mercury export.");
    return {std::atan2(2.0 * (value.w * value.x + value.y * value.z), 1.0 - 2.0 * (value.x * value.x + value.y * value.y)),
            std::asin(std::clamp(2.0 * (value.w * value.y - value.z * value.x), Real{-1.0}, Real{1.0})),
            std::atan2(2.0 * (value.w * value.z + value.x * value.y), 1.0 - 2.0 * (value.y * value.y + value.z * value.z))};
}

/** Writes one directed force record in the benchmark reference column order. */
void contactRecord(std::ostream& stream, const periodicShear& simulation, const shearContact& contact, bool reverse)
{
    const Real sign = reverse ? -1.0 : 1.0;
    const Real tangentialMagnitude = math::norm(contact.tangentialForce_);
    const Vec3 tangent = tangentialMagnitude > 0.0 ? contact.tangentialForce_ * (sign / tangentialMagnitude) : Vec3::zero();
    stream << simulation.time() << ' ' << (reverse ? contact.second_ : contact.first_) << ' ' << (reverse ? contact.first_ : contact.second_) << ' ';
    vectorFields(stream, simulation.wrapPosition(contact.point_));
    stream << ' ' << contact.overlap_ << ' ' << -math::norm(contact.tangentialSpring_) << ' ' << contact.normalForce_ << ' ' << tangentialMagnitude << ' ';
    vectorFields(stream, contact.normal_ * sign);
    stream << ' ';
    vectorFields(stream, tangent);
    stream << '\n';
}

} // namespace

initialConfiguration readInitialConfiguration(const std::filesystem::path& path, bool requireBenchmark)
{
    std::ifstream input(path);
    input.imbue(std::locale::classic());
    if (!input)
        throw std::runtime_error("Cannot open initial configuration: " + path.string());
    initialConfiguration result;
    std::istringstream bounds(nextRecord(input));
    bounds.imbue(std::locale::classic());
    const char* labels[] = {"xMin", "xMax", "yMin", "yMax", "zMin", "zMax"};
    Real* coordinates[] = {&result.minimum_.x, &result.maximum_.x, &result.minimum_.y, &result.maximum_.y, &result.minimum_.z, &result.maximum_.z};
    for (int index = 0; index < 6; ++index)
    {
        std::string label;
        if (!(bounds >> label >> *coordinates[index]) || label != labels[index] || !math::isFinite(*coordinates[index]))
            throw std::runtime_error("Invalid labeled bounds in initial configuration.");
    }
    requireEnd(bounds, "bounds");
    const Vec3 extent = result.maximum_ - result.minimum_;
    if (!math::isFinite(extent) || extent.x <= 0.0 || extent.y <= 0.0 || extent.z <= 0.0)
        throw std::runtime_error("Initial configuration bounds must be finite and ordered.");

    std::istringstream countRecord(nextRecord(input));
    std::string countLabel;
    int count = 0;
    if (!(countRecord >> countLabel >> count) || countLabel != "Particles" || count <= 0)
        throw std::runtime_error("Invalid initial configuration particle count.");
    requireEnd(countRecord, "particle count");
    if (requireBenchmark && (count != 3480 || !close(result.minimum_.x, 0.0) || !close(result.minimum_.y, 0.0) || !close(result.minimum_.z, 0.0) || !close(result.maximum_.x, 20.0) ||
                             !close(result.maximum_.y, 10.3923048454133) || !close(result.maximum_.z, 22.0)))
        throw std::runtime_error("The benchmark requires its original 3480-particle domain.");

    std::istringstream columns(nextRecord(input));
    const char* expectedColumns[] = {"id", "x-position", "y-position", "z-position", "radius"};
    for (const char* expected : expectedColumns)
    {
        std::string column;
        if (!(columns >> column) || column != expected)
            throw std::runtime_error("Unexpected initial configuration column labels.");
    }
    requireEnd(columns, "column labels");

    std::unordered_set<int> ids;
    for (int index = 0; index < count; ++index)
    {
        std::istringstream record(nextRecord(input));
        record.imbue(std::locale::classic());
        initialSphere sphere;
        if (!(record >> sphere.sourceId_ >> sphere.position_.x >> sphere.position_.y >> sphere.position_.z >> sphere.radius_) || sphere.sourceId_ <= 0 || !math::isFinite(sphere.position_) ||
            !math::isFinite(sphere.radius_) || sphere.radius_ <= 0.0 || !ids.insert(sphere.sourceId_).second)
            throw std::runtime_error("Invalid or duplicate initial particle at row " + std::to_string(index + 1) + '.');
        requireEnd(record, "particle");
        for (int axis = 0; axis < 3; ++axis)
        {
            if (sphere.position_[axis] < result.minimum_[axis] || sphere.position_[axis] > result.maximum_[axis])
                throw std::runtime_error("An initial particle center lies outside the source domain.");
        }
        if (requireBenchmark)
        {
            sphere.wall_ = index < 240 ? -1 : (index < 480 ? 1 : 0);
            if (sphere.sourceId_ != index + 1 || !close(sphere.radius_, 0.5) || (sphere.wall_ != 0 && !close(sphere.position_.z, sphere.wall_ < 0 ? 0.5 : 21.5)))
                throw std::runtime_error("Original IDs, unit diameters or wall ordering do not match the benchmark.");
        }
        else
        {
            sphere.wall_ = close(sphere.position_.z, result.minimum_.z + sphere.radius_) ? -1 : (close(sphere.position_.z, result.maximum_.z - sphere.radius_) ? 1 : 0);
        }
        result.particles_.push_back(sphere);
    }
    if (!nextRecord(input).empty())
        throw std::runtime_error("The initial configuration contains more particles than its declared count.");
    return result;
}

void writeMercuryFrame(const periodicShear& simulation, const std::filesystem::path& directory, int frame)
{
    if (frame < 0)
        throw std::invalid_argument("A Mercury output frame index cannot be negative.");
    std::ostringstream suffix;
    suffix << std::setw(4) << std::setfill('0') << frame;
    const auto dataPath = directory / ("ShearFlow.data." + suffix.str());
    const auto forcePath = directory / ("ShearFlow.fstat." + suffix.str());
    if (std::filesystem::exists(dataPath) || std::filesystem::exists(forcePath))
        throw std::runtime_error("Mercury export refuses to overwrite an existing frame.");
    std::filesystem::create_directories(directory);
    std::ofstream data(dataPath), forces(forcePath);
    if (!data || !forces)
        throw std::runtime_error("Cannot open Mercury frame output files.");
    data.imbue(std::locale::classic());
    forces.imbue(std::locale::classic());
    data << std::scientific << std::setprecision(16);
    forces << std::scientific << std::setprecision(16);
    const auto& particles = simulation.spheres().host();
    const Vec3& minimum = simulation.minimumBoundary();
    const Vec3& maximum = simulation.maximumBoundary();
    data << particles.size() << ' ' << simulation.time() << ' ';
    vectorFields(data, minimum);
    data << ' ';
    vectorFields(data, maximum);
    data << '\n';
    Real smallestRadius = particles.front().radius(), largestRadius = smallestRadius;
    for (const particle& sphere : particles)
    {
        vectorFields(data, simulation.wrapPosition(sphere.position()));
        data << ' ';
        vectorFields(data, sphere.velocity());
        data << ' ' << sphere.radius() << ' ';
        vectorFields(data, eulerAngles(sphere.orientation()));
        data << ' ';
        vectorFields(data, sphere.angularVelocity());
        data << " 0\n";
        smallestRadius = std::min(smallestRadius, sphere.radius());
        largestRadius = std::max(largestRadius, sphere.radius());
    }
    forces << "# " << simulation.time() << " 1\n# ";
    vectorFields(forces, minimum);
    forces << ' ';
    vectorFields(forces, maximum);
    forces << "\n# " << smallestRadius << ' ' << largestRadius << " 0 0 0 0\n";
    for (const shearContact& contact : simulation.contacts())
    {
        contactRecord(forces, simulation, contact, false);
        if (particles[contact.second_].inverseMass() > 0.0)
            contactRecord(forces, simulation, contact, true);
    }
    data.close();
    forces.close();
    if (!data || !forces)
        throw std::runtime_error("Writing a Mercury frame failed.");
}

} // namespace fundem::validation
