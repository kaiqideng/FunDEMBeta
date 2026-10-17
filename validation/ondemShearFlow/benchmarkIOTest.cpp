/** @file benchmarkIOTest.cpp @brief Bounded coordinate-import and Mercury-frame regressions. */
#include "benchmarkIO.h"
#include "../validationCommon.h"

#include <array>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
using namespace fundem;
using namespace fundem::validation;
namespace fs = std::filesystem;

/** Removes only the uniquely created, test-owned temporary directory. */
struct temporaryDirectory {
    fs::path path_;
    temporaryDirectory()
    {
        const auto token = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        path_ = fs::temp_directory_path() / ("fundem-shear-io-" + std::to_string(token));
        require(fs::create_directory(path_), "Could not create the test-owned temporary directory.");
    }
    ~temporaryDirectory()
    {
        std::error_code ignored;
        fs::remove_all(path_, ignored);
    }
};

/** Creates one small source-file fixture, never an application project or result. */
void fixture(const fs::path& path, const std::string& rows, int count = 4)
{
    std::ofstream stream(path);
    stream << "xMin 0 xMax 10 yMin 0 yMax 10 zMin 0 zMax 10\nParticles " << count << "\nid x-position y-position z-position radius\n" << rows;
    require(static_cast<bool>(stream), "Could not write the small input fixture.");
}

/** Requires malformed source input to fail before constructing a solver. */
void rejected(const fs::path& path)
{
    bool threw = false;
    try
    {
        (void)readInitialConfiguration(path, false);
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    require(threw, "Malformed initial coordinates were accepted.");
}

/** Parses an exact numeric row without silently accepting missing columns. */
std::vector<Real> values(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<Real> result;
    Real value;
    while (stream >> value)
        result.push_back(value);
    require(stream.eof(), "Non-numeric output field.");
    return result;
}
} // namespace

int main()
{
    try
    {
        temporaryDirectory temporary;
        const auto source = temporary.path_ / "small.data";
        fixture(source, "1 2 2 0.5 0.5\n2 2 2 1.4 0.5\n3 2.9 2 1.4 0.5\n4 8 8 9.5 0.5\n");
        const auto initial = readInitialConfiguration(source, false);
        require(initial.particles_.size() == 4 && initial.particles_[0].wall_ == -1 && initial.particles_[1].wall_ == 0 && initial.particles_[3].wall_ == 1,
                "Import changed source order or wall roles.");
        bool strictRejected = false;
        try
        {
            (void)readInitialConfiguration(source);
        }
        catch (const std::exception&)
        {
            strictRejected = true;
        }
        require(strictRejected, "Benchmark mode accepted a nonbenchmark particle count.");

        periodicShear simulation(initial.particles_, initial.minimum_, initial.maximum_);
        simulation.initialize();
        require(simulation.contacts().size() == 2, "Small export fixture should contain one grain/grain and one grain/wall contact.");
        const auto output = temporary.path_ / "frames";
        writeMercuryFrame(simulation, output, 0);
        std::ifstream data(output / "ShearFlow.data.0000");
        std::string line;
        std::getline(data, line);
        const auto header = values(line);
        require(header.size() == 8 && header[0] == 4 && header[1] == 0.0 && header[5] == 10.0, "Invalid Mercury data header.");
        int particleRows = 0;
        while (std::getline(data, line))
        {
            const auto row = values(line);
            require(row.size() == 14 && row[6] == 0.5, "Invalid Mercury particle fields.");
            ++particleRows;
        }
        require(particleRows == 4, "Wall particles must remain in the data file.");

        std::ifstream forces(output / "ShearFlow.fstat.0000");
        for (int index = 0; index < 3; ++index)
        {
            std::getline(forces, line);
            require(!line.empty() && line[0] == '#', "Mercury force metadata is missing.");
        }
        int forceRows = 0, wallRows = 0, pairRows = 0;
        while (std::getline(forces, line))
        {
            const auto row = values(line);
            require(row.size() == 16, "Invalid Mercury force-column count.");
            const int first = static_cast<int>(row[1]), second = static_cast<int>(row[2]);
            require(first >= 0 && first < 4 && second >= 0 && second < 4 && first != second, "Force endpoint is not a zero-based data-row index.");
            require(simulation.spheres().host()[first].inverseMass() > 0.0, "A wall must not be the first directed force endpoint.");
            require(row[6] > 0.0 && row[7] <= 0.0 && row[8] >= 0.0 && row[9] >= 0.0, "Force or spring signs differ from the export convention.");
            const Vec3 normal{row[10], row[11], row[12]}, tangent{row[13], row[14], row[15]};
            require(relativeError(math::norm(normal), 1.0) < 1.0e-12, "Exported normal is not a unit vector.");
            const Vec3 reconstructed = normal * row[8] + tangent * row[9];
            bool matched = false;
            for (const auto& contact : simulation.contacts())
            {
                const bool forward = contact.first_ == first && contact.second_ == second;
                const bool reverse = contact.first_ == second && contact.second_ == first;
                if (forward || reverse)
                {
                    const Vec3 expected = (contact.normal_ * contact.normalForce_ + contact.tangentialForce_) * (forward ? 1.0 : -1.0);
                    require(relativeError(reconstructed, expected) < 1.0e-12, "Force direction or reversed contact sign was lost in export.");
                    matched = true;
                }
            }
            require(matched, "Export introduced a contact absent from the solver.");
            simulation.spheres().host()[second].inverseMass() == 0.0 ? ++wallRows : ++pairRows;
            ++forceRows;
        }
        require(forceRows == 3 && wallRows == 1 && pairRows == 2, "Directed force multiplicity differs from the reference format.");
        bool overwriteRejected = false;
        try
        {
            writeMercuryFrame(simulation, output, 0);
        }
        catch (const std::exception&)
        {
            overwriteRejected = true;
        }
        require(overwriteRejected, "Mercury export overwrote a prior result frame.");

        fixture(source, "1 1 1 1 0.5\n1 2 2 2 0.5\n", 2);
        rejected(source);
        fixture(source, "1 1 1 1 0.5\n", 2);
        rejected(source);
        fixture(source, "1 1 1 1 0.5 extra\n", 1);
        rejected(source);
        fixture(source, "1 1 1 1 -0.5\n", 1);
        rejected(source);
        fixture(source, "1 11 1 1 0.5\n", 1);
        rejected(source);
        std::cout << "ON-DEM coordinate and Mercury IO: PASS; no time integration.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
