/**
 * @file main.cpp
 * @brief Runs the isolated CPU adaptation of the ON-DEM bumpy-wall shear benchmark.
 */
#include "benchmarkIO.h"
#include "data/vtuWriter.h"

#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
using namespace fundem;
using namespace fundem::validation;
namespace fs = std::filesystem;

/** Explicit run limits prevent accidentally starting the full 200-million-step benchmark. */
struct options {
    fs::path input_;
    fs::path output_;
    int steps_{0};
    int interval_{50000};
    math::Real timeStep_{7.2320e-5};
    shearParameters parameters_;
    bool writeOutput_{true};
    bool writeVTU_{false};
};

/** Parses one integer without accepting a truncated number or an overflow. */
int integerOption(const std::string& text, const char* option)
{
    std::size_t consumed = 0;
    const int result = std::stoi(text, &consumed);
    if (consumed != text.size() || result < 0)
        throw std::invalid_argument(std::string(option) + " requires a nonnegative integer.");
    return result;
}

/** Parses a finite coefficient without accepting trailing text or negative values. */
math::Real realOption(const std::string& text, const std::string& option)
{
    std::size_t consumed = 0;
    const math::Real result = std::stod(text, &consumed);
    if (consumed != text.size() || !math::isFinite(result) || result < 0.0)
        throw std::invalid_argument(option + " requires a finite nonnegative number.");
    return result;
}

/** Returns false only for the informational --help path. */
bool parseOptions(options& result, int argc, char** argv)
{
    const fs::path executableDirectory = fs::absolute(argv[0]).parent_path();
    result.input_ = executableDirectory / "ondemShearFlow_input" / "GranularShearFlow.data";
    result.output_ = executableDirectory / "ondemShearFlow_files";
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        if (argument == "--help")
        {
            std::cout << "ondemShearFlow [--input FILE] [--check | --steps N] [--output DIR]\n"
                         "               [--output-interval N] [--time-step DT] [--vtu] [--no-output]\n"
                         "               [--normal-damping CN] [--sliding-damping CT] [--sliding-stiffness KT]\n"
                         "Default: load/check the original 3480-particle configuration without advancing.\n"
                         "The published benchmark requires --steps 200000000 and interval 50000.\n"
                         "This validation-only CPU solver is not a Workbench project or GPU backend.\n";
            return false;
        }
        if (argument == "--check")
        {
            result.steps_ = 0;
            continue;
        }
        if (argument == "--vtu")
        {
            result.writeVTU_ = true;
            continue;
        }
        if (argument == "--no-output")
        {
            result.writeOutput_ = false;
            continue;
        }
        if (++index >= argc)
            throw std::invalid_argument("Missing value for " + argument);
        if (argument == "--input")
            result.input_ = fs::u8path(argv[index]);
        else if (argument == "--output")
            result.output_ = fs::u8path(argv[index]);
        else if (argument == "--steps")
            result.steps_ = integerOption(argv[index], "--steps");
        else if (argument == "--time-step")
        {
            result.timeStep_ = realOption(argv[index], argument);
            if (result.timeStep_ == 0.0)
                throw std::invalid_argument("--time-step requires a finite positive number.");
        }
        else if (argument == "--normal-damping")
            result.parameters_.normalDamping_ = realOption(argv[index], argument);
        else if (argument == "--sliding-damping")
            result.parameters_.slidingDamping_ = realOption(argv[index], argument);
        else if (argument == "--sliding-stiffness")
            result.parameters_.slidingStiffness_ = realOption(argv[index], argument);
        else if (argument == "--output-interval")
        {
            result.interval_ = integerOption(argv[index], "--output-interval");
            if (result.interval_ == 0)
                throw std::invalid_argument("--output-interval must be positive.");
        }
        else
            throw std::invalid_argument("Unknown option: " + argument);
    }
    return true;
}

/** Writes a viewable particle snapshot without routing through native nonperiodic contacts. */
void writeVTU(const periodicShear& simulation, const fs::path& directory, int frame)
{
    std::vector<math::Vec3> positions, velocities, angularVelocities, forces;
    std::vector<math::Real> radii;
    std::vector<int> ids, infiniteMass;
    const auto& particles = simulation.spheres().host();
    for (int index = 0; index < static_cast<int>(particles.size()); ++index)
    {
        const auto& value = particles[index];
        positions.push_back(value.position());
        velocities.push_back(value.velocity());
        angularVelocities.push_back(value.angularVelocity());
        forces.push_back(value.force());
        radii.push_back(value.radius());
        ids.push_back(simulation.sourceIds()[index]);
        infiniteMass.push_back(value.inverseMass() == 0.0 ? 1 : 0);
    }
    std::ostringstream name;
    name << "ShearFlow_" << std::setw(6) << std::setfill('0') << frame << ".vtu";
    vtuWriter writer((directory / name.str()).string());
    writer.setPoints(positions);
    writer.setVertexCells();
    writer.addPointData("velocity", velocities);
    writer.addPointData("angularVelocity", angularVelocities);
    writer.addPointData("force", forces);
    writer.addPointData("radius", radii);
    writer.addPointData("sourceId", ids);
    writer.addPointData("infiniteMass", infiniteMass);
    writer.write();
}
} // namespace

int main(int argc, char** argv)
{
    try
    {
        options settings;
        if (!parseOptions(settings, argc, argv))
            return 0;
        const auto initial = readInitialConfiguration(settings.input_);
        periodicShear simulation(initial.particles_, initial.minimum_, initial.maximum_, settings.parameters_);
        simulation.setTimeStep(settings.timeStep_);
        std::cout << "ON-DEM bumpy-wall shear adaptation: 3000 mobile + 480 prescribed wall spheres.\n"
                     "Periodic X/Y, rho=d=V=1, direct spring/dashpot contact coefficients.\n"
                     "Reference contact-law details and example timestamps need confirmation; see README.md.\n";
        const auto& parameters = simulation.parameters();
        std::cout << std::setprecision(17) << "Contact parameters: kn=" << parameters.normalStiffness_ << ", kt=" << parameters.slidingStiffness_
                  << ", cn=" << parameters.normalDamping_ << ", ct=" << parameters.slidingDamping_ << ", mu=" << parameters.friction_ << ".\n"
                  << "Time step: " << simulation.timeStep() << ".\n";
        simulation.initialize();
        if (settings.steps_ == 0)
        {
            std::cout << "CHECK PASS: configuration, stable IDs, periodic search and initial forces; no time integration.\n";
            return 0;
        }

        const fs::path output = fs::absolute(settings.output_);
        if (settings.writeOutput_)
        {
            if (fs::exists(output) && !fs::is_empty(output))
                throw std::runtime_error("Output directory is not empty; select a new directory to preserve existing results: " + output.string());
            fs::create_directories(output);
        }
        const auto started = std::chrono::steady_clock::now();
        int frame = 0;
        const auto writeFrame = [&]
        {
            if (!settings.writeOutput_)
                return;
            writeMercuryFrame(simulation, output, frame);
            if (settings.writeVTU_)
                writeVTU(simulation, output, frame);
            std::cout << "Step " << simulation.stepCount() << ", time " << std::setprecision(12) << simulation.time() << ", frame " << frame << '\n';
            ++frame;
        };
        writeFrame();
        for (int step = 0; step < settings.steps_; ++step)
        {
            simulation.step();
            if (simulation.stepCount() % settings.interval_ == 0)
                writeFrame();
            else if (simulation.stepCount() % 5000 == 0)
                std::cout << "Progress: step " << simulation.stepCount() << ", time " << std::setprecision(12) << simulation.time() << std::endl;
        }
        if (simulation.stepCount() % settings.interval_ != 0)
            writeFrame();
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        std::cout << "Completed " << simulation.stepCount() << " steps at time " << simulation.time() << " in " << seconds << " wall-clock seconds. This is not a steady-state validation result.\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << "ondemShearFlow: FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
