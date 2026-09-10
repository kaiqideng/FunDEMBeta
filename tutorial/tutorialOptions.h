/**
 * @file tutorialOptions.h
 * @brief Shared command-line options for backend-selectable tutorials.
 */
#pragma once

#include "solver/solver.h"

#include <stdexcept>
#include <string>

namespace fundem::tutorial
{

/** Runtime controls shared by tutorial2 and tutorial3. */
struct TutorialOptions {
    std::string outputDirectory_;            ///< Output directory resolved beside the executable when relative.
    int stepCount_{-1};                      ///< Optional step-count override; a negative value selects the case duration.
    int device_{0};                          ///< CUDA device used by GPU and Hybrid modes.
    executionMode mode_{executionMode::CPU}; ///< Selected CPU, GPU, or Hybrid execution mode.
    vtuFormat vtuFormat_{vtuFormat::binary}; ///< Binary output by default.
};

/** Parses the common backend, output, and step-count arguments. */
inline TutorialOptions parseTutorialOptions(int argc, char** argv, const char* defaultOutputDirectory)
{
    TutorialOptions result;
    result.outputDirectory_ = defaultOutputDirectory;
    for (int argumentIndex = 1; argumentIndex < argc; ++argumentIndex)
    {
        const std::string argument = argv[argumentIndex];
        const auto nextValue = [&]() -> const char*
        {
            if (++argumentIndex >= argc)
            {
                throw std::invalid_argument(argument + " requires a value.");
            }
            return argv[argumentIndex];
        };

        if (argument == "--ascii")
        {
            result.vtuFormat_ = vtuFormat::ascii;
        }
        else if (argument == "--steps")
        {
            result.stepCount_ = std::stoi(nextValue());
            if (result.stepCount_ < 0)
            {
                throw std::invalid_argument("--steps requires a non-negative integer.");
            }
        }
        else if (argument == "--device")
        {
            result.device_ = std::stoi(nextValue());
            if (result.device_ < 0)
            {
                throw std::invalid_argument("--device requires a non-negative integer.");
            }
        }
        else if (argument == "--output")
        {
            result.outputDirectory_ = nextValue();
            if (result.outputDirectory_.empty())
            {
                throw std::invalid_argument("--output requires a non-empty directory.");
            }
        }
        else if (argument == "--mode")
        {
            const std::string mode = nextValue();
            if (mode == "cpu")
            {
                result.mode_ = executionMode::CPU;
            }
            else if (mode == "gpu")
            {
                result.mode_ = executionMode::GPU;
            }
            else if (mode == "hybrid")
            {
                result.mode_ = executionMode::Hybrid;
            }
            else
            {
                throw std::invalid_argument("--mode must be cpu, gpu, or hybrid.");
            }
        }
        else
        {
            throw std::invalid_argument("Unknown option: " + argument);
        }
    }
    return result;
}

} // namespace fundem::tutorial
