/**
 * @file validationCommon.h
 * @brief Small output and comparison utilities shared by validation cases.
 */
#pragma once

#include "math/Vector3.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fundem::validation
{

using math::Real;

/** Returns a scale-aware scalar relative error. */
inline Real relativeError(Real value, Real reference) noexcept { return std::abs(value - reference) / std::max({Real{1.0}, std::abs(value), std::abs(reference)}); }

/** Returns a scale-aware Euclidean vector error. */
inline Real relativeError(const math::Vec3& value, const math::Vec3& reference) noexcept { return math::norm(value - reference) / std::max({Real{1.0}, math::norm(value), math::norm(reference)}); }

/** Throws a validation failure with a concise diagnostic. */
inline void require(bool condition, const std::string& message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

/** Resolves the case output directory beside the executable or from --output. */
inline std::filesystem::path outputDirectory(int argc, char** argv, const char* caseName)
{
    std::filesystem::path result;
    for (int argumentIndex = 1; argumentIndex < argc; ++argumentIndex)
    {
        const std::string argument = argv[argumentIndex];
        if (argument != "--output" || ++argumentIndex >= argc)
        {
            throw std::invalid_argument("Only --output DIR is supported.");
        }
        result = argv[argumentIndex];
    }

    std::error_code error;
    std::filesystem::path executable = std::filesystem::weakly_canonical(argv[0], error);
    if (error)
    {
        error.clear();
        executable = std::filesystem::absolute(argv[0], error);
    }
    const std::filesystem::path executableDirectory = error ? std::filesystem::current_path() : executable.parent_path();
    if (result.empty())
    {
        result = executableDirectory / (std::string(caseName) + "_files");
    }
    else if (result.is_relative())
    {
        result = executableDirectory / result;
    }
    std::filesystem::create_directories(result);
    return result;
}

/** Opens one scientific-notation DAT file and fails loudly on I/O errors. */
inline std::ofstream openDataFile(const std::filesystem::path& directory, const char* fileName)
{
    std::ofstream stream(directory / fileName, std::ios::trunc);
    if (!stream)
    {
        throw std::runtime_error("Cannot open validation output file: " + (directory / fileName).string());
    }
    stream << std::scientific << std::setprecision(16);
    return stream;
}

/** Prints a uniform one-line PASS summary. */
inline void printPass(const char* caseName, Real maximumError, Real tolerance)
{
    std::cout << std::scientific << std::setprecision(6) << caseName << ": PASS, maximum error=" << maximumError << ", tolerance=" << tolerance << '\n';
}

} // namespace fundem::validation
