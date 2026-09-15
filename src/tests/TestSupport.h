/**
 * @file TestSupport.h
 * @brief Defines lightweight assertions and CUDA checks shared by unit tests.
 */
#pragma once

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

namespace fundem::test
{

inline void require(bool condition, const char* expression, const char* file, int line)
{
    if (!condition)
    {
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": requirement failed: " + expression);
    }
}

inline void requireCUDA(cudaError_t error, const char* operation, const char* file, int line)
{
    if (error != cudaSuccess)
    {
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": " + operation + " failed: " + cudaGetErrorString(error));
    }
}

} // namespace fundem::test

#define FUNDEM_TEST_REQUIRE(expression) ::fundem::test::require(static_cast<bool>(expression), #expression, __FILE__, __LINE__)
#define FUNDEM_TEST_REQUIRE_CUDA(operation) ::fundem::test::requireCUDA((operation), #operation, __FILE__, __LINE__)
