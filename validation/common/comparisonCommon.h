#pragma once
#include "data/myLSObject.h"
#include "material/material.h"
#include "material/materialPair.h"
#include "particle/LSParticle.h"
#include "particle/particle.h"
#include "solver/SphereDEM.h"
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace fundem::comparison
{
struct Check {
    int step{};
    std::string metric, target, errorKind;
    math::Real expected{}, tolerance{};
};
struct CaseSettings {
    const char* id;
    const char* family;
    const char* historicalOutcome;
    math::Real timeStep;
    int steps;
    math::Vec3 gravity;
    std::vector<Check> checks;
};
struct Body {
    std::string id;
    bool sphere;
    int index;
};
struct RunOptions {
    executionMode mode{executionMode::CPU};
    bool halfStep{};
    int stepLimit{-1};
    std::filesystem::path output;
    bool help{};
};
RunOptions parseOptions(int argc, char** argv, const char* id);
void printUsage(const char* id);
int execute(SphereDEM& simulation, const CaseSettings& settings, const std::vector<Body>& bodies, const RunOptions& options);

// Preserve a reference wall's sampling phase using only core geometry types.
std::unique_ptr<levelset::LSInfo> prepareWall(const levelset::LSInfo& analytic,
    math::Vec3 origin, math::Real spacing, int3 size, math::Real boundingRadius,
    bool reverse, const std::vector<math::Vec3>& nodes, const std::vector<int3>& faces);
}
