/** Core-only output/check helpers. Each executable constructs its own scene in C++. */
#include "comparisonCommon.h"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

namespace fundem::comparison
{
namespace
{
namespace fs = std::filesystem;
using Vec3 = math::Vec3;
using Real = math::Real;
void need(bool ok, const std::string& error)
{
    if (!ok) throw std::runtime_error(error);
}
Real squared(const Vec3& value) { return math::dot(value, value); }
struct RuntimeBody {
    std::string id;
    bool sphere;
    int index;
    Vec3 initialPosition, initialVelocity;
};
const rigidBody& state(const SphereDEM& solver, const RuntimeBody& value)
{
    if (value.sphere) return solver.spheres().host().at(static_cast<std::size_t>(value.index));
    return solver.LSParticles().host().at(static_cast<std::size_t>(value.index));
}
struct Snapshot {
    Real force{}, overlap{}, elastic{}, translational{}, momentumX{};
};
Real measure(const Check& check, const SphereDEM& solver, const std::vector<RuntimeBody>& bodies, const Snapshot& frame)
{
    if (check.metric == "normal_force")
        return frame.force;
    if (check.metric == "overlap")
        return frame.overlap;
    if (check.metric == "translational_plus_contact_energy")
        return frame.translational + frame.elastic;
    if (check.metric == "momentum_x")
        return frame.momentumX;
    std::vector<const RuntimeBody*> spheres;
    for (const auto& value : bodies)
        if (value.sphere)
            spheres.push_back(&value);
    need(!spheres.empty(), "Metric requires a sphere.");
    if (check.metric == "restitution")
    {
        if (spheres.size() == 2)
        {
            const Vec3 normal = math::normalizedOrZero(spheres[1]->initialPosition - spheres[0]->initialPosition);
            const Real incoming = math::dot(spheres[0]->initialVelocity - spheres[1]->initialVelocity, normal);
            need(incoming > 0.0, "Restitution requires a positive approach speed.");
            return math::dot(state(solver, *spheres[1]).velocity() - state(solver, *spheres[0]).velocity(), normal) / incoming;
        }
        need(spheres.size() == 1, "Restitution fixture requires one sphere/wall or two spheres.");
        const Vec3 initial = spheres[0]->initialVelocity;
        need(squared(initial) > 0.0, "Restitution fixture has zero approach speed.");
        return -math::dot(state(solver, *spheres[0]).velocity(), initial) / squared(initial);
    }
    const RuntimeBody* target = nullptr;
    for (const auto& value : bodies)
        if (value.id == check.target)
            target = &value;
    need(target != nullptr, "Unknown check body: " + check.target);
    const auto& value = state(solver, *target);
    if (check.metric == "x")
        return value.position().x;
    if (check.metric == "z")
        return value.position().z;
    if (check.metric == "vx")
        return value.velocity().x;
    if (check.metric == "vy")
        return value.velocity().y;
    if (check.metric == "vz")
        return value.velocity().z;
    if (check.metric == "wy")
        return value.angularVelocity().y;
    if (check.metric == "wz")
        return value.angularVelocity().z;
    throw std::runtime_error("Unknown metric: " + check.metric);
}
void validateChecks(const CaseSettings& settings, const std::map<std::string, bool>& bodyIds, const std::vector<RuntimeBody>& bodies)
{
    for (const auto& check : settings.checks)
    {
        need(check.step >= 0 && check.step <= settings.steps, "Checkpoint step is outside the case duration.");
        need(math::isFinite(check.expected) && math::isFinite(check.tolerance) && check.tolerance >= 0.0,
             "Checkpoint expected value and nonnegative tolerance must be finite.");
        need(check.errorKind == "absolute" || check.errorKind == "relative", "Unknown checkpoint error kind: " + check.errorKind);
        need(check.errorKind != "relative" || check.expected != 0.0, "A relative checkpoint requires a nonzero expected value.");
        const bool aggregate = check.metric == "normal_force" || check.metric == "overlap" ||
                               check.metric == "translational_plus_contact_energy" || check.metric == "momentum_x";
        const bool body = check.metric == "x" || check.metric == "z" || check.metric == "vx" || check.metric == "vy" ||
                          check.metric == "vz" || check.metric == "wy" || check.metric == "wz";
        need(aggregate || body || check.metric == "restitution", "Unknown checkpoint metric: " + check.metric);
        need(!body || bodyIds.count(check.target) != 0, "Unknown checkpoint body: " + check.target);
        if (check.metric == "restitution")
        {
            std::vector<const RuntimeBody*> spheres;
            for (const auto& value : bodies)
                if (value.sphere)
                    spheres.push_back(&value);
            need(spheres.size() == 1 || spheres.size() == 2, "Restitution fixture requires one sphere/wall or two spheres.");
            if (spheres.size() == 2)
            {
                const Vec3 normal = math::normalizedOrZero(spheres[1]->initialPosition - spheres[0]->initialPosition);
                const Real incoming = math::dot(spheres[0]->initialVelocity - spheres[1]->initialVelocity, normal);
                need(math::isFinite(incoming) && incoming > 0.0, "Restitution requires a finite positive approach speed.");
            }
            else
            {
                const Real speedSquared = squared(spheres[0]->initialVelocity);
                need(math::isFinite(speedSquared) && speedSquared > 0.0, "Restitution fixture has invalid/zero approach speed.");
            }
        }
    }
}
void validateBodyState(const rigidBody& body, const std::string& id)
{
    need(math::isFinite(body.position()) && math::isFinite(body.orientation()) && math::isFinite(body.velocity()) &&
         math::isFinite(body.angularVelocity()) && math::isFinite(body.force()) && math::isFinite(body.torque()) &&
         math::isFinite(body.inverseMass()) && body.inverseMass() >= 0.0, "Invalid/non-finite body state: " + id);
}
void writeVector(std::ostream& output, const Vec3& value) { output << ',' << value.x << ',' << value.y << ',' << value.z; }

} // namespace

std::unique_ptr<levelset::LSInfo> prepareWall(const levelset::LSInfo& analytic,
    math::Vec3 origin, math::Real spacing, int3 size, math::Real radius,
    bool reverse, const std::vector<math::Vec3>& nodes, const std::vector<int3>& faces)
{
    need(math::isFinite(origin) && math::isFinite(spacing) && spacing > 0 &&
         size.x >= 2 && size.y >= 2 && size.z >= 2 && !nodes.empty(), "Invalid prepared wall.");
    std::uint64_t count=1;
    for(int dimension : {size.x,size.y,size.z})
    {
        need(std::uint64_t(dimension) <= 50000000/count, "Prepared wall exceeds grid budget.");
        count*=std::uint64_t(dimension);
    }
    std::vector<math::Real> distances;
    distances.reserve(static_cast<std::size_t>(count));
    for(int z=0;z<size.z;++z) for(int y=0;y<size.y;++y) for(int x=0;x<size.x;++x)
    {
        const auto value=analytic.signedDistance(origin+math::Vec3{math::Real(x),math::Real(y),math::Real(z)}*spacing);
        distances.push_back(reverse ? -value : value);
    }
    return levelset::LSInfo::fromPreparedData(origin,spacing,size,distances,nodes,faces,radius,0.0,math::Mat3::zero());
}

void printUsage(const char* id)
{
    std::cout << "Usage: " << id << " [--mode cpu|gpu|hybrid] [--half-step] [--steps N] [--output NEW_DIRECTORY]\n"
                 "Scene and reference criteria are compiled from this case's C++ source.\n";
}

RunOptions parseOptions(int argc,char** argv,const char*)
{
    RunOptions options;
    for(int i=1;i<argc;++i)
    {
        const std::string argument=argv[i];
        const auto next=[&]() -> std::string {
            need(i+1<argc,"Missing option value: "+argument);
            return argv[++i];
        };
        if(argument=="--help" || argument=="-h") options.help=true;
        else if(argument=="--half-step") options.halfStep=true;
        else if(argument=="--output") { options.output=next(); need(!options.output.empty(),"Empty output directory."); }
        else if(argument=="--steps")
        {
            const auto text=next(); int value=0;
            const auto parsed=std::from_chars(text.data(),text.data()+text.size(),value);
            need(parsed.ec==std::errc{} && parsed.ptr==text.data()+text.size() && value>0,"--steps requires a positive integer.");
            options.stepLimit=value;
        }
        else if(argument=="--mode")
        {
            const auto value=next();
            need(value=="cpu" || value=="gpu" || value=="hybrid","Unknown backend: "+value);
            options.mode=value=="cpu" ? executionMode::CPU : value=="gpu" ? executionMode::GPU : executionMode::Hybrid;
        }
        else throw std::invalid_argument("Unknown option: "+argument);
    }
    return options;
}

int execute(SphereDEM& solver,const CaseSettings& settings,const std::vector<Body>& definitions,const RunOptions& options)
{
    const auto mode=options.mode;
    const bool half=options.halfStep;
    const int limit=options.stepLimit;
    const fs::path folder=options.output.empty() ? fs::current_path()/(std::string(settings.id)+"-results") : options.output;
    need(settings.steps>0 && math::isFinite(settings.timeStep) && settings.timeStep>0,"Invalid case duration.");
    need(solver.mode()==mode && solver.timeStep()==settings.timeStep/(half?2.0:1.0),"Case execution settings disagree with its C++ scene.");
    std::map<std::string,bool> ids;
    std::vector<RuntimeBody> bodies;
    for(const auto& value:definitions)
    {
        need(value.index>=0 && ids.emplace(value.id,true).second,"Invalid/duplicate body registration.");
        RuntimeBody body{value.id,value.sphere,value.index,{},{}};
        const auto& initial=state(solver,body);
        body.initialPosition=initial.position(); body.initialVelocity=initial.velocity();
        bodies.push_back(body);
    }
    need(settings.steps <= std::numeric_limits<int>::max() / (half ? 2 : 1), "Case step count overflows.");
    validateChecks(settings, ids, bodies);
    const int total = settings.steps * (half ? 2 : 1), steps = limit < 0 ? total : std::min(limit, total);
    // Strip empty final components without changing symlink/parent traversal semantics.
    fs::path reservation = folder;
    while (reservation.has_relative_path() && !reservation.has_filename())
        reservation = reservation.parent_path();
    need(reservation.filename() != "." && reservation.filename() != "..", "Output must name a new directory leaf, not . or ..");
    if (!reservation.parent_path().empty())
        fs::create_directories(reservation.parent_path());
    need(fs::create_directory(reservation), "Output directory already exists; choose a new directory.");
    std::ofstream trajectory(folder / "trajectory.csv"), contacts(folder / "contacts.csv"), comparisons(folder / "comparison.csv");
    need(bool(trajectory) && bool(contacts) && bool(comparisons), "Cannot create result files.");
    trajectory << std::setprecision(17);
    contacts << std::setprecision(17);
    comparisons << std::setprecision(17);
    trajectory << "step,time,id,kind,x,y,z,vx,vy,vz,wx,wy,wz,fx,fy,fz,tx,ty,tz,translational_kinetic_energy\n";
    contacts << "step,time,family,master,slave,overlap,normal_force,fx,fy,fz,elastic_energy\n";
    comparisons << "step,time,quantity,target,expected,actual,error_kind,error,tolerance,passed\n";
    bool passes = true;
    std::size_t evaluated = 0;
    const auto record = [&]
    {
        solver.observeCurrentState(
            [&](const fundem::solver&)
            {
                Snapshot frame;
                need(math::isFinite(solver.time()), "Non-finite simulation time during export.");
                for (const auto& value : bodies)
                {
                    const auto& body = state(solver, value);
                    validateBodyState(body, value.id);
                    const Real kinetic = body.inverseMass() > 0.0 ? 0.5 * squared(body.velocity()) / body.inverseMass() : 0.0;
                    need(math::isFinite(kinetic), "Non-finite translational kinetic energy: " + value.id);
                    trajectory << solver.stepCount() << ',' << solver.time() << ',' << value.id << ',' << (value.sphere ? "sphere" : "wall");
                    writeVector(trajectory, body.position());
                    writeVector(trajectory, body.velocity());
                    writeVector(trajectory, body.angularVelocity());
                    writeVector(trajectory, body.force());
                    writeVector(trajectory, body.torque());
                    trajectory << ',' << kinetic << '\n';
                    frame.translational += kinetic;
                    if (body.inverseMass() > 0.0)
                        frame.momentumX += body.velocity().x / body.inverseMass();
                }
                const auto collect = [&](const contactContainer& values, const char* family)
                {
                    for (const auto& value : values.host())
                    {
                        need(math::isFinite(value.normal()) && math::isFinite(value.force()) && math::isFinite(value.overlap()) &&
                             math::isFinite(value.normalElasticEnergy()) && math::isFinite(value.slidingElasticEnergy()) &&
                             math::isFinite(value.rollingElasticEnergy()) && math::isFinite(value.torsionalElasticEnergy()),
                             "Non-finite contact state during export.");
                        const Real normal = math::dot(value.force(), value.normal());
                        const Real energy = value.normalElasticEnergy() + value.slidingElasticEnergy() + value.rollingElasticEnergy() + value.torsionalElasticEnergy();
                        need(math::isFinite(normal) && math::isFinite(energy), "Non-finite contact observable during export.");
                        contacts << solver.stepCount() << ',' << solver.time() << ',' << family << ',' << value.masterParticleIndex() << ',' << value.slaveParticleIndex() << ',' << value.overlap()
                                 << ',' << normal;
                        writeVector(contacts, value.force());
                        contacts << ',' << energy << '\n';
                        frame.force += std::abs(normal);
                        frame.overlap = std::max(frame.overlap, value.overlap());
                        frame.elastic += energy;
                    }
                };
                collect(solver.sphereInteractions().contacts(), "sphere-sphere");
                collect(solver.sphereLSInteractions().contacts(), "sphere-LS");
                need(math::isFinite(frame.force) && math::isFinite(frame.overlap) && math::isFinite(frame.elastic) &&
                     math::isFinite(frame.translational) && math::isFinite(frame.momentumX), "Non-finite aggregate observable during export.");
                for (const auto& check : settings.checks)
                    if (check.step * (half ? 2 : 1) == solver.stepCount())
                    {
                        const Real actual = measure(check, solver, bodies, frame);
                        need(math::isFinite(actual), "Non-finite checkpoint observable: " + check.metric);
                        Real error = std::abs(actual - check.expected);
                        if (check.errorKind == "relative")
                        {
                            need(check.expected != 0.0, "Zero expected relative check.");
                            error /= std::abs(check.expected);
                        }
                        need(math::isFinite(error), "Non-finite checkpoint error: " + check.metric);
                        const bool passed = math::isFinite(actual) && error <= check.tolerance;
                        passes = passes && passed;
                        ++evaluated;
                        comparisons << solver.stepCount() << ',' << solver.time() << ',' << check.metric << ',' << check.target << ',' << check.expected << ',' << actual << ',' << check.errorKind
                                    << ',' << error << ',' << check.tolerance << ',' << (passed ? "true" : "false") << '\n';
                    }
            });
    };
    solver.initialize();
    record();
    for (int step = 0; step < steps; ++step)
    {
        solver.step();
        record();
    }
    trajectory.flush();
    contacts.flush();
    comparisons.flush();
    need(bool(trajectory) && bool(contacts) && bool(comparisons), "Failed writing comparison results.");
    const bool complete = steps == total;
    need(!complete || evaluated == settings.checks.size(), "Some checkpoints were not evaluated.");
    std::ofstream summary(folder / "summary.txt");
    need(bool(summary), "Cannot create summary.");
    summary << "Case: " << settings.id << "\nHistorical reported outcome: " << settings.historicalOutcome << "\nNative backend: "
            << (mode == executionMode::CPU   ? "cpu"
                : mode == executionMode::GPU ? "gpu"
                                             : "hybrid")
            << "\nCompleted steps: " << steps << " / " << total << "\nFresh configured checks evaluated: " << evaluated << " / " << settings.checks.size() << "\nFresh configured check outcome: "
            << (!complete              ? "INCOMPLETE (short run)"
                : settings.checks.empty() ? "NOT EVALUATED (use retained source criteria)"
                : passes               ? "PASS"
                                       : "FAIL")
            << "\nHistorical outcomes are preserved, not inferred from execution success.\n";
    summary.flush();
    need(bool(summary), "Failed writing summary.");
    std::cout << settings.id << ": " << (complete ? "COMPLETED" : "INCOMPLETE") << ", " << steps << " steps; historical=" << settings.historicalOutcome << "; configured checks=" << evaluated << '/'
              << settings.checks.size() << (complete && !settings.checks.empty() ? (passes ? " PASS" : " FAIL") : " (no full acceptance claim)") << '\n';
    return complete && !passes ? 2 : 0;
}
} // namespace fundem::comparison
