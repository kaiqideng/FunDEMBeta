#include "solver.h"

#include "data/energyOutput.h"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
#include "data/HostAoSDeviceSoA.h"
#endif

namespace fundem
{
namespace
{

namespace fs = std::filesystem;

fs::path executableDirectory()
{
#if defined(_WIN32)
    std::wstring path(32768, L'\0');
    const DWORD size = GetModuleFileNameW(nullptr, path.data(), static_cast<DWORD>(path.size()));
    if (size > 0 && size < path.size())
    {
        path.resize(size);
        return fs::path(path).parent_path();
    }
#elif defined(__APPLE__)
    std::uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::vector<char> path(size);
    if (_NSGetExecutablePath(path.data(), &size) == 0)
    {
        return fs::weakly_canonical(fs::path(path.data())).parent_path();
    }
#elif defined(__linux__)
    std::error_code error;
    const fs::path path = fs::read_symlink("/proc/self/exe", error);
    if (!error && !path.empty())
    {
        return path.parent_path();
    }
#endif
    return fs::current_path();
}

fs::path resolveOutputDirectory(const std::string& directory)
{
    fs::path path(directory);
    if (path.is_relative())
    {
        path = executableDirectory() / path;
    }
    return path.lexically_normal();
}

bool isSolverOutputFile(const fs::path& path) { return path.extension() == ".vtu" || path.extension() == ".dat"; }

} // namespace

#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA

struct solver::gpuState {
    explicit gpuState(int device) : device_(device)
    {
        host_device_detail::checkCuda(cudaSetDevice(device_), "cudaSetDevice");
        host_device_detail::checkCuda(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking), "cudaStreamCreateWithFlags");
    }

    ~gpuState()
    {
        if (stream_ != nullptr)
        {
            cudaSetDevice(device_);
            cudaStreamDestroy(stream_);
        }
    }

    int device_{0};
    cudaStream_t stream_{nullptr};
};

#else

struct solver::gpuState {};

#endif

solver::solver(int device, executionMode mode) : device_(device), executionMode_(mode)
{
    if (device_ < 0)
    {
        throw std::invalid_argument("The device index cannot be negative.");
    }
}

solver::~solver() = default;

void solver::setExecutionMode(executionMode value)
{
    if (initialized_)
    {
        throw std::logic_error("The execution mode cannot be changed while the solver is initialized.");
    }
    validateExecutionMode(value);
    if (value != executionMode::CPU)
    {
        validateGPUEnable();
    }
    if (value == executionMode_)
    {
        return;
    }
    executionMode_ = value;
    executionModeChanged(value);
    advanceTemplate_ = nullptr;
    invalidateContinuation();
}

void solver::setBoundary(const Vec3& minimumBoundary, const Vec3& maximumBoundary)
{
    if (initialized_)
    {
        throw std::logic_error("The solver boundary cannot be changed after initialization.");
    }
    if (!math::isFinite(minimumBoundary) || !math::isFinite(maximumBoundary) || maximumBoundary.x <= minimumBoundary.x || maximumBoundary.y <= minimumBoundary.y ||
        maximumBoundary.z <= minimumBoundary.z)
    {
        throw std::invalid_argument("Invalid solver boundary.");
    }
    if (minimumBoundary == minimumBoundary_ && maximumBoundary == maximumBoundary_)
    {
        return;
    }
    minimumBoundary_ = minimumBoundary;
    maximumBoundary_ = maximumBoundary;
    invalidateContinuation();
}

void solver::setGravity(const Vec3& value)
{
    if (initialized_)
    {
        throw std::logic_error("Gravity cannot be changed while the solver is initialized.");
    }
    if (!math::isFinite(value))
    {
        throw std::invalid_argument("Gravity must be finite.");
    }
    if (value == gravity_)
    {
        return;
    }
    gravity_ = value;
    invalidateContinuation();
}

void solver::setTimeStep(Real value)
{
    if (initialized_)
    {
        throw std::logic_error("The time step cannot be changed while the solver is initialized.");
    }
    if (!math::isFinite(value) || value <= 0.0)
    {
        throw std::invalid_argument("The time step must be finite and positive.");
    }
    if (value == timeStep_)
    {
        return;
    }
    timeStep_ = value;
    invalidateContinuation();
}

void solver::setOutputDirectory(std::string directory)
{
    if (directory.empty())
    {
        throw std::invalid_argument("Solver output directory must not be empty.");
    }
    const std::string resolvedDirectory = resolveOutputDirectory(directory).string();
    if (resolvedDirectory != outputDirectory_)
    {
        outputFrameCount_ = 0;
        nextOutputStep_ = stepCount_;
    }
    outputDirectory_ = resolvedDirectory;
}

void solver::setOutputStepInterval(int stepInterval)
{
    if (stepInterval <= 0)
    {
        throw std::invalid_argument("Output step interval must be positive.");
    }
    if (outputDirectory_.empty())
    {
        throw std::logic_error("Set the solver output directory before enabling output.");
    }
    outputStepInterval_ = stepInterval;
    nextOutputStep_ = stepCount_;
}

std::string solver::makeVTUFileName(const char* group, const char* name, int frameIndex) const
{
    std::ostringstream fileName;
    fileName << name << '_' << std::setw(6) << std::setfill('0') << frameIndex << ".vtu";
    return (std::filesystem::path(outputDirectory_) / group / fileName.str()).string();
}

void solver::writeOutput()
{
    if (outputDirectory_.empty())
    {
        throw std::logic_error("Output has not been configured.");
    }
    ensureOutputDirectory();
    observeCurrentState([this](const solver&)
    {
        const int frameIndex = outputFrameCount_;
        writeSystemVTU(frameIndex);
        appendEnergyDAT((std::filesystem::path(outputDirectory_) / "energy.dat").string(), time_, systemEnergy(), frameIndex > 0);
        std::cout << "[Solver] Calculated steps: " << stepCount_ << ", current time: " << time_ << ", current frame: " << frameIndex << ", device memory: " << std::fixed << std::setprecision(6)
                  << deviceMemoryGB() << " GB" << std::defaultfloat << std::endl;
    });
    ++outputFrameCount_;
    nextOutputStep_ = stepCount_ + outputStepInterval_;
}

void solver::observeCurrentState(const std::function<void(const solver&)>& observer)
{
    if (!observer)
    {
        throw std::invalid_argument("A current-state observer is required.");
    }
    if (observingCurrentState_)
    {
        observer(*this);
        return;
    }
    activateGPUDevice();
    observingCurrentState_ = true;
    try
    {
        beginOutputSnapshot();
        synchronizeHostState();
        observer(*this);
    }
    catch (...)
    {
        observingCurrentState_ = false;
        endOutputSnapshot();
        throw;
    }
    observingCurrentState_ = false;
    endOutputSnapshot();
}

void solver::writeOutputIfDue()
{
    if (outputStepInterval_ > 0 && stepCount_ >= nextOutputStep_)
    {
        writeOutput();
    }
}

void solver::requireNoCurrentStateObservation() const
{
    if (observingCurrentState_)
        throw std::logic_error("The solver cannot initialize or advance during a current-state observation.");
}

void solver::solve(int numberOfSteps)
{
    requireNoCurrentStateObservation();
    if (numberOfSteps < 0)
    {
        throw std::invalid_argument("The number of steps cannot be negative.");
    }
    if (firstSolve_)
    {
        if (!outputDirectory_.empty())
        {
            clearOutputFiles();
        }
        firstSolve_ = false;
    }
    if (!initialized_)
    {
        initialize();
    }
    writeOutputIfDue();
    for (int stepIndex = 0; stepIndex < numberOfSteps; ++stepIndex)
    {
        step();
        writeOutputIfDue();
    }
    activateGPUDevice();
    try
    {
        beginSolveCompletionSnapshot();
        synchronizeHostState();
    }
    catch (...)
    {
        endSolveCompletionSnapshot();
        throw;
    }
    endSolveCompletionSnapshot();
    continuationValid_ = true;
    initialized_ = false;
}

void solver::initialize()
{
    requireNoCurrentStateObservation();
    if (initialized_)
    {
        return;
    }
    validateConfiguration();
    if (usesDevice())
    {
        initializeDeviceContext();
    }
    if (continuationValid_ && advanceTemplate_ != nullptr)
    {
        continuationValid_ = false;
        resumeContinuation();
        synchronize();
        initialized_ = true;
        return;
    }
    continuationValid_ = false;
    switch (executionMode_)
    {
    case executionMode::CPU:
        initialize(cpu::mode{});
        advanceTemplate_ = &solver::advanceCPU;
        break;
    case executionMode::GPU:
        initialize(gpu::mode{}, gpuStream());
        advanceTemplate_ = &solver::advanceGPU;
        break;
    case executionMode::Hybrid:
        initialize(hybrid::mode{}, gpuStream());
        advanceTemplate_ = &solver::advanceHybrid;
        break;
    }
    synchronize();
    initialized_ = true;
}

void solver::ensureOutputDirectory() const { fs::create_directories(outputDirectory_); }

void solver::clearOutputFiles()
{
    ensureOutputDirectory();
    const fs::path directory(outputDirectory_);
    for (const fs::directory_entry& entry : fs::recursive_directory_iterator(directory))
    {
        if (entry.is_regular_file() && isSolverOutputFile(entry.path()))
        {
            fs::remove(entry.path());
        }
    }
}

void solver::step()
{
    requireNoCurrentStateObservation();
    if (!initialized_)
    {
        initialize();
    }
    activateGPUDevice();
    (this->*advanceTemplate_)(timeStep_);
    time_ += timeStep_;
    ++stepCount_;
}

void solver::validateConfiguration() const
{
    if (!math::isFinite(timeStep_) || timeStep_ <= 0.0 || !math::isFinite(gravity_) || !math::isFinite(minimumBoundary_) || !math::isFinite(maximumBoundary_) ||
        maximumBoundary_.x <= minimumBoundary_.x || maximumBoundary_.y <= minimumBoundary_.y || maximumBoundary_.z <= minimumBoundary_.z)
    {
        throw std::invalid_argument("Invalid solver configuration.");
    }
    validateExecutionMode(executionMode_);
    validateSystemConfiguration();
}

void solver::invalidateContinuation() noexcept
{
    continuationValid_ = false;
    invalidateDeferredState();
}

void solver::validateGPUEnable() const
{
    if (initialized_)
    {
        throw std::logic_error("GPU execution cannot be enabled after solver initialization.");
    }
#if !defined(FUNDEM_HAS_CUDA) || !FUNDEM_HAS_CUDA
    throw std::runtime_error("The GPU backend is unavailable in this FunDEM build.");
#endif
}

void solver::initializeDeviceContext()
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (!gpu_)
    {
        gpu_ = std::make_unique<gpuState>(device_);
    }
    else
    {
        host_device_detail::checkCuda(cudaSetDevice(device_), "cudaSetDevice");
    }
#else
    throw std::runtime_error("The GPU backend is unavailable in this FunDEM build.");
#endif
}

cudaStream_t solver::gpuStream() const noexcept
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    return gpu_ ? gpu_->stream_ : nullptr;
#else
    return nullptr;
#endif
}

void solver::activateGPUDevice() const
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (gpu_)
    {
        host_device_detail::checkCuda(cudaSetDevice(device_), "cudaSetDevice");
    }
#endif
}

void solver::synchronize()
{
#if defined(FUNDEM_HAS_CUDA) && FUNDEM_HAS_CUDA
    if (gpu_)
    {
        host_device_detail::checkCuda(cudaSetDevice(device_), "cudaSetDevice");
        host_device_detail::synchronize(gpu_->stream_);
    }
#endif
}

void solver::copyDeviceToHost()
{
    if (!initialized_ || !usesDevice())
    {
        return;
    }
    activateGPUDevice();
    copySystemDeviceToHost(gpuStream());
    synchronize();
}

void solver::synchronizeHostState()
{
    if (usesDevice() && (initialized_ || continuationValid_))
    {
        // solve() leaves a resumable solver uninitialized. Its observation
        // projection still owns valid device state that must reach the host;
        // public copyDeviceToHost() intentionally skips uninitialized models.
        copySystemDeviceToHost(gpuStream());
    }
    synchronize();
}

double solver::deviceMemoryGB() const noexcept { return systemDeviceMemoryGB(); }

} // namespace fundem
