cmake_minimum_required(VERSION 3.24)

project(
    FunDEM
    VERSION 1.0.2
    DESCRIPTION "CPU/GPU discrete element method framework"
    LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckLanguage)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(FUNDEM_ENABLE_CUDA "Enable the CUDA backend when available" ON)
option(FUNDEM_ENABLE_CUDA_DEVICE_DEBUG "Enable CUDA device debug flags" OFF)
option(FUNDEM_ENABLE_CUDA_LINEINFO "Add CUDA line information" ON)
option(FUNDEM_ENABLE_OPENMP "Enable OpenMP for host parallelism when available" ON)
option(FUNDEM_INSTALL "Generate FunDEM install and package targets" ON)
option(FUNDEM_BUILD_TUTORIALS "Build FunDEM tutorials" ON)
option(FUNDEM_BUILD_EXAMPLES "Build the native FunDEM example template" ${PROJECT_IS_TOP_LEVEL})
option(FUNDEM_BUILD_VALIDATIONS "Build FunDEM validation cases" ON)

set(
    FUNDEM_CUDA_ARCHITECTURES
    ""
    CACHE STRING
    "CUDA architectures to build; empty uses the CMake/compiler default"
)

set(FUNDEM_HAS_CUDA OFF)
if(FUNDEM_ENABLE_CUDA)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
        set(FUNDEM_HAS_CUDA ON)
    else()
        message(
            STATUS
            "CUDA compiler not found: the CPU modules will still be built."
        )
    endif()
endif()

set(FUNDEM_HAS_OPENMP OFF)
if(FUNDEM_ENABLE_OPENMP)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        include(CheckCXXSourceCompiles)
        include(CMakePushCheckState)
        cmake_push_check_state(RESET)
        set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_CXX)
        # Recheck after compiler/runtime option changes in an existing build.
        unset(FUNDEM_OPENMP_REDUCTIONS_SUPPORTED CACHE)
        check_cxx_source_compiles([=[
            #include <omp.h>
            int main() {
                double minimum = 100.0, maximum = -100.0;
                #pragma omp parallel for reduction(min:minimum) reduction(max:maximum)
                for (int i = 0; i < 16; ++i) {
                    if (i < minimum) minimum = i;
                    if (i > maximum) maximum = i;
                }
                return minimum == 0.0 && maximum == 15.0 ? 0 : 1;
            }
        ]=] FUNDEM_OPENMP_REDUCTIONS_SUPPORTED)
        cmake_pop_check_state()
        if(FUNDEM_OPENMP_REDUCTIONS_SUPPORTED)
            set(FUNDEM_HAS_OPENMP ON)
        else()
            message(STATUS "OpenMP lacks required min/max reductions: host loops will use serial execution.")
        endif()
    else()
        message(STATUS "OpenMP not found: host loops will use serial execution.")
    endif()
endif()

set(FUNDEM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

add_subdirectory(src)

add_library(FunDEM INTERFACE)
add_library(FunDEM::FunDEM ALIAS FunDEM)

target_include_directories(
    FunDEM
    INTERFACE "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_compile_features(FunDEM INTERFACE cxx_std_17)
target_compile_definitions(
    FunDEM
    INTERFACE FUNDEM_HAS_CUDA=$<BOOL:${FUNDEM_HAS_CUDA}>
)
target_link_libraries(
    FunDEM
    INTERFACE
        FunDEM::math
        FunDEM::dataLayout
        FunDEM::data
        FunDEM::dataIO
        FunDEM::material
        FunDEM::geometry
        FunDEM::particle
        FunDEM::interaction
        FunDEM::execution
        FunDEM::solver
)

if(FUNDEM_HAS_CUDA)
    target_link_libraries(FunDEM INTERFACE FunDEM::cuda)
endif()

set_target_properties(FunDEM PROPERTIES EXPORT_NAME FunDEM)

if(FUNDEM_BUILD_TUTORIALS)
    add_subdirectory(tutorial)
endif()

if(FUNDEM_BUILD_EXAMPLES)
    add_subdirectory(example)
endif()

if(FUNDEM_BUILD_VALIDATIONS)
    add_subdirectory(validation)
endif()

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(src/tests)
endif()

if(FUNDEM_INSTALL)
    install(TARGETS FunDEM EXPORT FunDEMTargets)

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FunDEMConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/FunDEMConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FunDEM"
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/FunDEMConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}"
        COMPATIBILITY SameMajorVersion
    )

    install(
        FILES
            "${CMAKE_CURRENT_BINARY_DIR}/FunDEMConfig.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/FunDEMConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FunDEM"
    )

    install(
        EXPORT FunDEMTargets
        FILE FunDEMTargets.cmake
        NAMESPACE FunDEM::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FunDEM"
    )
endif()
