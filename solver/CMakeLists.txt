set(FUNDEM_SOLVER_HEADERS solver.h LSDEM.h SphereDEM.h SPHDEM.h solvers.h "${FUNDEM_ROOT}/data/vtuOutput.h" "${FUNDEM_ROOT}/data/energyOutput.h")
set(
    FUNDEM_SOLVER_SOURCES
    solver.cpp
    LSDEM.cpp
    SphereDEM.cpp
    SPHDEM.cpp
    "${FUNDEM_ROOT}/data/vtuOutput.cpp"
    "${FUNDEM_ROOT}/data/energyOutput.cpp"
)

add_library(FunDEM_solver ${FUNDEM_SOLVER_SOURCES})
add_library(FunDEM::solver ALIAS FunDEM_solver)

target_sources(
    FunDEM_solver
    PUBLIC
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES ${FUNDEM_SOLVER_HEADERS}
)

target_include_directories(
    FunDEM_solver
    PUBLIC
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_solver PUBLIC cxx_std_17)
target_compile_definitions(FunDEM_solver PRIVATE FUNDEM_HAS_CUDA=$<BOOL:${FUNDEM_HAS_CUDA}>)
target_link_libraries(FunDEM_solver PUBLIC FunDEM::math FunDEM::dataIO FunDEM::material FunDEM::geometry FunDEM::particle FunDEM::interaction)

if(FUNDEM_HAS_OPENMP)
    target_link_libraries(FunDEM_solver PRIVATE OpenMP::OpenMP_CXX)
endif()

if(FUNDEM_HAS_CUDA)
    target_link_libraries(FunDEM_solver PUBLIC FunDEM::cuda CUDA::cudart)
endif()

set_target_properties(FunDEM_solver PROPERTIES EXPORT_NAME solver)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_solver
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
