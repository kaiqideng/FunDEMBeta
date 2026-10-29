add_library(FunDEM_execution cpu/pointGrid.cpp)
add_library(FunDEM::execution ALIAS FunDEM_execution)

target_sources(
    FunDEM_execution
    PUBLIC
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_SOURCE_DIR}"
        FILES
            motionIntegration.h
            contactDetection.h
            contactFunctions.h
            contactMaterialFunctions.h
            bondFunctions.h
            sphFunctions.h
            sphCouplingFunctions.h
            sphJetFunctions.h
            sphSourceDiscretization.h
            sphSpatialGridFunctions.h
            cpu/parallelPolicy.h
            cpu/particleFunctions.h
            cpu/pointGrid.h
)

target_include_directories(
    FunDEM_execution
    PUBLIC
        "$<BUILD_INTERFACE:${FUNDEM_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_execution PUBLIC cxx_std_17)
target_link_libraries(FunDEM_execution PUBLIC FunDEM::math FunDEM::data FunDEM::material)
if(FUNDEM_HAS_OPENMP)
    target_link_libraries(FunDEM_execution PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(FunDEM_execution PROPERTIES EXPORT_NAME execution)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_execution
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()

if(FUNDEM_HAS_CUDA)
    add_subdirectory(cuda)
endif()
