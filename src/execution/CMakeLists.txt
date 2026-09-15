add_library(FunDEM_execution INTERFACE)
add_library(FunDEM::execution ALIAS FunDEM_execution)

target_sources(
    FunDEM_execution
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_SOURCE_DIR}"
        FILES
            motionIntegration.h
            contactDetection.h
            contactFunctions.h
            bondFunctions.h
            sphFunctions.h
            sphCouplingFunctions.h
            sphJetFunctions.h
            sphSourceDiscretization.h
            cpu/parallelPolicy.h
            cpu/particleFunctions.h
)

target_include_directories(
    FunDEM_execution
    INTERFACE
        "$<BUILD_INTERFACE:${FUNDEM_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_execution INTERFACE cxx_std_17)
target_link_libraries(FunDEM_execution INTERFACE FunDEM::math FunDEM::data)
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
