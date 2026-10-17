function(fundem_add_validation target source)
    add_executable("${target}" "${source}")
    target_link_libraries("${target}" PRIVATE FunDEM::FunDEM)
    target_include_directories("${target}" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    set_target_properties(
        "${target}"
        PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/validation"
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options("${target}" PRIVATE -Wall -Wextra -Wpedantic)
    elseif(MSVC)
        target_compile_options("${target}" PRIVATE /W4)
    endif()
endfunction()

fundem_add_validation(freeFall freeFall.cpp)
fundem_add_validation(normalCollision normalCollision.cpp)
fundem_add_validation(frictionLaw frictionLaw.cpp)
fundem_add_validation(bondModes bondModes.cpp)
fundem_add_validation(LSCentralCompression LSCentralCompression.cpp)
fundem_add_validation(backendConsistency backendConsistency.cpp)

# This benchmark has its own explicitly periodic CPU pair law. Its full run is
# opt-in and is never part of runValidations or CTest.
add_library(fundem_ondem_shear_model STATIC
    ondemShearFlow/periodicContactForce.cpp
    ondemShearFlow/periodicShear.cpp)
target_link_libraries(fundem_ondem_shear_model PUBLIC FunDEM::FunDEM)
target_compile_features(fundem_ondem_shear_model PUBLIC cxx_std_17)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fundem_ondem_shear_model PRIVATE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(fundem_ondem_shear_model PRIVATE /W4)
endif()
fundem_add_validation(ondemShearFlow ondemShearFlow/main.cpp)
target_sources(ondemShearFlow PRIVATE ondemShearFlow/benchmarkIO.cpp)
target_link_libraries(ondemShearFlow PRIVATE fundem_ondem_shear_model)
fundem_add_validation(ondemShearFlowCheck ondemShearFlow/periodicShearTest.cpp)
target_link_libraries(ondemShearFlowCheck PRIVATE fundem_ondem_shear_model)
fundem_add_validation(ondemShearFlowIOCheck ondemShearFlow/benchmarkIOTest.cpp)
target_sources(ondemShearFlowIOCheck PRIVATE ondemShearFlow/benchmarkIO.cpp)
target_link_libraries(ondemShearFlowIOCheck PRIVATE fundem_ondem_shear_model)
fundem_add_validation(ondemShearFlowExternalForceCheck ondemShearFlow/periodicContactForceTest.cpp)
target_link_libraries(ondemShearFlowExternalForceCheck PRIVATE fundem_ondem_shear_model)

add_custom_target(checkOndemShearFlow
    COMMAND "$<TARGET_FILE:ondemShearFlowCheck>"
    COMMAND "$<TARGET_FILE:ondemShearFlowIOCheck>"
    COMMAND "$<TARGET_FILE:ondemShearFlowExternalForceCheck>"
    DEPENDS ondemShearFlowCheck ondemShearFlowIOCheck ondemShearFlowExternalForceCheck
    USES_TERMINAL VERBATIM)

add_custom_target(
    runValidations
    COMMAND "$<TARGET_FILE:freeFall>"
    COMMAND "$<TARGET_FILE:normalCollision>"
    COMMAND "$<TARGET_FILE:frictionLaw>"
    COMMAND "$<TARGET_FILE:bondModes>"
    COMMAND "$<TARGET_FILE:LSCentralCompression>"
    COMMAND "$<TARGET_FILE:backendConsistency>"
    COMMAND "$<TARGET_FILE:ondemShearFlowCheck>"
    COMMAND "$<TARGET_FILE:ondemShearFlowIOCheck>"
    COMMAND "$<TARGET_FILE:ondemShearFlowExternalForceCheck>"
    DEPENDS
        freeFall
        normalCollision
        frictionLaw
        bondModes
        LSCentralCompression
        backendConsistency
        ondemShearFlowCheck
        ondemShearFlowIOCheck
        ondemShearFlowExternalForceCheck
    USES_TERMINAL
    VERBATIM
)
