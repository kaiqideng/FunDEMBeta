# The full reference run remains an explicit command.
add_library(fundem_ondem_shear_model STATIC periodicContactForce.cpp periodicShear.cpp)
target_link_libraries(fundem_ondem_shear_model PUBLIC FunDEM::FunDEM)
target_compile_features(fundem_ondem_shear_model PUBLIC cxx_std_17)

add_executable(ondemShearFlow main.cpp benchmarkIO.cpp)
add_executable(ondemShearFlowCheck periodicShearTest.cpp)
add_executable(ondemShearFlowIOCheck benchmarkIOTest.cpp benchmarkIO.cpp)
add_executable(ondemShearFlowExternalForceCheck periodicContactForceTest.cpp)

set(ondem_programs ondemShearFlow ondemShearFlowCheck ondemShearFlowIOCheck ondemShearFlowExternalForceCheck)
foreach(target IN LISTS ondem_programs)
    target_link_libraries(${target} PRIVATE fundem_ondem_shear_model)
    set_target_properties(${target} PROPERTIES
        CXX_EXTENSIONS NO RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/validation")
endforeach()
foreach(target IN ITEMS fundem_ondem_shear_model ${ondem_programs})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    elseif(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    endif()
endforeach()

add_custom_target(checkOndemShearFlow
    COMMAND "$<TARGET_FILE:ondemShearFlowCheck>"
    COMMAND "$<TARGET_FILE:ondemShearFlowIOCheck>"
    COMMAND "$<TARGET_FILE:ondemShearFlowExternalForceCheck>"
    DEPENDS ondemShearFlowCheck ondemShearFlowIOCheck ondemShearFlowExternalForceCheck
    USES_TERMINAL VERBATIM)
