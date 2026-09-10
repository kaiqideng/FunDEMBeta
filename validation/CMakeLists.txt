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

add_custom_target(
    runValidations
    COMMAND "$<TARGET_FILE:freeFall>"
    COMMAND "$<TARGET_FILE:normalCollision>"
    COMMAND "$<TARGET_FILE:frictionLaw>"
    COMMAND "$<TARGET_FILE:bondModes>"
    COMMAND "$<TARGET_FILE:LSCentralCompression>"
    COMMAND "$<TARGET_FILE:backendConsistency>"
    DEPENDS
        freeFall
        normalCollision
        frictionLaw
        bondModes
        LSCentralCompression
        backendConsistency
    USES_TERMINAL
    VERBATIM
)
