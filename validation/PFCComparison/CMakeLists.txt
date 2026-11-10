set(comparison_cases
    pfc-linear-restitution-md0-beta02
    pfc-linear-restitution-md0-beta06
    pfc-linear-restitution-md0-sweep
)
set(comparison_targets)
foreach(case_id IN LISTS comparison_cases)
    string(REPLACE "-" "_" target "pfc_${case_id}")
    add_executable(${target} "cases/${case_id}/${case_id}.cpp")
    target_link_libraries(${target} PRIVATE fundem_comparison_runtime)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/validation/PFCComparison")
    list(APPEND comparison_targets ${target})
endforeach()
add_custom_target(pfcComparison DEPENDS ${comparison_targets})
