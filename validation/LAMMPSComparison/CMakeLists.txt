set(comparison_cases
    dem01-hertz-tsuji
    dem01-hertz-twosphere
    dem01-hooke-2d-lj
    dem01-hooke-massvel
    dem01-hooke-twosphere
    dem01-hooke-twosphere-e02
    dem01-hooke-twosphere-e08
    dem01-hooke-veldamp
    dem01-hookehistory-twosphere
    dem01-legacy-hooke-twosphere
    dem01-legacy-hookehistory-twosphere
    dem02-hertz-peak
    dem02-hertz-twosphere
    dem02-hertzmaterial-twosphere
    dem02-legacy-hertzhistory-twosphere
    dem03-hertz-peak-wall
    dem03-hertz-tsuji-wall
    dem03-hertz-wall
    dem03-hertzmaterial-wall
    dem03-legacy-hertzhistory-wall
    dem04-hertz-oblique
    dem04-hookehistory-oblique
    dem05-hookehistory-slip
    dem05-mindlin-slip
    dem05-slip-mu02
    dem05-sliptime
    dem06-hookehistory-spin
    dem06-spin-identical
    dem07-hooke-rolling
    dem09-viscous-linear
    dem10-freefall
    dem10-stack
    dem11-regionwall
    dem13-twist-sds
    dem15-hooke-oblique
)
set(comparison_targets)
foreach(case_id IN LISTS comparison_cases)
    string(REPLACE "-" "_" target "lammps_${case_id}")
    add_executable(${target} "cases/${case_id}/${case_id}.cpp")
    target_link_libraries(${target} PRIVATE fundem_comparison_runtime)
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/validation/LAMMPSComparison")
    list(APPEND comparison_targets ${target})
endforeach()
add_custom_target(lammpsComparison DEPENDS ${comparison_targets})
