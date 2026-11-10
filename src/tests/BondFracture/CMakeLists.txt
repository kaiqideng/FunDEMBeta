set(bond_fracture_targets)
foreach(case_id tension shear bending torsion)
    add_executable(bond_fracture_${case_id} ${case_id}.cpp)
    target_link_libraries(bond_fracture_${case_id} PRIVATE FunDEM::FunDEM)
    fundem_enable_warnings(bond_fracture_${case_id})
    set_target_properties(bond_fracture_${case_id} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    list(APPEND bond_fracture_targets bond_fracture_${case_id})
endforeach()
add_executable(bond_fracture_bending_zero_torque bendingZeroTorque.cpp)
target_link_libraries(bond_fracture_bending_zero_torque PRIVATE FunDEM::FunDEM)
fundem_enable_warnings(bond_fracture_bending_zero_torque)
set_target_properties(bond_fracture_bending_zero_torque PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
list(APPEND bond_fracture_targets bond_fracture_bending_zero_torque)
add_custom_target(bondFractureTests DEPENDS ${bond_fracture_targets})

# Cheap prescribed-pose completion checks are distinct from the strict
# h/h2/h4 campaign. Each invocation retains a fresh build-local output.
foreach(target IN LISTS bond_fracture_targets)
    foreach(backend cpu gpu)
        if(backend STREQUAL "gpu" AND NOT FUNDEM_HAS_CUDA)
            continue()
        endif()
        string(REPLACE "bond_fracture_" "" case_id "${target}")
        string(TOUPPER "${backend}" backend_label)
        set(backend_test_label "${backend}")
        if(backend STREQUAL "gpu")
            set(backend_test_label cuda)
        endif()
        set(test_name "FunDEM.BondFracture.Smoke.${case_id}.${backend_label}")
        add_test(NAME "${test_name}" COMMAND "${CMAKE_COMMAND}"
            "-DBOND_EXECUTABLE=$<TARGET_FILE:${target}>"
            "-DBOND_BACKEND=${backend}"
            "-DBOND_CASE=${case_id}"
            "-DBOND_OUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR}/ctest/$<CONFIG>"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/runSmoke.cmake")
        set_tests_properties("${test_name}" PROPERTIES
            LABELS "unit;constitutive;smoke;${backend_test_label}" TIMEOUT 90
            ENVIRONMENT "OMP_NUM_THREADS=2")
        if(backend STREQUAL "gpu")
            set_tests_properties("${test_name}" PROPERTIES RESOURCE_LOCK fundem_cuda_device)
        endif()
    endforeach()
endforeach()
