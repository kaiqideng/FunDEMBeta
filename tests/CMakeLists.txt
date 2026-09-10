function(fundem_enable_warnings target)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(
            "${target}"
            PRIVATE
                "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>"
                "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall,-Wextra>"
        )
    elseif(MSVC)
        target_compile_options(
            "${target}"
            PRIVATE
                "$<$<COMPILE_LANGUAGE:CXX>:/W4>"
        )
    endif()
endfunction()

function(fundem_add_cpp_test target test_name source)
    add_executable("${target}" "${source}")
    target_link_libraries("${target}" PRIVATE FunDEM::FunDEM)
    fundem_enable_warnings("${target}")
    set_target_properties(
        "${target}"
        PROPERTIES
            CXX_STANDARD 17
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
    )
    add_test(NAME "${test_name}" COMMAND "${target}")
    set_tests_properties("${test_name}" PROPERTIES LABELS "unit;cpu")
endfunction()

fundem_add_cpp_test(
    fundem_physics_test
    FunDEM.Physics
    PhysicsTest.cpp
)

fundem_add_cpp_test(
    fundem_sph_cpu_state_test
    FunDEM.SPHState.CPU
    SPHCPUStateTest.cpp
)

fundem_add_cpp_test(fundem_sph_robustness_test FunDEM.SPHRobustness SPHRobustnessTest.cpp)

fundem_add_cpp_test(
    fundem_sph_jet_test
    FunDEM.SPHJet
    SPHJetTest.cpp
)

fundem_add_cpp_test(
    fundem_sph_coupling_test
    FunDEM.SPHCoupling
    SPHCouplingTest.cpp
)

fundem_add_cpp_test(
    fundem_fiber_test
    FunDEM.Fiber
    FiberTest.cpp
)
set_tests_properties(FunDEM.Fiber PROPERTIES LABELS "integration;cpu" TIMEOUT 300)

if(FUNDEM_HAS_CUDA)
    function(fundem_add_cuda_test target test_name source)
        add_executable("${target}" "${source}")
        target_link_libraries("${target}" PRIVATE FunDEM::FunDEM)
        fundem_enable_warnings("${target}")

        set_target_properties(
            "${target}"
            PROPERTIES
                CXX_STANDARD 17
                CXX_STANDARD_REQUIRED YES
                CXX_EXTENSIONS NO
                CUDA_STANDARD 17
                CUDA_STANDARD_REQUIRED YES
                CUDA_EXTENSIONS NO
        )

        if(FUNDEM_CUDA_ARCHITECTURES)
            set_target_properties(
                "${target}"
                PROPERTIES
                    CUDA_ARCHITECTURES "${FUNDEM_CUDA_ARCHITECTURES}"
            )
        endif()

        if(FUNDEM_ENABLE_CUDA_DEVICE_DEBUG)
            target_compile_options(
                "${target}"
                PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-G;-g;-O0>"
            )
        elseif(FUNDEM_ENABLE_CUDA_LINEINFO)
            target_compile_options(
                "${target}"
                PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>"
            )
        endif()

        add_test(NAME "${test_name}" COMMAND "${target}")
        set_tests_properties("${test_name}" PROPERTIES LABELS "unit;cuda")
    endfunction()

    fundem_add_cuda_test(
        fundem_math_test
        FunDEM.Math
        MathTest.cu
    )

    fundem_add_cuda_test(
        fundem_host_device_test
        FunDEM.HostAoSDeviceSoA
        HostAoSDeviceSoATest.cu
    )

    fundem_add_cuda_test(
        fundem_sph_state_test
        FunDEM.SPHState
        SPHStateTest.cu
    )
    fundem_add_cuda_test(
        fundem_sph_robustness_gpu_test
        FunDEM.SPHRobustness.GPU
        SPHRobustnessTest.cu
    )
else()
    message(STATUS "CUDA tests are skipped because CUDA is unavailable.")
endif()
