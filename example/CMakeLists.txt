add_subdirectory(structureIceLoad)

if(NOT FUNDEM_HAS_CUDA)
    message(STATUS "The SPH examples are skipped because CUDA is unavailable.")
    return()
endif()

function(add_fundem_cuda_example target)
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} PRIVATE FunDEM::FunDEM)
    target_compile_features(${target} PRIVATE cxx_std_17)
    set_source_files_properties(${target}.cpp PROPERTIES LANGUAGE CUDA)
    set_target_properties(
        ${target}
        PROPERTIES
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED YES
            CUDA_EXTENSIONS NO
    )
    if(FUNDEM_CUDA_ARCHITECTURES)
        set_target_properties(${target} PROPERTIES CUDA_ARCHITECTURES "${FUNDEM_CUDA_ARCHITECTURES}")
    endif()
endfunction()

add_fundem_cuda_example(zhang2017ObstacleRWCSPH)
add_fundem_cuda_example(canelas2016StackedCubesRWCSPH)
