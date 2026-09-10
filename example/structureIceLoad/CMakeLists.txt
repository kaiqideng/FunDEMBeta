set(
    STRUCTURE_ICE_LOAD_SOURCES
    structureIceLoad.cpp
    structureIceLoadSolver.cpp
)

if(FUNDEM_HAS_CUDA)
    list(APPEND STRUCTURE_ICE_LOAD_SOURCES structureIceLoadKernel.cu)
endif()

add_executable(structureIceLoad ${STRUCTURE_ICE_LOAD_SOURCES})
target_link_libraries(structureIceLoad PRIVATE FunDEM::FunDEM)
target_include_directories(structureIceLoad PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(structureIceLoad PRIVATE cxx_std_17)

if(FUNDEM_HAS_CUDA)
    set_target_properties(
        structureIceLoad
        PROPERTIES
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED YES
            CUDA_EXTENSIONS NO
    )
    if(FUNDEM_CUDA_ARCHITECTURES)
        set_target_properties(structureIceLoad PROPERTIES CUDA_ARCHITECTURES "${FUNDEM_CUDA_ARCHITECTURES}")
    endif()
    if(FUNDEM_ENABLE_CUDA_DEVICE_DEBUG)
        target_compile_options(structureIceLoad PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-G;-g;-O0>")
    elseif(FUNDEM_ENABLE_CUDA_LINEINFO)
        target_compile_options(structureIceLoad PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>")
    endif()
endif()
