add_library(
    FunDEM_cuda
    integrationKernel.cu
    contactForceKernel.cu
    bondForceKernel.cu
    contactDetectionKernel.cu
    spatialGridKernel.cu
    sphInteractionKernel.cu
    sphJetKernel.cu
    atomicKernel.cuh
)
add_library(FunDEM::cuda ALIAS FunDEM_cuda)

target_sources(
    FunDEM_cuda
    PUBLIC
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES
            integrationKernel.cuh
            contactForceKernel.cuh
            bondForceKernel.cuh
            contactDetectionKernel.cuh
            spatialGridKernel.cuh
            sphInteractionKernel.cuh
            sphJetKernel.cuh
)

target_include_directories(
    FunDEM_cuda
    PUBLIC
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_cuda PUBLIC cxx_std_17)
target_link_libraries(FunDEM_cuda PUBLIC FunDEM::execution FunDEM::data FunDEM::interaction CUDA::cudart)

set_target_properties(
    FunDEM_cuda
    PROPERTIES
        EXPORT_NAME cuda
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED YES
        CUDA_EXTENSIONS NO
)

if(FUNDEM_CUDA_ARCHITECTURES)
    set_target_properties(FunDEM_cuda PROPERTIES CUDA_ARCHITECTURES "${FUNDEM_CUDA_ARCHITECTURES}")
endif()

if(FUNDEM_ENABLE_CUDA_DEVICE_DEBUG)
    target_compile_options(FunDEM_cuda PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-G;-g;-O0>")
elseif(FUNDEM_ENABLE_CUDA_LINEINFO)
    target_compile_options(FunDEM_cuda PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>")
endif()

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_cuda
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
