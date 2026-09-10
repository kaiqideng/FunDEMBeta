add_library(FunDEM_interaction contact.cpp bond.cpp contactSearch.cpp)
add_library(FunDEM::interaction ALIAS FunDEM_interaction)

if(FUNDEM_HAS_CUDA)
    target_sources(FunDEM_interaction PRIVATE interactionContainer.cu)
endif()

target_sources(
    FunDEM_interaction
    PUBLIC
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES
            bond.h
            contact.h
            contactSearch.h
            neighborRange.h
            surfaceNodeMapping.h
            interactionContainer.h
)

target_include_directories(
    FunDEM_interaction
    PUBLIC
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_interaction PUBLIC cxx_std_17)
target_link_libraries(FunDEM_interaction PUBLIC FunDEM::math FunDEM::data FunDEM::particle)
target_link_libraries(FunDEM_interaction PRIVATE FunDEM::execution)

if(FUNDEM_HAS_OPENMP)
    target_link_libraries(FunDEM_interaction PRIVATE OpenMP::OpenMP_CXX)
endif()

if(FUNDEM_HAS_CUDA)
    set_target_properties(FunDEM_interaction PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED YES CUDA_EXTENSIONS NO)
endif()

set_target_properties(FunDEM_interaction PROPERTIES EXPORT_NAME interaction)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_interaction
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
