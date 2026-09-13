add_executable(tutorial0 tutorial0.cpp)
target_link_libraries(tutorial0 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial0 PRIVATE cxx_std_17)

add_executable(tutorial1 tutorial1.cpp)
target_link_libraries(tutorial1 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial1 PRIVATE cxx_std_17)

add_executable(tutorial2 tutorial2.cpp)
target_link_libraries(tutorial2 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial2 PRIVATE cxx_std_17)
if(FUNDEM_HAS_CUDA)
    target_sources(tutorial2 PRIVATE tutorial2DampingKernel.cu tutorial2DampingKernel.cuh)
    set_target_properties(
        tutorial2
        PROPERTIES
            CUDA_STANDARD 17
            CUDA_STANDARD_REQUIRED YES
            CUDA_EXTENSIONS NO
    )
    if(FUNDEM_CUDA_ARCHITECTURES)
        set_target_properties(tutorial2 PROPERTIES CUDA_ARCHITECTURES "${FUNDEM_CUDA_ARCHITECTURES}")
    endif()
endif()

add_executable(tutorial3 tutorial3.cpp)
target_link_libraries(tutorial3 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial3 PRIVATE cxx_std_17)

add_custom_command(
    TARGET tutorial0
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_LIST_DIR}/Gomboc.obj" "$<TARGET_FILE_DIR:tutorial0>/Gomboc.obj"
)
