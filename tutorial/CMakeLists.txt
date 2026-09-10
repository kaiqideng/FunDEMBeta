add_executable(tutorial0 tutorial0.cpp)
target_link_libraries(tutorial0 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial0 PRIVATE cxx_std_17)

add_executable(tutorial1 tutorial1.cpp)
target_link_libraries(tutorial1 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial1 PRIVATE cxx_std_17)

add_executable(tutorial2 tutorial2.cpp)
target_link_libraries(tutorial2 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial2 PRIVATE cxx_std_17)

add_executable(tutorial3 tutorial3.cpp)
target_link_libraries(tutorial3 PRIVATE FunDEM::FunDEM)
target_compile_features(tutorial3 PRIVATE cxx_std_17)

add_custom_command(
    TARGET tutorial0
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_LIST_DIR}/Gomboc.obj" "$<TARGET_FILE_DIR:tutorial0>/Gomboc.obj"
)
