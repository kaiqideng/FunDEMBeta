foreach(index RANGE 0 3)
    add_executable(tutorial${index} tutorial${index}.cpp)
    target_link_libraries(tutorial${index} PRIVATE FunDEM::FunDEM)
endforeach()
