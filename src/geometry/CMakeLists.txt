add_library(FunDEM_geometry INTERFACE)
add_library(FunDEM::geometry ALIAS FunDEM_geometry)

target_sources(
    FunDEM_geometry
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_SOURCE_DIR}"
        FILES
            LSGeometry.h
)

target_include_directories(
    FunDEM_geometry
    INTERFACE
        "$<BUILD_INTERFACE:${FUNDEM_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_geometry INTERFACE cxx_std_17)
target_link_libraries(FunDEM_geometry INTERFACE FunDEM::math FunDEM::data)

set_target_properties(FunDEM_geometry PROPERTIES EXPORT_NAME geometry)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_geometry
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
