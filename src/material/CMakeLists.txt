add_library(FunDEM_material INTERFACE)
add_library(FunDEM::material ALIAS FunDEM_material)

target_sources(
    FunDEM_material
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_SOURCE_DIR}"
        FILES
            LSMaterial.h
            material.h
)

target_include_directories(
    FunDEM_material
    INTERFACE
        "$<BUILD_INTERFACE:${FUNDEM_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_material INTERFACE cxx_std_17)
target_link_libraries(FunDEM_material INTERFACE FunDEM::math FunDEM::data)

set_target_properties(FunDEM_material PROPERTIES EXPORT_NAME material)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_material
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
