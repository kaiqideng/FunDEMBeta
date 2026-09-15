add_library(FunDEM_particle pointMass.cpp rigidBody.cpp)
add_library(FunDEM::particle ALIAS FunDEM_particle)

target_sources(
    FunDEM_particle
    PUBLIC
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_SOURCE_DIR}"
        FILES
            pointMass.h
            rigidBody.h
            particle.h
            LSParticle.h
            SPHParticle.h
            SPHJet.h
            virtualParticle.h
            spatialGrid.h
)

target_include_directories(
    FunDEM_particle
    PUBLIC
        "$<BUILD_INTERFACE:${FUNDEM_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_particle PUBLIC cxx_std_17)
target_link_libraries(
    FunDEM_particle
    PUBLIC
        FunDEM::math
        FunDEM::data
        FunDEM::material
        FunDEM::geometry
        FunDEM::execution
)

set_target_properties(FunDEM_particle PROPERTIES EXPORT_NAME particle)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_particle
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
