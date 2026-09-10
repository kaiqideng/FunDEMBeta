add_library(FunDEM_math INTERFACE)
add_library(FunDEM::math ALIAS FunDEM_math)

target_sources(
    FunDEM_math
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES
            Constants.h
            Indexing.h
            Vector3.h
            Matrix3.h
            Quaternion.h
)

target_include_directories(
    FunDEM_math
    INTERFACE
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_math INTERFACE cxx_std_17)
set_target_properties(FunDEM_math PROPERTIES EXPORT_NAME math)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_math
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
