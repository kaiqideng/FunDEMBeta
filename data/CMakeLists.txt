add_library(FunDEM_dataLayout INTERFACE)
add_library(FunDEM::dataLayout ALIAS FunDEM_dataLayout)

target_sources(
    FunDEM_dataLayout
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES
            ContainerObject.h
            DeviceLayout.h
)

target_include_directories(
    FunDEM_dataLayout
    INTERFACE
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_dataLayout INTERFACE cxx_std_17)
set_target_properties(FunDEM_dataLayout PROPERTIES EXPORT_NAME dataLayout)

add_library(FunDEM_data INTERFACE)
add_library(FunDEM::data ALIAS FunDEM_data)

target_sources(
    FunDEM_data
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES
            CudaTypes.h
            HostAoSDeviceSoA.h
)

target_include_directories(
    FunDEM_data
    INTERFACE
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_data INTERFACE cxx_std_17)
target_compile_definitions(FunDEM_data INTERFACE FUNDEM_HAS_CUDA=$<BOOL:${FUNDEM_HAS_CUDA}>)
target_link_libraries(FunDEM_data INTERFACE FunDEM::dataLayout)
if(FUNDEM_HAS_CUDA)
    target_link_libraries(FunDEM_data INTERFACE CUDA::cudart)
endif()

set_target_properties(FunDEM_data PROPERTIES EXPORT_NAME data)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_dataLayout
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )

    install(
        TARGETS FunDEM_data
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()

add_library(FunDEM_dataIO vtuWriter.cpp datWriter.cpp myLSObject.cpp)
add_library(FunDEM::dataIO ALIAS FunDEM_dataIO)

target_sources(
    FunDEM_dataIO
    PUBLIC
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS "${FUNDEM_ROOT}"
        FILES
            vtuWriter.h
            datWriter.h
            myLSObject.h
)

target_include_directories(
    FunDEM_dataIO
    PUBLIC
        "$<BUILD_INTERFACE:${FUNDEM_ROOT}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FunDEM>"
)

target_compile_features(FunDEM_dataIO PUBLIC cxx_std_17)
target_link_libraries(FunDEM_dataIO PUBLIC FunDEM::math FunDEM::data)
if(FUNDEM_HAS_OPENMP)
    target_link_libraries(FunDEM_dataIO PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(FunDEM_dataIO PROPERTIES EXPORT_NAME dataIO)

if(FUNDEM_INSTALL)
    install(
        TARGETS FunDEM_dataIO
        EXPORT FunDEMTargets
        FILE_SET public_headers
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/FunDEM"
    )
endif()
