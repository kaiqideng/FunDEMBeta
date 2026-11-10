# Download the author's unmodified initial coordinates without vendoring an
# unlicensed third-party dataset or silently changing an existing local file.
cmake_minimum_required(VERSION 3.24)

if(NOT DEFINED DESTINATION)
    set(DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/build/validation/ondemShearFlow_input/GranularShearFlow.data")
endif()
get_filename_component(DESTINATION "${DESTINATION}" ABSOLUTE)
get_filename_component(destination_directory "${DESTINATION}" DIRECTORY)
set(expected_sha256 "947fec877c67e9ba7e8d1d46b66f47f66df6d713e0b33014556f967ce7ca099d")
set(source_url "https://raw.githubusercontent.com/dalilavescovi/GranularShearFlows/56b7cfdfd774cf5ae057a1c5c90e5191a0deb97b/InitialGeometricalConfiguration/GranularShearFlow.data")

if(EXISTS "${DESTINATION}")
    file(SHA256 "${DESTINATION}" existing_sha256)
    if(NOT existing_sha256 STREQUAL expected_sha256)
        message(FATAL_ERROR "A different file already exists at ${DESTINATION}; it was not overwritten.")
    endif()
    message(STATUS "Verified existing official initial configuration: ${DESTINATION}")
    return()
endif()

file(MAKE_DIRECTORY "${destination_directory}")
file(DOWNLOAD "${source_url}" "${DESTINATION}.download"
    TLS_VERIFY ON
    EXPECTED_HASH "SHA256=${expected_sha256}"
    STATUS download_status
    TIMEOUT 60
    SHOW_PROGRESS)
list(GET download_status 0 download_code)
if(NOT download_code EQUAL 0)
    message(FATAL_ERROR "Initial-configuration download failed: ${download_status}")
endif()
file(RENAME "${DESTINATION}.download" "${DESTINATION}" NO_REPLACE RESULT rename_result)
if(NOT rename_result STREQUAL "0")
    message(FATAL_ERROR "Could not install downloaded data without replacement: ${rename_result}")
endif()
message(STATUS "Fetched and SHA256-verified official initial configuration: ${DESTINATION}")
