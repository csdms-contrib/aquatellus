cmake_minimum_required (VERSION 2.6)

project (aquatellus)

set (AQUATELLUS_EXE ${CMAKE_CURRENT_BINARY_DIR}/src/aquatellus_exe)

include (CTest)
enable_testing ()

add_subdirectory (testing/test1)
add_subdirectory (src)

