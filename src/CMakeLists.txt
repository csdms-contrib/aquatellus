cmake_minimum_required(VERSION 2.6)

project (aquatellus)

set (BUILD_SHARED_LIBS ON)

set_source_files_properties (${aquatellus_srcs} PROPERTIES LANGUAGE CXX)

### The aquatellus library

set (aquatellus_lib_srcs
     Erf.c
     Gammq.c
     expon_int.cpp
     Gcf.c
     output.cpp
     del_lat.cpp
     flowpath.cpp
     Gser.c
     inputcontrols.cpp
     Space.cpp
     dynamics.cpp
     fluvlat2.cpp
     ei.cpp
     Gammln.c
     storage2.cpp
     aquatellus_init_func.cpp
     run_aquatellus.cpp
     aquatellus_finalize.cpp
     aquatellusreadinput.cpp
     Gammp.c)

add_library (aquatellus ${aquatellus_lib_srcs})
install (TARGETS aquatellus DESTINATION lib COMPONENT aquatellus)
install (FILES bmi.h DESTINATION include COMPONENT aquatellus RENAME aquatellus_bmi.h)

### The aquatellus executable

set (aquatellus_exe_srcs aquatellus.cpp)
add_executable (aquatellus_exe ${aquatellus_exe_srcs})
target_link_libraries (aquatellus_exe aquatellus)
install (PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/aquatellus_exe DESTINATION bin RENAME aquatellus
  COMPONENT aquatellus)


set (MAJOR_VERSION 1)
set (MINOR_VERSION 1)
set (RELEASE_VERSION 0)

set (CPACK_PACKAGE_NAME aquatellus)
set (CPACK_PACKAGE_VERSION_MAJOR ${MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${RELEASE_VERSION})
set (CPACK_SOURCE_GENERATOR "TGZ;ZIP")
set (CPACK_SOURCE_PACKAGE_FILE_NAME
  "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set (CPACK_SOURCE_IGNORE_FILES "/_build/;/_inst/;/.git/")

include (CPack)

