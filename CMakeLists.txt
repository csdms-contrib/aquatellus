cmake_minimum_required(VERSION 2.6)

project (aquatellus)

set (BUILD_SHARED_LIBS OFF)

set (aquatellus_srcs
     aquatellus.cpp
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
     Gammp.c)

add_executable (aquatellus ${aquatellus_srcs})

install (TARGETS aquatellus DESTINATION bin)

set (MAJOR_VERSION 1)
set (MINOR_VERSION 0)
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

