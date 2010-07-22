cmake_minimum_required(VERSION 2.6)

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

