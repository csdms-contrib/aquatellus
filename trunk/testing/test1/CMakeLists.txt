add_test (test1 ${CMAKE_CURRENT_BINARY_DIR}/test1.sh)

configure_file (${CMAKE_CURRENT_SOURCE_DIR}/test1.sh.in
                ${CMAKE_CURRENT_BINARY_DIR}/test1.sh)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/aquatellusinput.rtf.in
                ${CMAKE_CURRENT_BINARY_DIR}/aquatellusinput.rtf)

