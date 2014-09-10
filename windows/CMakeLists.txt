configure_file(${CMAKE_CURRENT_SOURCE_DIR}/package.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/package.sh)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/wizard.nsi.in
    ${CMAKE_CURRENT_BINARY_DIR}/wizard.nsi @ONLY)

