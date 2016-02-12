configure_file(${CMAKE_CURRENT_SOURCE_DIR}/package.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/package.sh)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/wizard.nsi.in
    ${CMAKE_CURRENT_BINARY_DIR}/wizard.nsi @ONLY)

set(NPGE_RC_FILENAME "npge.exe")
set(NPGE_RC_FILEDESCR " (command line tool)")
set(NPGE_RC_INTERNAL "npge")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/npge.rc.in
    ${CMAKE_CURRENT_BINARY_DIR}/npge.rc @ONLY)
set(NPGE_RC_FILENAME "qnpge.exe")
set(NPGE_RC_FILEDESCR " (GUI tool)")
set(NPGE_RC_INTERNAL "qnpge")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/npge.rc.in
    ${CMAKE_CURRENT_BINARY_DIR}/qnpge.rc @ONLY)
unset(NPGE_RC_FILENAME)
unset(NPGE_RC_FILEDESCR)
unset(NPGE_RC_INTERNAL)
