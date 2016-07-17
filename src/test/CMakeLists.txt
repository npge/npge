include_directories(../util/ ../model/ ../algo/)

set(exe_suffix "")
IF (MINGW)
    set(exe_suffix ".exe")
endif ()

file(GLOB sources "*.cpp")

add_executable(npge_test ${sources})
target_link_libraries(npge_test ${COMMON_LIBS})

add_executable(rand_seq rand_seq.cxx)

add_test(npge_test npge_test${exe_suffix} --log_level=warning)

# TODO (for zer0main): install blast binary for Mac in Travis
# and run meta_test on OSX
IF (NOT APPLE)
    add_executable(meta_test meta_test.cxx)
    target_link_libraries(meta_test ${COMMON_LIBS})
    add_test(meta_test
        meta_test${exe_suffix} ${PROJECT_SOURCE_DIR}/test-script)
endif ()

add_executable(lua_test lua_test.cxx)
target_link_libraries(lua_test ${COMMON_LIBS})
add_test(lua_test
    lua_test${exe_suffix} ${PROJECT_SOURCE_DIR}/test-lua)

add_custom_target(test-create-all
    ${PROJECT_BINARY_DIR}/src/tool/npge${exe_suffix}
    ${CMAKE_CURRENT_SOURCE_DIR}/create_all.npge)

if(CMAKE_HOST_UNIX)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_repeats.sh.in
        ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
    add_test(find_anchors ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
endif(CMAKE_HOST_UNIX)

