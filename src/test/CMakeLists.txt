include_directories(../util/ ../model/ ../algo/)

file(GLOB sources "*.cpp")

add_executable(npge_test ${sources})
target_link_libraries(npge_test npge_util)
target_link_libraries(npge_test npge_model)
target_link_libraries(npge_test npge_algo)
target_link_libraries(npge_test ${Boost_LIBRARIES})

add_executable(rand_seq rand_seq.cxx)

add_test(npge_test npge_test --log_level=warning)

add_executable(meta_test meta_test.cxx)
target_link_libraries(meta_test npge_util)
target_link_libraries(meta_test npge_model)
target_link_libraries(meta_test npge_algo)
target_link_libraries(meta_test ${Boost_LIBRARIES})
add_test(meta_test meta_test ${PROJECT_SOURCE_DIR}/test-script)

add_executable(lua_test lua_test.cxx)
target_link_libraries(lua_test npge_util)
target_link_libraries(lua_test npge_model)
target_link_libraries(lua_test npge_algo)
target_link_libraries(lua_test ${Boost_LIBRARIES})
target_link_libraries(lua_test ${LUA_LIBRARIES} luabind)
add_test(lua_test lua_test ${PROJECT_SOURCE_DIR}/test-lua)

add_executable(create_all create_all.cxx)
target_link_libraries(create_all npge_util)
target_link_libraries(create_all npge_model)
target_link_libraries(create_all npge_algo)
target_link_libraries(create_all ${Boost_LIBRARIES})

if(CMAKE_HOST_UNIX)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_repeats.sh.in
        ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
    add_test(find_anchors ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
endif(CMAKE_HOST_UNIX)

