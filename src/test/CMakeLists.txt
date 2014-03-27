include_directories(../util/ ../model/ ../algo/)

file(GLOB sources "*.cpp")

add_executable(bloomrepeats_test ${sources})
target_link_libraries(bloomrepeats_test bloomrepeats_util)
target_link_libraries(bloomrepeats_test bloomrepeats_model)
target_link_libraries(bloomrepeats_test bloomrepeats_algo)

add_executable(rand_seq rand_seq.cxx)

add_test(bloomrepeats_test bloomrepeats_test --log_level=warning)
add_test(find_anchors_no_file
    ${PROJECT_BINARY_DIR}/src/tool/bloomrepeats.br --in-blocks nofile)

add_executable(meta_test meta_test.cxx)
target_link_libraries(meta_test bloomrepeats_util)
target_link_libraries(meta_test bloomrepeats_model)
target_link_libraries(meta_test bloomrepeats_algo)
target_link_libraries(meta_test ${Boost_LIBRARIES})
add_test(meta_test meta_test ${PROJECT_SOURCE_DIR}/test-script)

add_executable(create_all create_all.cxx)
target_link_libraries(create_all bloomrepeats_util)
target_link_libraries(create_all bloomrepeats_model)
target_link_libraries(create_all bloomrepeats_algo)
target_link_libraries(create_all ${Boost_LIBRARIES})

if(CMAKE_HOST_UNIX)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_repeats.sh.in
        ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
    add_test(find_anchors ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
endif(CMAKE_HOST_UNIX)

