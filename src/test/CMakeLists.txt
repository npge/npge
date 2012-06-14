file(GLOB sources "*.cpp")

add_executable(bloomrepeats_test ${sources})
target_link_libraries(bloomrepeats_test bloomrepeats_static)

add_executable(rand_seq rand_seq.cxx)

add_test(bloomrepeats_test bloomrepeats_test)

if(CMAKE_HOST_UNIX)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_repeats.sh.in
        ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
    add_test(find_anchors ${CMAKE_CURRENT_BINARY_DIR}/test_repeats.sh)
endif(CMAKE_HOST_UNIX)

