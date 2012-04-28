file(GLOB sources "*.cpp")

add_executable(bloomrepeats_test ${sources})
target_link_libraries(bloomrepeats_test bloomrepeats_static)

add_test(bloomrepeats_test bloomrepeats_test)

