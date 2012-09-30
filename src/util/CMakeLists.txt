file(GLOB sources "*.cpp")

add_library(bloomrepeats_util STATIC ${sources})
target_link_libraries(bloomrepeats_util boost_program_options)

