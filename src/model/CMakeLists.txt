include_directories(../util/)

file(GLOB sources "*.cpp")

add_library(bloomrepeats_model STATIC ${sources})
target_link_libraries(bloomrepeats_model bloomrepeats_util)
target_link_libraries(bloomrepeats_model ${Boost_LIBRARIES})

