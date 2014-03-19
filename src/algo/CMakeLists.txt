include_directories(../util/ ../model/)

file(GLOB sources "*.cpp")

add_library(bloomrepeats_algo STATIC ${sources})
target_link_libraries(bloomrepeats_algo ${Boost_LIBRARIES})
target_link_libraries(bloomrepeats_algo bloomrepeats_util)
target_link_libraries(bloomrepeats_algo bloomrepeats_model)

