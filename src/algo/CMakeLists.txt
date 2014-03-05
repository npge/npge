include_directories(../util/ ../model/)

file(GLOB sources "*.cpp")

add_library(bloomrepeats_algo STATIC ${sources})
target_link_libraries(bloomrepeats_algo boost_program_options)
target_link_libraries(bloomrepeats_algo boost_thread)
target_link_libraries(bloomrepeats_algo boost_date_time)
target_link_libraries(bloomrepeats_algo bloomrepeats_util)
target_link_libraries(bloomrepeats_algo bloomrepeats_model)

