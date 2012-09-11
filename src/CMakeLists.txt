include(FindBoost)
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

file(GLOB sources "*.cpp")

add_library(bloomrepeats_static STATIC ${sources})
target_link_libraries(bloomrepeats_static boost_thread)
target_link_libraries(bloomrepeats_static boost_program_options)

add_executable(bloomrepeats bloomrepeats.cxx)
target_link_libraries(bloomrepeats bloomrepeats_static)

add_executable(rest rest.cxx)
target_link_libraries(rest bloomrepeats_static)

add_executable(add_rest add_rest.cxx)
target_link_libraries(add_rest bloomrepeats_static)

add_executable(find_anchors find_anchors.cxx)
target_link_libraries(find_anchors bloomrepeats_static)

include_directories(.)
subdirs(test)

