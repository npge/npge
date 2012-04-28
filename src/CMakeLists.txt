include(FindBoost)
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

file(GLOB sources "*.cpp")

add_library(bloomrepeats_static STATIC ${sources})
add_executable(bloomrepeats bloomrepeats.cxx)
target_link_libraries(bloomrepeats bloomrepeats_static)

include_directories(.)
subdirs(test)

