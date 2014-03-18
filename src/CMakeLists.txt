include(FindBoost)
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp)

include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

subdirs(util)
subdirs(model)
subdirs(algo)
subdirs(test)
subdirs(tool)

option(BUILD_GUI "Build graphical user interface" ON)
if(BUILD_GUI)
    subdirs(gui)
endif()

