find_package(Threads REQUIRED)
include(FindBoost)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
FIND_PACKAGE(Boost 1.42 COMPONENTS program_options REQUIRED)
FIND_PACKAGE(Boost 1.42 COMPONENTS filesystem REQUIRED)
FIND_PACKAGE(Boost 1.42 COMPONENTS system REQUIRED)
IF (MINGW)
    FIND_PACKAGE(Boost 1.42 COMPONENTS thread_win32 REQUIRED)
    add_definitions(-DBOOST_THREAD_USE_LIB)
ELSE (MINGW)
    FIND_PACKAGE(Boost 1.42 COMPONENTS thread REQUIRED)
ENDIF (MINGW)
FIND_PACKAGE(Boost 1.42 COMPONENTS date_time REQUIRED)

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

