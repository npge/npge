if (${NPGE_STATIC_LINUX})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
    set(QT_IS_STATIC 1)
endif ()

find_package(Threads REQUIRED)
include(FindBoost)
IF (MINGW)
    SET(B_THREAD "thread_win32")
    add_definitions(-DBOOST_THREAD_USE_LIB)
ELSE (MINGW)
    SET(B_THREAD "thread")
ENDIF (MINGW)
FIND_PACKAGE(Boost 1.42
    COMPONENTS program_options filesystem system date_time ${B_THREAD}
    REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include(FindLua51)
find_package(Lua51 REQUIRED)
include_directories(${LUA_INCLUDE_DIR})
add_definitions(-DLUA_COMPAT_ALL)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR})

find_package(Luabind REQUIRED)
include_directories(${LUABIND_INCLUDE_DIRS})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp)

include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(COMMON_LIBS
 npge_algo
 npge_model
 npge_util
 ${Boost_LIBRARIES}
 ${LUABIND_LIBRARIES}
 ${LUA_LIBRARIES})
IF (MINGW)
    set(COMMON_LIBS ${COMMON_LIBS} -lws2_32)
ENDIF (MINGW)
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(COMMON_LIBS ${COMMON_LIBS} dl)
ENDIF ()

subdirs(util)
subdirs(model)
subdirs(algo)
subdirs(test)
subdirs(tool)

FIND_PACKAGE(Qt4)
option(BUILD_GUI "Build graphical user interface" ${QT4_FOUND})
if(BUILD_GUI)
    subdirs(gui)
endif()

