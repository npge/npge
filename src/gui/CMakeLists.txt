cmake_minimum_required(VERSION 2.6)
PROJECT(brgui)

include_directories(../util/ ../model/ ../algo/)

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

file(GLOB uis "*.ui")
file(GLOB headers "*.hpp")
file(GLOB sources "*.cpp")

QT4_WRAP_UI(UISrcs ${uis})
QT4_WRAP_CPP(MOCSrcs ${headers})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

ADD_EXECUTABLE(brgui ${sources} ${MOCSrcs} ${UISrcs})
TARGET_LINK_LIBRARIES(brgui ${QT_LIBRARIES})

target_link_libraries(brgui ${Boost_LIBRARIES})
target_link_libraries(brgui bloomrepeats_util)
target_link_libraries(brgui bloomrepeats_model)
target_link_libraries(brgui bloomrepeats_algo)
IF (MINGW)
    target_link_libraries(brgui lzma lcms2)
ENDIF (MINGW)

