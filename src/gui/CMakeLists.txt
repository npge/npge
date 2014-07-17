cmake_minimum_required(VERSION 2.6)
PROJECT(qnpge)

include(CheckCXXSourceCompiles)

include_directories(../util/ ../model/ ../algo/)

INCLUDE(${QT_USE_FILE})

file(GLOB uis "*.ui")
file(GLOB headers "*.hpp")
file(GLOB sources "*.cpp")

QT4_WRAP_UI(UISrcs ${uis})
QT4_WRAP_CPP(MOCSrcs ${headers})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

ADD_EXECUTABLE(qnpge ${sources} ${MOCSrcs} ${UISrcs})
TARGET_LINK_LIBRARIES(qnpge ${QT_LIBRARIES})

target_link_libraries(qnpge ${Boost_LIBRARIES})
target_link_libraries(qnpge npge_util)
target_link_libraries(qnpge npge_model)
target_link_libraries(qnpge npge_algo)
IF (MINGW)
    target_link_libraries(qnpge lzma)
    check_cxx_source_compiles("#include <lcms2.h>\n
    int main() { }" LCMS2_FOUND)
    IF (LCMS2_FOUND)
        target_link_libraries(qnpge lcms2)
    ELSE (LCMS2_FOUND)
        target_link_libraries(qnpge lcms)
    ENDIF (LCMS2_FOUND)
    SET(CMAKE_EXE_LINKER_FLAGS -mwindows)
ENDIF (MINGW)

