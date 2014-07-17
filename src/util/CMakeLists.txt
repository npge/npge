file(GLOB sources "*.cpp")

add_library(npge_util STATIC ${sources})
target_link_libraries(npge_util ${Boost_LIBRARIES})
IF (MINGW)
    target_link_libraries(npge_util ws2_32)
ENDIF (MINGW)

