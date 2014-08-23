include_directories(../util/)

file(GLOB sources "*.cpp")

add_library(npge_model STATIC ${sources})
target_link_libraries(npge_model npge_util)
target_link_libraries(npge_model ${Boost_LIBRARIES})
target_link_libraries(npge_model ${LUA_LIBRARIES} luabind)

