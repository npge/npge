include_directories(../util/ ../model/ ../algo/)

file(GLOB sources "*.cxx")

foreach(source ${sources})
    get_filename_component(tool ${source} NAME_WE)
    add_executable(${tool} ${source})
    target_link_libraries(${tool} ${Boost_LIBRARIES})
    target_link_libraries(${tool} npge_util)
    target_link_libraries(${tool} npge_model)
    target_link_libraries(${tool} npge_algo)
endforeach()

file(GLOB scripts "*.npge")
set(NPGE ${CMAKE_CURRENT_BINARY_DIR}/npge)
foreach(script ${scripts})
    get_filename_component(tool ${script} NAME)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${tool}
        ${CMAKE_CURRENT_BINARY_DIR}/${tool})
endforeach()

