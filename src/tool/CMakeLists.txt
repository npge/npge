include_directories(../util/ ../model/ ../algo/)

file(GLOB sources "*.cxx")

foreach(source ${sources})
    get_filename_component(tool ${source} NAME_WE)
    add_executable(${tool} ${source})
    target_link_libraries(${tool} boost_program_options)
    target_link_libraries(${tool} boost_filesystem boost_system)
    target_link_libraries(${tool} bloomrepeats_util)
    target_link_libraries(${tool} bloomrepeats_model)
    target_link_libraries(${tool} bloomrepeats_algo)
endforeach()

file(GLOB scripts "*.br")
set(BRMETA ${CMAKE_CURRENT_BINARY_DIR}/brmeta)
foreach(script ${scripts})
    get_filename_component(tool ${script} NAME)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${tool}
        ${CMAKE_CURRENT_BINARY_DIR}/${tool})
endforeach()

