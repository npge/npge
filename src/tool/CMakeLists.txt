include_directories(../util/ ../model/ ../algo/)

file(GLOB sources "*.cxx")

foreach(source ${sources})
    get_filename_component(tool ${source} NAME_WE)
    add_executable(${tool} ${source})
    target_link_libraries(${tool} boost_program_options)
    target_link_libraries(${tool} bloomrepeats_util)
    target_link_libraries(${tool} bloomrepeats_model)
    target_link_libraries(${tool} bloomrepeats_algo)
endforeach()

