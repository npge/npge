if(CMAKE_HOST_UNIX)

set(BRUCELLA_TEST_TABLE 2genomes CACHE STRING
    "Table with identifiers of sequences to run tests with")

add_custom_target(test-brucella-targets
    ${CMAKE_CURRENT_SOURCE_DIR}/make_target_list.sh
        ${BRUCELLA_TEST_TABLE} ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}
)

add_custom_target(test-brucella
    ${CMAKE_CURRENT_SOURCE_DIR}/brucella.make
    `cat ${CMAKE_CURRENT_BINARY_DIR}/targets`
        TARGET=${BRUCELLA_TEST_TABLE}
        PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
)

add_dependencies(test-brucella test-brucella-targets)

endif(CMAKE_HOST_UNIX)

