if(CMAKE_HOST_UNIX)

    set(BRUCELLA_TEST_TABLE 2genomes CACHE STRING
        "Table with identifiers of sequences to run tests with")

    add_custom_target(test-brucella
        ${CMAKE_CURRENT_SOURCE_DIR}/brucella.make
        TARGET=${BRUCELLA_TEST_TABLE}
        PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
    )

endif(CMAKE_HOST_UNIX)

