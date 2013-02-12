if(CMAKE_HOST_UNIX)

    set(BRUCELLA_TEST_TABLE 2genomes CACHE STRING
        "Table with identifiers of sequences to run tests with")

    add_custom_target(test-brucella
        ${CMAKE_CURRENT_SOURCE_DIR}/brucella.sh
        ${BRUCELLA_TEST_TABLE} ${PROJECT_BINARY_DIR}/src/tool/
    )

endif(CMAKE_HOST_UNIX)

