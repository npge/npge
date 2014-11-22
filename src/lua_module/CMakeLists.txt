
# ????* excludes "gui" (util, model, algo)
file(GLOB s1 "${PROJECT_SOURCE_DIR}/src/????*/*.cpp")
file(GLOB s2 "${PROJECT_BINARY_DIR}/src/????*/*.cpp")
file(GLOB s3 "${PROJECT_SOURCE_DIR}/src/????*/*.c")
file(GLOB s4 "${PROJECT_BINARY_DIR}/src/????*/*.c")

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/util)
include_directories(${PROJECT_SOURCE_DIR}/src/model)
include_directories(${PROJECT_SOURCE_DIR}/src/algo)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/src/util)
include_directories(${PROJECT_BINARY_DIR}/src/model)
include_directories(${PROJECT_BINARY_DIR}/src/algo)

add_library(npge_module MODULE
    ${s1};${s2};${s3};${s4})
target_link_libraries(npge_module ${COMMON_LIBS})

