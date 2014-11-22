
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

IF (MINGW)
    SET(B_THREAD "thread_win32")
ELSE (MINGW)
    SET(B_THREAD "thread")
ENDIF (MINGW)
FIND_PACKAGE(Boost 1.42
    COMPONENTS program_options filesystem system date_time ${B_THREAD}
    REQUIRED)

add_library(npge_module MODULE
    ${s1};${s2};${s3};${s4})
target_link_libraries(npge_module
    -lluabind
    -lboost_program_options
    -lboost_filesystem
    -lboost_system
    -lboost_date_time
    -lboost_${B_THREAD})

