set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR})
FIND_PACKAGE(Readline)

include(CheckCSourceCompiles)

CHECK_C_SOURCE_COMPILES("#include <stdio.h>\n int main() {
    vasprintf(0, 0, 0); }" HAVE_ASPRINTF)

CHECK_C_SOURCE_COMPILES("#include <wordexp.h>\n int main;"
    HAVE_WORDEXP)
CHECK_C_SOURCE_COMPILES("#include <sys/ioctl.h>\n int main;"
    HAVE_IOCTL)
CHECK_C_SOURCE_COMPILES("#include <term.h>\n int main;"
    HAVE_TERM_H)

if (READLINE_FOUND AND HAVE_TERM_H)
    set(USE_LUAPROMPT TRUE)
endif()

include_directories(../util/ ../model/ ../algo/)

file(GLOB sources "*.cxx")

foreach(source ${sources})
    get_filename_component(tool ${source} NAME_WE)
    set(this_sources ${source})
    if (USE_LUAPROMPT AND ${tool} STREQUAL "lnpge")
        set(this_sources ${this_sources};luaprompt/prompt.c)
    endif()
    add_executable(${tool} ${this_sources})
    if (USE_LUAPROMPT AND ${tool} STREQUAL "lnpge")
        include_directories(${READLINE_INCLUDE_DIR})
        include_directories(${CMAKE_CURRENT_SOURCE_DIR}/luaprompt)
        add_definitions(-DLUAPROMPT)
        add_definitions(-D_GNU_SOURCE
            -DHAVE_LIBREADLINE -DHAVE_READLINE_READLINE_H
            -DHAVE_READLINE_HISTORY -DHAVE_READLINE_HISTORY_H)
        if (HAVE_ASPRINTF)
            add_definitions(-DHAVE_ASPRINTF)
        endif()
        if (HAVE_IOCTL)
            add_definitions(-DHAVE_IOCTL)
        endif()
        if (HAVE_WORDEXP)
            add_definitions(-DHAVE_WORDEXP)
        endif()
        target_link_libraries(${tool} -lreadline -lhistory
            -lcurses)
    endif()
    target_link_libraries(${tool} ${Boost_LIBRARIES})
    target_link_libraries(${tool} ${LUA_LIBRARIES} luabind)
    target_link_libraries(${tool} npge_util)
    target_link_libraries(${tool} npge_model)
    target_link_libraries(${tool} npge_algo)
endforeach()

