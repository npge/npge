include(FindBoost)
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(.)

subdirs(util)
subdirs(model)
subdirs(algo)
subdirs(test)
subdirs(tool)

