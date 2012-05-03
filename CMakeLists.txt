cmake_minimum_required(VERSION 2.6)
project(bloomrepeats)

enable_testing()

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

subdirs(src)
subdirs(doc)

