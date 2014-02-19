cmake_minimum_required(VERSION 2.6)
project(bloomrepeats)

enable_testing()

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(MIN_LENGTH 100 CACHE STRING "Minimal acceptable length of fragment")
set(MIN_IDENTITY 0.9 CACHE STRING "Minimal acceptable identity of block")
set(MAX_EVALUE 0.001 CACHE STRING "Maximal acceptable e-value")
set(MAX_SPREADING 0.2 CACHE STRING
    "Maximal acceptable distance spreading of fragments in block")
set(MAX_GAPS 0.2 CACHE STRING "Maximal acceptable percentage of gap positions")
set(EXPANDER_BATCH 100 CACHE STRING "Pair alignment expander batch size")
set(EXPANDER_MAX_OVERLAP 200 CACHE STRING "Neighbors overlaps when expanding")
set(ALIGNER_GAP_RANGE 5 CACHE STRING
    "Max distance from main diagonal of considered states of pair alignment")
set(ALIGNER_GAP_PENALTY 2 CACHE STRING "Gap penalty for aligner")
set(BOUNDARIES_MIN_DISTANCE 100 CACHE STRING
    "Min distance between fragment boundaries (overlaps resolver 2)")
set(JOINER_MAX_DIST 100 CACHE STRING "Maximal distance between joint fragments")
set(JOINER_RATIO_TO_FRAGMENT 0.5 CACHE STRING
    "Max allowed gap length to fragment length ratio")
set(JOINER_GAP_RATIO 1.5 CACHE STRING
    "Max allowed ratio of gaps' lengths (inside a block)")

subdirs(src)
subdirs(brucella)

option(BUILD_DOC "Build documentation" ON)
if(BUILD_DOC)
    subdirs(doc)
endif()

