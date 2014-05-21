cmake_minimum_required(VERSION 2.6)
project(bloomrepeats)

enable_testing()

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(MIN_LENGTH 100 CACHE STRING "Minimum acceptable length of fragment")
set(MIN_IDENTITY 0.9 CACHE STRING "Minimum acceptable identity of block")
set(MAX_EVALUE 0.001 CACHE STRING "Maximum acceptable e-value")
set(MAX_SPREADING 0.2 CACHE STRING
    "Maximum acceptable distance spreading of fragments in block")
set(MAX_GAPS 0.2 CACHE STRING "Maximum acceptable percentage of gap positions")
set(EXPANDER_BATCH 100 CACHE STRING "Pair alignment expander batch size")
set(EXPANDER_MAX_OVERLAP 200 CACHE STRING "Neighbors overlaps when expanding")
set(ALIGNER_GAP_RANGE 15 CACHE STRING
    "Max distance from main diagonal of considered states of pair alignment")
set(ALIGNER_GAP_PENALTY 2 CACHE STRING "Gap penalty for aligner")
set(ALIGNER_MISMATCH_PENALTY 1 CACHE STRING "Mismatch penalty for aligner")
set(BOUNDARIES_MIN_DISTANCE 30 CACHE STRING
    "Min distance between fragment boundaries (overlaps resolver 2)")
set(JOINER_MAX_DIST 100 CACHE STRING "Maximum distance between joint fragments")
set(JOINER_RATIO_TO_FRAGMENT 0.5 CACHE STRING
    "Max allowed gap length to fragment length ratio")
set(JOINER_GAP_RATIO 1.5 CACHE STRING
    "Max allowed ratio of gaps' lengths (inside a block)")
set(SPLIT_REPEATS_MIN_MUTATIONS 4 CACHE STRING
    "Min number of mutations in candidate block to be splited")
set(SPLIT_REPEATS_MIN_DIAGNOSTIC_MUTATIONS 2 CACHE STRING
    "Min number of diagnostic mutations in part of splitted block")
set(MISMATCH_CHECK 1 CACHE STRING
    "Min number of equal columns after single mismatch")
set(GAP_CHECK 2 CACHE STRING
    "Min number of equal columns after single gap")
set(ALIGNED_CHECK 10 CACHE STRING
    "Min equal aligned part")
set(BLAST_PLUS 0 CACHE STRING
    "Use blast+ instead of blast")
set(BLAST_EVALUE 0.001 CACHE STRING
    "E-value filter for blast")
set(BLAST_DUST false CACHE STRING
    "E-value filter out low complexity regions")

option(BR_ASSERTS "Enable asserts" ON)

subdirs(src)
subdirs(brucella)

option(BUILD_DOC "Build documentation" ON)
if(BUILD_DOC)
    subdirs(doc)
endif()

