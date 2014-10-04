cmake_minimum_required(VERSION 2.6)
project(NPG-explorer)

enable_testing()

set(VERSION_MAJOR 0)
set(VERSION_MINOR 2)
set(VERSION_PATCH 5)

set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(LOG_TO :stderr CACHE STRING "Where to write log messages")
if (WIN32)
    set(DEV_NULL NUL CACHE STRING "Null device")
else (WIN32)
    set(DEV_NULL /dev/null CACHE STRING "Null device")
endif (WIN32)
set(WORKERS -1 CACHE STRING "Number of threads (-1 = number of cores)")
set(BLOCKS_IN_GROUP 10 CACHE STRING
    "Number of blocks processing at once (BlocksJobs)")
set(TIMING 0 CACHE STRING "Log begin/end of calls and final time summary")
set(MIN_LENGTH 100 CACHE STRING "Minimum acceptable length of fragment")
set(MIN_IDENTITY 0.9 CACHE STRING "Minimum acceptable identity of block")
set(ANCHOR_SIZE 20 CACHE STRING "Anchor size (AnchorFinder)")
set(ANCHOR_FP 0.1 CACHE STRING
    "Probability of false positive in Bloom filter")
set(MAX_ANCHOR_FRAGMENTS 100000 CACHE STRING
    "Maximum number of anchors fragments to return")
set(ALIGNER "similar" CACHE STRING "Aligner implementation")
set(ALIGNER_GAP_RANGE 15 CACHE STRING
    "Max distance from main diagonal of considered states of pair alignment")
set(ALIGNER_GAP_PENALTY 2 CACHE STRING "Gap penalty for aligner")
set(ALIGNER_MISMATCH_PENALTY 1 CACHE STRING "Mismatch penalty for aligner")
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
    "Use blast+ instead of blast (may be ignored if blast tools are found in application dir)")
set(BLAST_EVALUE 0.001 CACHE STRING
    "E-value filter for blast")
set(BLAST_DUST false CACHE STRING
    "E-value filter out low complexity regions")
set(MAX_NS 3 CACHE STRING
    "Maximum number of subsequent N's in consensus")

set(MAFFT_CMD
    "{MAFFT} --quiet --retree 1 --maxiterate 1 %1% > %2%"
    CACHE STRING "Mafft command template")
set(MUSCLE_CMD
    "{MUSCLE} -quiet -in %1% -out %2% -maxiters 1 -diags1"
    CACHE STRING "Muscle command template")
set(FORMATDB_CMD
    "{FORMATDB} -l {nul} -p F -i {in} -n {bank}"
    CACHE STRING "formatdb command template (legacy blast)")
set(MAKEBLASTDB_CMD
    "{MAKEBLASTDB} -dbtype nucl -out {bank} -in {in} -logfile {nul}"
    CACHE STRING "makeblastdb command template (blast+)")
set(BLASTALL_CMD
    "{BLASTALL} -p blastn -m 8 -d {bank} -i {in} -e {evalue} -a {workers} -F {F} > {out}"
    CACHE STRING "blastall command template (legacy blast)")
set(BLASTN_CMD
    "{BLASTN} -task blastn -outfmt 6 -db {bank} -query {in} -evalue {evalue} -num_threads {workers} -dust {F} > {out}"
    CACHE STRING "blastn command template (blast+)")

set(MAX_TAIL 3 CACHE STRING "Max length of tail (MoveGaps)")
set(MAX_TAIL_TO_GAP 1.0 CACHE STRING
    "Max tail length to gap length ratio (MoveGaps)")
set(MAX_MATCHES 1000 CACHE STRING
    "Max number of matches (FragmentFinder)")

set(UPSTREAM_LENGTH 20 CACHE STRING
    "Length of upstream")

set(CONFIG0 "LOCAL_CONF" CACHE STRING "Config file 0")
set(CONFIG1 "{app_dir}/npge.conf" CACHE STRING "Config file 1")
set(CONFIG2 "/etc/npge.conf" CACHE STRING "Config file 2")
set(CONFIG3 "" CACHE STRING "Config file 3")
set(CONFIG4 "~/.npge.conf" CACHE STRING "Config file 4")
set(CONFIG5 "" CACHE STRING "Config file 5")
set(CONFIG6 "ENV" CACHE STRING "Config file 6")
set(CONFIG7 "" CACHE STRING "Config file 7")
set(CONFIG8 "LOCAL_CONF" CACHE STRING "Config file 8")
set(CONFIG9 "" CACHE STRING "Config file 9")

set(LOCAL_CONF "npge.conf" CACHE STRING
    "Local config file name")

option(NPGE_ASSERTS "Enable asserts" ON)
set(NPGE_DEBUG 0 CACHE STRING "Debug mode")

subdirs(src)
subdirs(windows)
subdirs(linux)

option(BUILD_DOC "Build documentation" ON)
if(BUILD_DOC)
    subdirs(doc)
endif()

set(PANDOC_OPTS
    "--toc -s -T NPG-explorer --css http://maxcdn.bootstrapcdn.com/bootstrap/3.2.0/css/bootstrap.min.css --include-before-body ${PROJECT_SOURCE_DIR}/doc/header.html --include-in-header ${PROJECT_SOURCE_DIR}/doc/head.html | sed -e '/#TOC/s!</a>!!g' -e 's!<a href=\"#TOC\">!!g' -e 's!NPG-explorer - </title>!NPG-explorer</title>!'")

