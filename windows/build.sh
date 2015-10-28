#!/bin/bash

set -xue

# create build environment and cross compile to static EXE file

# get path to dir windows
script=$0
if ! (echo $0 | grep -q '/windows') then
    # maybe the command is "sh .../windows/build.sh"
    script=$1
fi
WINDOWSPATH=$(cd `dirname $script`; pwd)

npge_src=$WINDOWSPATH/../

MXE_DIR=/usr/lib/mxe
MXE_CMAKE=/share/cmake/mxe-conf.cmake

if [ -z "$NOWINDOWS32" ]; then
    mkdir -p npge-build-windows32
    cd npge-build-windows32
    toolchain=$(echo $MXE_DIR/usr/i686*static/$MXE_CMAKE)
    cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain \
        -DBLAST_PLUS=1 $npge_src
    make VERBOSE=1
    cd ..
fi

if [ -z "$NOWINDOWS64" ]; then
    mkdir -p npge-build-windows64
    cd npge-build-windows64
    toolchain=$(echo $MXE_DIR/usr/x86_64*static/$MXE_CMAKE)
    cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain \
        -DBLAST_PLUS=1 $npge_src
    make VERBOSE=1
    cd ..
fi
