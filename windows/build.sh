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

$npge_src/src/init_lua-npge.sh

MXE_DIR=/usr/lib/mxe

MXE_TARGET="${MXE_TARGET-i686-w64-mingw32.static}"

if [ "$MXE_TARGET" = "i686-w64-mingw32.static" ]; then
    mkdir -p npge-build-windows32
    cd npge-build-windows32
    $MXE_DIR/usr/bin/i686-w64-mingw32.static-cmake \
        -DUSE_LUAJIT=OFF \
        -DBLAST_PLUS=1 $npge_src
    make VERBOSE=1
    cd ..
fi

if [ "$MXE_TARGET" = "x86_64-w64-mingw32.static" ]; then
    mkdir -p npge-build-windows64
    cd npge-build-windows64
    $MXE_DIR/usr/bin/x86_64-w64-mingw32.static-cmake \
        -DUSE_LUAJIT=OFF \
        -DBLAST_PLUS=1 $npge_src
    make VERBOSE=1
    cd ..
fi
