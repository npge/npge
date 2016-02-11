#!/bin/bash

set -xue

MXE_TARGET="${MXE_TARGET-i686-w64-mingw32.static}"

if [ "$MXE_TARGET" = "i686-w64-mingw32.static" ]; then
    cd npge-build-windows32
    export BLASTARCH='ia32'
    export NPGEARCH='32'
    ./windows/package.sh
    cd ..
fi

if [ "$MXE_TARGET" = "x86_64-w64-mingw32.static" ]; then
    cd npge-build-windows64
    export BLASTARCH='x64'
    export NPGEARCH='64'
    ./windows/package.sh
    cd ..
fi
