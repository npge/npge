#!/bin/bash

set -x

if [ -z "$NO_WINDOWS32" ]; then
    cd npge-build-windows32
    export BLASTARCH='ia32'
    export NPGEARCH='32'
    ./windows/package.sh
    cd ..
fi

if [ -z "$NO_WINDOWS64" ]; then
    cd npge-build-windows64
    export BLASTARCH='x64'
    export NPGEARCH='64'
    ./windows/package.sh
    cd ..
fi
