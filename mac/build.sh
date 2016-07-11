#!/bin/bash

set -xue

# get path to dir mac
script=$0
if ! (echo $0 | grep -q '/mac') then
    # maybe the command is "sh .../mac/build.sh"
    script=$1
fi
MACPATH=$(cd `dirname $script`; pwd)

npge_src=$MACPATH/../

$npge_src/src/init_lua-npge.sh

mkdir -p npge-build-mac
cd npge-build-mac
cmake . -Bbuild-dir
cmake --build build-dir --config Release
