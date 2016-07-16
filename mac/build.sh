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
cmake $npge_src -Bbuild-dir -DNPGE_LUA_CMD=lua5.1 -DUSE_LUAJIT=OFF
cmake --build build-dir --config Release
