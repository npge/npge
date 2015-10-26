#!/bin/bash

set -x

# get path to dir linux
script=$0
if ! (echo $0 | grep -q '/linux') then
    # maybe the command is "sh .../linux/build.sh"
    script=$1
fi
LINUXPATH=$(cd `dirname $script`; pwd)

npge_src=$LINUXPATH/../

mkdir -p npge-build-linux
cd npge-build-linux
cmake -DNPGE_STATIC_LINUX:BOOL=1 -DCMAKE_BUILD_TYPE=Release \
    -DBLAST_PLUS=1 -DNPGE_LUA_CMD=luajit $npge_src
make

