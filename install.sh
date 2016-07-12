#!/bin/bash

set -xue

./src/init_lua-npge.sh

if [ "$NPGE_TARGET" != "mac" ]; then
    sudo apt-get update
    ./linux/requirements.sh
fi

./${NPGE_TARGET}/requirements.sh

if [ "$NPGE_TARGET" != "mac" ]; then
    sudo ln -s /usr/bin/luajit-* /bin/luajit
    sudo apt-get --yes install ncbi-blast+
fi
