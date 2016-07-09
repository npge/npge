#!/bin/bash

set -xue

./src/init_lua-npge.sh

brew install \
    qt4 cmake libz boost \
    luajit luabind readline ncurses5 \
    pandoc tar coreutils binutils

