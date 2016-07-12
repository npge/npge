#!/bin/bash

set -xue

./src/init_lua-npge.sh

brew install \
    qt4 cmake libzip libzzip boost \
    luajit luabind readline ncurses \
    pandoc tar coreutils binutils

