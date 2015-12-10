#!/bin/bash

set -xue

sudo apt-get update

sudo apt-get --yes install \
    g++ make cmake \
    libqt4-dev \
    libz-dev \
    libboost-iostreams-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-filesystem-dev \
    libboost-system-dev \
    libboost-thread-dev \
    luajit libluabind-dev \
    libreadline-dev libncurses5-dev \
    pandoc \
    upx-ucl tar coreutils binutils

