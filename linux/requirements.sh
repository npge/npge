#!/bin/bash

set -x

echo 'deb http://ftp.debian.org/debian/ wheezy-backports main' \
    | sudo tee /etc/apt/sources.list.d/wheezy-backports.list

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

