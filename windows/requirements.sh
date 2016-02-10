#!/bin/bash

set -xue

# pre-built MXE packages
echo "deb http://pkg.mxe.cc/repos/apt/debian wheezy main" \
    | sudo tee /etc/apt/sources.list.d/mxeapt.list
sudo apt-key adv --keyserver x-hkp://keys.gnupg.net \
    --recv-keys D43A795B73B16ABE9643FE1AFD8FFF16DB45C6AB

sudo apt-get update

sudo apt-get --yes install \
    lua5.2 \
    curl \
    tar p7zip-full \
    coreutils \
    gawk \
    binutils \
    upx-ucl pandoc zip

# lua5.2 serves as NPGE_LUA_CMD
# curl is used to download blast TODO use wget
# tar is used to unpack blast
# coreutils for sha256sum
# gawk to convert Unix line-ends to Windows line-ends
# binutils for strip

# install MXE packages from http://pkg.mxe.cc/

sudo apt-get clean
if [ -z "$MXE_TARGET" ]; then MXE_TARGET=i686-w64-mingw32.static; fi
# x86_64 -> x86-64
MXE2_TARGET=$(echo "$MXE_TARGET" | sed 's/_/-/g')
sudo apt-get --yes install \
    mxe-$MXE2_TARGET-{qt,boost,luabind}
sudo apt-get clean
# MXE doesn't have 64bit NSIS
sudo apt-get --yes install \
    mxe-i686-w64-mingw32.static-nsis
sudo apt-get clean
