#!/bin/bash

set -x

# pre-built MXE packages
echo "deb http://mxe.redjohn.tk/repos/apt/debian wheezy main" \
    | sudo tee /etc/apt/sources.list.d/mxeapt.list
cat windows/mxeapt.gpg | sudo apt-key add -

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
# coreutils for sha1sum
# gawk to convert Unix line-ends to Windows line-ends
# binutils for strip

# install MXE packages from http://mxe.redjohn.tk/
sudo apt-get clean
if [ -z "$NOWINDOWS32" ]; then
    sudo apt-get --yes install \
        mxe-i686-w64-mingw32.static-{qt,boost,luabind,nsis}
    sudo apt-get clean
fi
if [ -z "$NOWINDOWS64" ]; then
    sudo apt-get --yes install \
        mxe-x86-64-w64-mingw32.static-{qt,boost,luabind,nsis}
    sudo apt-get clean
    # MXE doesn't have 64bit NSIS
    sudo apt-get --yes install \
        mxe-i686-w64-mingw32.static-nsis
    sudo apt-get clean
fi
