#!/bin/bash

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
    binfmt-support \
    upx-ucl pandoc zip nsis

# lua5.2 serves as NPGE_LUA_CMD
# curl is used to download blast TODO use wget
# tar is used to unpack blast
# coreutils for sha1sum
# gawk to convert Unix line-ends to Windows line-ends
# binutils for strip
# binfmt-support allows to run EXE files transparently

# Rebuild NSIS with larger NSIS_MAX_STRLEN
# https://community.openvpn.net/openvpn/wiki/%3ABuildingMakeNSIS

sudo apt-get --yes install dpkg-dev devscripts
sudo apt-get --yes build-dep nsis
apt-get source nsis
cd nsis-*
sed '/^SCONSOPTS_\(NSIS\|COMMON\)/s@SKIPUTILS@NSIS_MAX_STRLEN=8192 SKIPUTILS@' \
    -i debian/rules
DEB_BUILD_OPTIONS=nocheck debuild
cd ..
sudo dpkg -i nsis*.deb

# install MXE packages from http://mxe.redjohn.tk/
sudo apt-get install \
    mxe-{i686,x86-64}-w64-mingw32.static-{qt,boost,luabind}
