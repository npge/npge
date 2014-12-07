sudo dpkg --add-architecture i386
echo 'deb http://ftp.debian.org/debian/ wheezy-backports main' \
    | sudo tee /etc/apt/sources.list.d/wheezy-backports.list
echo 'deb-src http://ftp.debian.org/debian/ wheezy main' \
    | sudo tee /etc/apt/sources.list.d/wheezy-src.list
sudo apt-get update

# !Debian: http://mxe.cc/#requirements

# MXE
sudo apt-get --yes install autoconf automake bash bison bzip2 \
    cmake flex gettext git g++ intltool \
    libffi-dev libtool libltdl-dev libssl-dev \
    libxml-parser-perl make openssl patch perl \
    pkg-config scons sed unzip wget xz-utils \
    autopoint gperf python ruby \
    libluajit-5.1-dev luajit \
    g++-multilib \
    wine-bin:i386 \
    curl \
    tar \
    coreutils \
    gawk \
    binutils \
    binfmt-support \
    upx-ucl pandoc zip nsis

# wine-bin to generate AllProcessors.html
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
debuild
cd ..
sudo dpkg -i nsis*.deb

