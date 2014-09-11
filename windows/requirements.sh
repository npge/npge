# !Debian: http://mxe.cc/#requirements

# MXE
sudo apt-get install autoconf automake bash bison bzip2 \
    cmake flex gettext git g++ intltool \
    libffi-dev libtool libltdl-dev libssl-dev \
    libxml-parser-perl make openssl patch perl \
    pkg-config scons sed unzip wget xz-utils \
    mercurial \
    wine-bin:i386 \
    curl \
    tar \
    coreutils \
    gawk \
    binutils \
    upx pandoc zip nsis

# wine-bin to generate AllProcessors.html
# curl is used to download blast TODO use wget
# tar is used to unpack blast
# coreutils for sha1sum
# gawk to convert Unix line-ends to Windows line-ends
# binutils for strip

