sudo dpkg --add-architecture amd64
sudo dpkg --add-architecture i386
sudo apt-get update

sudo apt-get install \
    gcc-multilib g++-multilib \
    g++ make cmake cmake-curses-gui \
    libqt4-dev:amd64 libboost-all-dev:amd64 \
    libqt4-dev:i386 libboost-all-dev:i386 \
    liblua5.1-dev:amd64 libluabind-dev:amd64 \
    liblua5.1-dev:i386 libluabind-dev:i386 \
    libreadline-dev:amd64 libncurses5-dev:amd64 \
    libreadline-dev:i386 libncurses5-dev:i386 \
    pandoc doxygen

