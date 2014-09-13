sudo apt-get --yes install \
    g++ make cmake cmake-curses-gui \
    libqt4-dev libboost-all-dev \
    liblua5.1-dev libluabind-dev \
    libreadline-dev libncurses5-dev \
    pandoc doxygen

sudo apt-get --yes install dpkg-dev devscripts
sudo apt-get --yes build-dep qt4-x11
apt-get source qt4-x11
cd qt4-x11-*
QO="-static -qt-zlib -no-gif -qt-libpng -qt-libmng -qt-libjpeg"
sed "s@./configure@./configure $QO@" -i debian/rules
debuild
cd ..
sudo dpkg -i *.deb
# ^^ TODO names of packages

