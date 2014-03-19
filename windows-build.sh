# create build environment and cross compile to static EXE file

git clone -b stable https://github.com/mxe/mxe.git
cd mxe
make boost
make qt
cd ..

hg clone https://bitbucket.org/starius/bloomrepeats
mkdir -p br-build
cd br-build
cmake -DCMAKE_TOOLCHAIN_FILE=../mxe/usr/i686-pc-mingw32.static/share/cmake/mxe-conf.cmake ../bloomrepeats/
make

