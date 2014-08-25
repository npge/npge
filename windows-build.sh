# create build environment and cross compile to static EXE file

git clone https://github.com/starius/mxe.git
cd mxe
make boost
make qt
make lua luabind
make readline
cd ..

hg clone https://bitbucket.org/starius/npg-explorer
mkdir -p npge-build
cd npge-build
toolchain=$(echo ../mxe/usr/*ming*/share/cmake/mxe-conf.cmake)
cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain ../npg-explorer/
make

