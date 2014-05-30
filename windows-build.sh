# create build environment and cross compile to static EXE file

git clone -b stable https://github.com/mxe/mxe.git
cd mxe
make boost
make qt
cd ..

hg clone https://bitbucket.org/starius/npg-explorer
mkdir -p npge-build
cd npge-build
toolchain=$(echo ../mxe/usr/*ming*/share/cmake/mxe-conf.cmake)
cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain ../npg-explorer/
make

