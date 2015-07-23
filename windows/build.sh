# create build environment and cross compile to static EXE file

# get path to dir windows
script=$0
if ! (echo $0 | grep -q '/windows') then
    # maybe the command is "sh .../windows/build.sh"
    script=$1
fi
WINDOWSPATH=$(cd `dirname $script`; pwd)

MXE_TARGETS='i686-w64-mingw32.static x86_64-w64-mingw32.static'
# replace with https://github.com/mxe/mxe.git
# when https://github.com/mxe/mxe/pull/579
# is merged
git clone https://github.com/starius/mxe.git
cd mxe
make boost qt luajit luabind readline zlib \
    "MXE_TARGETS=$MXE_TARGETS" LUA=luajit
cd ..

npge_src=$WINDOWSPATH/../

mkdir -p npge-build-windows32
cd npge-build-windows32
toolchain=$(echo ../mxe/usr/i686*static/share/cmake/mxe-conf.cmake)
cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain \
    -DBLAST_PLUS=1 $npge_src
make
cd ..

mkdir -p npge-build-windows64
cd npge-build-windows64
toolchain=$(echo ../mxe/usr/x86_64*static/share/cmake/mxe-conf.cmake)
cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain \
    -DBLAST_PLUS=1 $npge_src
make
cd ..

