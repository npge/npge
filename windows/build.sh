# create build environment and cross compile to static EXE file

# get path to dir windows
script=$0
if ! (echo $0 | grep -q '/windows') then
    # maybe the command is "sh .../windows/build.sh"
    script=$1
fi
WINDOWSPATH=$(cd `dirname $script`; pwd)

git clone https://github.com/mxe/mxe.git
cd mxe
make boost qt luajit luabind readline LUA=luajit
cd ..

npge_src=$WINDOWSPATH/../

mkdir -p npge-build-windows
cd npge-build-windows
toolchain=$(echo ../mxe/usr/*ming*/share/cmake/mxe-conf.cmake)
cmake -DCMAKE_TOOLCHAIN_FILE=$toolchain \
    -DBLAST_PLUS=1 $npge_src
make

