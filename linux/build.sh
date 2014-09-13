# get path to dir linux
script=$0
if ! (echo $0 | grep -q '/linux') then
    # maybe the command is "sh .../linux/build.sh"
    script=$1
fi
LINUXPATH=`dirname $script`

npge_src=$LINUXPATH/../

mkdir -p npge-build-linux64
cd npge-build-linux64
cmake -DNPGE_STATIC_LINUX:BOOL=1 -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=$LINUXPATH/toolchain-linux64.cmake \
    $npge_src
make
cd ..

mkdir -p npge-build-linux32
cd npge-build-linux32
cmake -DNPGE_STATIC_LINUX:BOOL=1 -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE=$LINUXPATH/toolchain-linux32.cmake \
    $npge_src
make
cd ..

