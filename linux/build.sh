# get path to dir linux
script=$0
if ! (echo $0 | grep -q '/linux') then
    # maybe the command is "sh .../linux/build.sh"
    script=$1
fi
LINUXPATH=`dirname $script`

npge_src=$LINUXPATH/../

mkdir -p npge-build-linux
cd npge-build-linux
cmake -DNPGE_STATIC_LINUX:BOOL=1 -DCMAKE_BUILD_TYPE=Release \
    $npge_src
make

