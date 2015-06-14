sudo apt-get update
sudo apt-get --yes install git > npge-1.log 2>&1
git clone https://github.com/npge/npge \
    > npge-2.log 2>&1
cd npge
git pull
sudo ./windows/requirements.sh > npge-3.log 2>&1
./windows/build.sh > npge-4.log 2>&1
cd npge-build-windows32
export BLASTARCH='ia32'
export NPGEARCH='32'
./windows/package.sh > npge-5-32.log 2>&1
cd ..
cd npge-build-windows64
export BLASTARCH='x64'
export NPGEARCH='64'
./windows/package.sh > npge-5-64.log 2>&1
