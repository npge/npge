sudo apt-get update
sudo apt-get --yes install git > npge-1.log 2>&1
git clone https://github.com/npge/npge \
    > npge-2.log 2>&1
cd npge
git pull
sudo ./linux/requirements.sh > npge-3.log 2>&1
./linux/build.sh > npge-4.log 2>&1
cd npge-build-linux
./linux/package.sh > npge-5.log 2>&1

