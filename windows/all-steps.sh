sudo apt-get --yes install mercurial > npge-1.log 2>&1
hg clone https://bitbucket.org/starius/npg-explorer \
    > npge-2.log 2>&1
cd npg-explorer
hg pull -u
sudo ./windows/requirements.sh > npge-3.log 2>&1
./windows/build.sh > npge-4.log 2>&1
cd npge-build-windows
./windows/package.sh > npge-5.log 2>&1

