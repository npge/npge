sudo apt-get --yes install mercurial
hg clone https://bitbucket.org/starius/npg-explorer
cd npg-explorer
sudo ./windows/requirements.sh
./windows/build.sh
cd npge-build-windows
./windows/package.sh

