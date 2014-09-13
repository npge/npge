sudo apt-get --yes install mercurial &>> npge.log
hg clone https://bitbucket.org/starius/npg-explorer &>> npge.log
cd npg-explorer
sudo ./windows/requirements.sh &>> npge.log
./windows/build.sh &>> npge.log
cd npge-build-windows
./windows/package.sh &>> npge.log

