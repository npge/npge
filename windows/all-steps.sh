sudo apt-get --yes install mercurial &> npge-1.log
hg clone https://bitbucket.org/starius/npg-explorer &> npge-2.log
cd npg-explorer
hg pull -u
sudo ./windows/requirements.sh &> npge-3.log
./windows/build.sh &> npge-4.log
cd npge-build-windows
./windows/package.sh &> npge-5.log

