cd npge-build-windows32
export BLASTARCH='ia32'
export NPGEARCH='32'
./windows/package.sh
cd ..

cd npge-build-windows64
export BLASTARCH='x64'
export NPGEARCH='64'
./windows/package.sh
cd ..
