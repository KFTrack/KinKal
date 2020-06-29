yum -y install make base-devel glibc-devel freetype-devel

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup mu2e
setup root v6_18_04d -f Linux64bit+3.10-2.17 -q e19:prof
setup cmake v3_9_0

rm -rf build || echo ""
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

make -j 8

make test
