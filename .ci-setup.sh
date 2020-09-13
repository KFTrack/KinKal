yum -y install make base-devel glibc-devel freetype-devel

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup mu2e
setup root v6_18_04d -f Linux64bit+3.10-2.17 -q e19:prof
scons v3_1_2 -q p383b 

rm -rf build || echo ""
mkdir build && cd build

source ../scripts/newBuild.sh prof
source setup.sh

scons -j 4

scons test
