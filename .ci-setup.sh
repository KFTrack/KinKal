echo -en 'travis_fold:start:ContainerSetup\r'

echo "set up centos7 container with CERN LCG_98"
yum -y install which make glibc-devel glibc git

source /cvmfs/sft.cern.ch/lcg/views/LCG_98/x86_64-centos7-gcc8-opt/setup.sh

echo -en 'travis_fold:end:ContainerSetup\r'

#echo -en 'travis_fold:start:Build\r'
#echo "Build"

rm -rf build || echo ""
mkdir build && cd build

cmake .. -DCMAKE_BUILD_TYPE=Release

make -j 8

#echo -en 'travis_fold:end:Build\r'

env CTEST_OUTPUT_ON_FAILURE=1 make test 

