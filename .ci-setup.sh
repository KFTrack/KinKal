

#yum -y install http://ftp.scientificlinux.org/linux/scientific/7x/contexts/fermilab/x86_64/yum-conf-context-fermilab-1.0-6.el7.noarch.rpm
#sleep 1
#yum -y upgrade
yum -y install make base-devel glibc-devel freetype-devel # redhat-lsb-core

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup mu2e
setup root v6_18_04d -f Linux64bit+3.10-2.17 -q e19:prof
setup cmake v3_9_0

rm -rf build || echo "build exists"
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

make -j 8

make test
