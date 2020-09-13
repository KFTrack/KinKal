if [ $TRAVIS_OS_NAME = 'osx' ]; then
  wget https://root.cern/download/root_v6.22.02.macosx64-10.13-clang100.pkg
  installer -pkg root_v6.22.02.macosx64-10.13-clang100.pkg -target CurrentUserHomeDirectory
  source ~/Applications/root_v6.22.02/bin/thisroot.sh  
else
  wget https://root.cern/download/root_v6.22.02.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
  tar xzvf root_v6.22.02.Linux-ubuntu20-x86_64-gcc9.3.tar.gz
  source ~/root/bin/thisroot.sh

fi

python -m pip install --user scons
