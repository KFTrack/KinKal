if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    mkdir build && cd build

    source ../../scripts/newBuild.sh prof 
    source setup.sh

    scons -j 4
    scons test

else

    cd .. 

    docker run --name KinKalCI -v /cvmfs:/cvmfs:ro,shared -v "$(pwd)"/KinKal:/KinKal -it -d scientificlinux/sl:7
    docker exec -ti KinKalCI /bin/bash -c "cd /KinKal && source .ci-setup.sh"

fi

