if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    mkdir build && cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    
    make -j 8
    env CTEST_OUTPUT_ON_FAILURE=1 make test
else
    # set up SL7 docker container
    cd .. 

    docker run --name KinKalCI -v /cvmfs:/cvmfs:ro,shared -v "$(pwd)"/KinKal:/KinKal -it -d scientificlinux/sl:7
    docker exec -ti KinKalCI /bin/bash -c "cd /KinKal && source .ci-setup.sh"
fi

