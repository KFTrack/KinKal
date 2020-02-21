# Kinematic Kalman fit code package

To build, you must have root (see https://root.cern.ch/) installed, and ROOTSYS defined to the base of that.
You must also have the python-basee SCONS build tool installed (https://scons.org/).

There are 2 build choices: debug or profile.  To build 

1. Create an empty directory to hold the github clone; lets call it KFTrack
2. cd to KFTrack, and clone the repo:
> cd KFTrack
> git clone https://github.com/KFTrack/KinKal.git
3. make a build directory; lets say build_prof (for a profile build).  cd there.
> cd build_prof
4. create the scripts
> source ../KinKal/scripts/newBuild.sh prof
5. run scons:
> scons

Test programs will be built in the bin directory under build_prof.  run them with --help to get a list of run parameters
