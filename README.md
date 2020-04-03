# Kinematic Kalman fit code package

To build, you must have ROOT (see https://root.cern.ch/) installed, and ROOTSYS defined to the base of that.
You must also have the python-based SCONS build tool installed (https://scons.org/).

There are 2 build choices: debug or profile.  To build 

1. Create an empty directory to hold the github clone; lets call it KFTrack
2. cd to KFTrack, and clone the repo:
`> cd KFTrack`
`> git clone https://github.com/KFTrack/KinKal.git`
3. make a build directory; lets say build_prof (for a profile build).  cd there.
`> cd build_prof`
4. create the setup and run it
`> source ../KinKal/scripts/newBuild.sh prof`
`> source setup.sh`
5. run scons:
`> scons`

Test programs will be built in the bin directory under `build_prof`. Run them with `--help` in the `build_prof` directory to get a list of run parameters.

### Brew-specific instructions
The `setup.sh` script will export the `ROOT_INC` environment variable which points to the ROOT includes needed for the compilation. If ROOT was installed with `brew` (https://brew.sh) on macOS, this line must be changed to:

```
export ROOT_INC=${ROOTSYS}/root/include
```
