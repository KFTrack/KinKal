[![Build Status](https://travis-ci.org/KFTrack/KinKal.svg?branch=master)](https://travis-ci.org/KFTrack/KinKal)

# Kinematic Kalman filter track fit code package

  KinKal implements a kinematic Kalman filter track fit (future ref to CTD/pub).
  The primary class of KinKal is Track, which shares the state describing
  the fit inputs (hits, material interactions, BField corrections, etc), and owns the result of the fit,
  and the methods for computing it.  The fit result is expressed as a piecewise kinematic covariant
  particle trajectory, providing 4-vector position and  momentum information about the particle with covariance
  as a function of physical time.

  The Track is templated on a simple kinematic trajectory class representing the 1-dimensional path and
  momentum of a particle traveling through empty space in a constant magnetic field, as a function of physical time.
  The piecewise kinematic trajectory fit result is expressed as a time sequence of these simple trajectory objects.
  Material effects and spatial variations of magnetic fields are modeled through changes between adjacent simple
  trajectories.  The simple kinematic trajectory class must satisfy the interface defined in the Track header, including
  the 4-vector position and 4-vector momentum of the particle as a function of time.  The simple kinematic trajectory class must
  provide access to its parameterization, and the derivatives of those parameters WRT physical effects of material interactions
  and magnetic field inhomogeneity.  KinKal can be instantiated with any simple kinematic trajectory which satisfies the interface.
  KinKal provides 3 fully-implemented and tested examples
   * LoopHelix = low-momentum looping helix parameterized in terms of curvature radius and longitudinal wavelength
   * CentralHelix = high-momentum helix parameterized in terms of inverse curvature radius and initial direction
   * KLine = linear path with no geometric momentum interpretation

  Measurements provided to Track constructor must provide a calculation of the Residual (difference between measurement
  and kinematic trajectory prediction) and the derivatives of that residual WRT the simple kinematic trajectory parameters,
  for the given choice of simple kinematic trajectory class.  KinKal provides several fully-implemented measurement classes,
  based on Residuals computed from the calculation of the TPOCA (Time and Position of Closest Approach) between the measurement
  and the simple kinematic trajectory.  The fully-functional examples include:
   * WireHit = wire chamber hit, with incomplete material an drift field properties
   * StrawHit = WireHit subclass with specific material and drift field properties
   * ScintHit = rectilinear sensor with a prompt (light-based) signal, such as a scintillating crystal or plastic extrusion

  The underlying processing model used in KinKal is a progressive BLUE fit first used in the geometric track fit implementation used by the BaBar
  collaboration, described in "D.N. Brown, E.A. Charles, D.A. Roberts, The BABAR track fitting algorithm, Proceedings of CHEP 2000, Padova, Italy, 2000"

  Track is constructed from a configuration object which can be shared between many instances, and a unique set of hit and
  material interactions.  The configuration object includes a convergence schedule, which holds the parameters used in successive
  meta-iteration steps.  Each meta-iteration is described by the temperature of the simulated annealing applied to the measurement
  uncertainties, convergence critera, and flags to specify which physical effects (like material interactions) should be updated
  based on the current complete kinematic trajectory estimate.  KinKal iterates each meta-iteration to algebraic convergence,
  by re-evaluating the extended Kalman filter derivatives, holding the physical paramters of the fit fixed.
  The fit is performed on construction.

  KinKal uses the root SVector and SMatrix classes for algebraic manipulation, and GenVector classes to represent geometric and
  kinematic vectors, both part of the root Math package.  These are described on the [root website](https://root.cern.ch/root/html608/namespaceROOT_1_1Math.html)

  The KinKal package includes a number of unit test programs to verify the individual components of the fit (simple kinematic
  trajectory class, hits, etc), as well as test program with an embedded toy Monte Carlo simulation which exercises the fit,
  and verifies its performance.

  The KinKal package is licensed under Apache v2, and is hosted at [GitHub](https://github.com/KFTrack/KinKal.git)

  David N. Brown, Lawrence Berkeley National Lab

## Installation

To build, you must have ROOT (see https://root.cern.ch/) installed, and the `bin/` directory should be on your `PATH`.
You must also have the python-based SCONS build tool installed (https://scons.org/).


There are 2 build configurations: *Debug* or *Release*.  To build 

1. First, clone this repo

```bash
git clone https://github.com/KFTrack/KinKal.git
```

2. Set up a new build directory; lets say `build_profile` for a profile build (or `build_debug` for a debug build)
```bash
mkdir build_profile
cd build_profile

```

### CMake

3.1. Generate the build system, and build

```bash

cmake ../KinKal  -DCMAKE_BUILD_TYPE=[Release/Debug]

make -j <jobs to run>
```
4.1 Optionally, run unit tests

```bash
make test
```

## SCons (retired)

3.2. Setup the build environment and run `scons`

```bash
source ../KinKal/scripts/newBuild.sh prof
source setup.sh
scons -j <jobs to run>

```
4.2 Optionally run tests; NB, this doesn't work under MacOS

```bash
scons test
```

Test programs will be built in the `bin/` directory. Run them with `--help` in the `build` directory to get a list of run parameters.

### Build FAQ
#### (MacOS) Brew not working
The build tries to find ROOT with the `root-config` executable. You should ensure before building that `brew` added the ROOT `bin/` directory correctly to the `$PATH` environment variable. Sometimes re-installing the package can fix the issue.

#### Problems building against a ROOT binary or manually compiled release
You should make sure to source the `<root_location>/bin/thisroot.sh` shell script before building. This sets all the necessary environment variables needed by KinKal.
```bash

source <ROOT>/bin/thisroot.sh

```

The ROOT binaries need to be compiled with C++17 (`-std=c++17`) support.
