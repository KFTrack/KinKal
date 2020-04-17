# Kinematic Kalman filter track fit code package

  KinKal implements a kinematic Kalman filter track fit (future ref to CTD/pub).
  The primary class of KinKal is KKTrk, which shares the state describing
  the fit inputs (hits, material interactions, BField corrections, etc), and owns the result of the fit,
  and the methods for computing it.  The fit result is expressed as a piecewise kinematic covariant
  particle trajectory, providing 4-vector position and  momentum information about the particle with covariance
  as a function of physical time.

  KKTrk is templated on a simple kinematic trajectory class representing the 1-dimensional path and
  momentum of a particle traveling through empty space in a constant magnetic field, as a function of physical time.
  The piecewise kinematic trajectory fit result is expressed as a time sequence of these simple trajectory objects.
  Material effects and spatial variations of magnetic fields are modeled through changes between adjacent simple
  trajectories.  The simple kinematic trajectory class must satisfy the interface defined in the KKTrk header, including
  the 4-vector position and 4-vector momentum of the particle as a function of time.  The simple kinematic trajectory class must
  provide access to its parameterization, and the derivatives of those parameters WRT physical effects of material interactions
  and magnetic field inhomogeneity.  KinKal can be instantiated with any simple kinematic trajectory which satisfies the interface.
  KinKal provides 3 fully-implemented and tested examples
   * LHelix = low-momentum looping helix parameterized in terms of curvature radius and longitudinal wavelength
   * IPHelix = high-momentum helix parameterized in terms of inverse curvature radius and initial direction
   * KTLine = linear path with no geometric momentum interpretation

  Measurements provided to KKTrk constructor must provide a calculation of the Residual (difference between measurement
  and kinematic trajectory prediction) and the derivatives of that residual WRT the simple kinematic trajectory parameters,
  for the given choice of simple kinematic trajectory class.  KinKal provides several fully-implemented measurement classes,
  based on Residuals computed from the calculation of the TPOCA (Time and Position of Closest Approach) between the measurement
  and the simple kinematic trajectory.  The fully-functional examples include:
   * WireHit = wire chamber hit, with incomplete material an drift field properties
   * StrawHit = WireHit subclass with specific material and drift field properties
   * LightHit = rectilinear sensor with a prompt (light-based) signal, such as a scintillating crystal or plastic extrusion

  The underlying processing model used in KinKal is a progressive BLUE fit first used in the geometric track fit implementation used by the BaBar
  collaboration, described in "D.N. Brown, E.A. Charles, D.A. Roberts, The BABAR track fitting algorithm, Proceedings of CHEP 2000, Padova, Italy, 2000"

  KKTrk is constructed from a configuration object which can be shared between many instances, and a unique set of hit and
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

  The KinKal package is licensed under Adobe v2, and is hosted at [GitHub](https://github.com/KFTrack/KinKal.git)

  David N. Brown, Lawrence Berkeley National Lab

## Installation

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
