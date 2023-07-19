/*
  Original Author: S Middleton 2020
*/
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Tests/Trajectory.hh"
int main(int argc, char **argv) {
  KinKal::DVEC sigmas(0.5, 0.004, 0.5, 0.002, 0.4, 0.05); // expected parameter sigmas
  return TrajectoryTest<KinematicLine>(argc,argv,sigmas);
}
