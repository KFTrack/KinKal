#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Tests/Trajectory.hh"

int main(int argc, char **argv) {
  KinKal::DVEC sigmas(0.5, 0.5, 0.5, 0.5, 0.005, 0.5); // expected parameter sigmas
  return TrajectoryTest<LoopHelix>(argc,argv,sigmas);
}
