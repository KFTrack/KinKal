#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Tests/Trajectory.hh"
int main(int argc, char **argv){
  KinKal::DVEC sigmas(0.5, 0.003, 0.00001, 3.0 , 0.004, 0.1); // expected parameter sigmas
  return TrajectoryTest<CentralHelix>(argc,argv,sigmas);
}
