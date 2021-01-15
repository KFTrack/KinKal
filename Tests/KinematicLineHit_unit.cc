#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Tests/HitTest.hh"
int main(int argc, char **argv) {
  vector<double> delpars { 1.0,0.001,1,0.001,0.1,1.0}; //d0_=0,phi0_=1,z0_=2,cost_=3,t0_=4, mom_=5
  return HitTest<KinematicLine>(argc,argv, delpars);
}
