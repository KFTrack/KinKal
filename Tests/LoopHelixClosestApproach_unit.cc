#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Tests/ClosestApproachTest.hh"
int main(int argc, char **argv) {
  KinKal::DVEC pchange(0.5, 0.5, 0.5, 0.5, 0.001, 0.5); // range for parameter variation
  return ClosestApproachTest<LoopHelix>(argc,argv,pchange);
}
