#include "KinKal/KTLine.hh"
#include "UnitTests/HitTest.hh"
int main(int argc, char **argv) {
  vector<double> delpars { 1,0.5,1,0.5,0.005, 5.0}; //d0_=0,phi0_=1,z0_=2,cost_=3,t0_=4
  return HitTest<KTLine>(argc,argv, delpars);
}
