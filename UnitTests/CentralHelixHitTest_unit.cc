#include "KinKal/CentralHelix.hh"
#include "UnitTests/HitTest.hh"
int main(int argc, char **argv) {
  vector<double> delpars {0.5,0.01,0.00001,0.5,0.001,0.1};
  return HitTest<CentralHelix>(argc,argv,delpars); // small parameter changes for derivative calcs
}
