#include "KinKal/LHelix.hh"
#include "UnitTests/HitTest.hh"
int main(int argc, char **argv) {
  vector<double> delpars { 0.5, 0.1, 0.5, 0.5, 0.005, 5.0};
  return HitTest<LHelix>(argc,argv,delpars);
}
