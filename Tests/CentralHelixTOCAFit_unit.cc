#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Tests/FitTest.hh"
// CentralHelix fit with hits and material crossings on the TOCA reference piece strategy
int main(int argc, char **argv) {
  KinKal::DVEC sigmas(0.5, 0.003, 0.00001, 3.0 , 0.004, 0.1); // expected parameter sigmas
  std::vector<std::string> arguments = {argv[0], "--Schedule", "MCAmbigFixedField.txt", "--HitRefPiece", "TOCA", "--XingRefPiece", "TOCA"};
  for(int iarg=1; iarg < argc; ++iarg) arguments.push_back(argv[iarg]);
  std::vector<char*> myargv;
  for (const auto& arg : arguments)
    myargv.push_back((char*)arg.data());
  myargv.push_back(nullptr);
  return FitTest<CentralHelix>(myargv.size()-1,myargv.data(),sigmas);
}
