#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Tests/FitTest.hh"
int main(int argc, char **argv) {
  KinKal::DVEC sigmas(0.5, 0.5, 0.5, 0.5, 0.002, 0.5); // expected parameter sigmas
  if(argc == 1){
    cout << "Testing gradient field, correction 2" << endl;
    std::vector<std::string> arguments;
    arguments.push_back(argv[0]);
    arguments.push_back("--Bgrad");
    arguments.push_back("-0.036"); // mu2e-like field gradient
    arguments.push_back("--bfcor");
    arguments.push_back("2"); // local field correction (BField rotation)
    std::vector<char*> myargv;
    for (const auto& arg : arguments)
      myargv.push_back((char*)arg.data());
    myargv.push_back(nullptr);
    return FitTest<LoopHelix>(myargv.size()-1,myargv.data(),sigmas);
  } else
  return FitTest<LoopHelix>(argc,argv,sigmas);
}
