#include "Trajectory/LoopHelix.hh"
#include "UnitTests/FitTest.hh"
int main(int argc, char **argv) {
  vector<double> sigmas = { 0.5, 0.5, 0.5, 0.5, 0.002, 0.5, 0.1}; // expected parameter sigmas: the last is momentum 
  if(argc == 1){
    cout << "Testing gradient field, correction 2" << endl;
    std::vector<std::string> arguments;
    arguments.push_back(argv[0]);
    arguments.push_back("--Bgrad");
    arguments.push_back("-0.036"); // mu2e-like field gradient
    arguments.push_back("--bfcor");
    arguments.push_back("2"); // local field correction (BField rotation)
    arguments.push_back("--tolerance");
    arguments.push_back("0.01"); // currently required as tolerance doesn't take into account rotation lever arm FIXME!
    std::vector<char*> myargv;
    for (const auto& arg : arguments)
      myargv.push_back((char*)arg.data());
    myargv.push_back(nullptr);
    return FitTest<LoopHelix>(myargv.size()-1,myargv.data(),sigmas);
  } else 
  return FitTest<LoopHelix>(argc,argv,sigmas);
}
