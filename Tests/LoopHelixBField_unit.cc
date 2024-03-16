#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Tests/BFieldMapTest.hh"
int main(int argc, char **argv) {
  if(argc == 1){
    cout << "Testing gradient field, correction 2" << endl;
    std::vector<std::string> arguments;
    arguments.push_back(argv[0]);
    arguments.push_back("--Bgrad");
    arguments.push_back("-0.036"); // mu2e-like field gradient
    arguments.push_back("--Tol");
    arguments.push_back("1.0e-4");
    arguments.push_back("1");
    std::vector<char*> myargv;
    for (const auto& arg : arguments)
      myargv.push_back((char*)arg.data());
    myargv.push_back(nullptr);
    return BFieldMapTest<LoopHelix>(myargv.size()-1,myargv.data());
  } else
    return BFieldMapTest<LoopHelix>(argc,argv);
}
