/*
   Original Author: S Middleton 2020
 */
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Tests/FitTest.hh"
int main(int argc, char *argv[]){
  KinKal::DVEC sigmas(0.5, 0.004, 0.5, 0.002, 0.4, 0.05); // expected parameter sigmas
  if(argc == 1){
    cout << "Adding momentum constraint" << endl;
    std::vector<std::string> arguments;
    arguments.push_back(argv[0]);
    arguments.push_back("--constrainpar");
    arguments.push_back("5");
    arguments.push_back("--Bz");
    arguments.push_back("0.0");
    arguments.push_back("--Schedule");
    arguments.push_back("driftfit.txt");
    std::vector<char*> myargv;
    for (const auto& arg : arguments)
      myargv.push_back((char*)arg.data());
    myargv.push_back(nullptr);
    return FitTest<KinematicLine>(myargv.size()-1,myargv.data(),sigmas);
  } else
    return FitTest<KinematicLine>(argc,argv,sigmas);
}

