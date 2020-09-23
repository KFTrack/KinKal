/*
   Original Author: S Middleton 2020
 */
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Tests/FitTest.hh"
int main(int argc, char *argv[]){
  vector<double> sigmas = { 3.0, 0.002, 3.5, 0.001, 0.5, 0.05, 0.05}; // expected parameter sigmas: the last is momentum 
  if(argc == 1){
    cout << "Adding momentum constraint" << endl;
    std::vector<std::string> arguments;
    arguments.push_back(argv[0]);
    arguments.push_back("--constrainpar");
    arguments.push_back("5");
    std::vector<char*> myargv;
    for (const auto& arg : arguments)
      myargv.push_back((char*)arg.data());
    myargv.push_back(nullptr);
    return FitTest<KinematicLine>(myargv.size()-1,myargv.data(),sigmas);
  } else 
    return FitTest<KinematicLine>(argc,argv,sigmas);
}

