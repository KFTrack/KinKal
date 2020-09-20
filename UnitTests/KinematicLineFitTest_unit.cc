/*
   Original Author: S Middleton 2020
 */
#include "Trajectory/KinematicLine.hh"
#include "UnitTests/FitTest.hh"
int main(int argc, char *argv[]){
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
    return FitTest<KinematicLine>(myargv.size()-1,myargv.data());
  } else 
    return FitTest<KinematicLine>(argc,argv);
}

