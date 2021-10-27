#include "KinKal/MatEnv/SimpleFileFinder.hh"
#include <cstdlib>
#include <stdexcept>
#include <iostream>
namespace MatEnv {
  SimpleFileFinder::SimpleFileFinder() : SimpleFileFinder(std::string("KINKAL_SOURCE_DIR"),std::string("/MatEnv/"))
  {}

  SimpleFileFinder::SimpleFileFinder(std::string const& project, std::string const& dir) :
    elemlist_("ElementsList.data"),
    isolist_("IsotopesList.data"),
    matlist_("MaterialsList.data")
  {
    const char* src = getenv(project.c_str());
    if(!src) throw std::invalid_argument("KINKAL_SOURCE_DIR not set: did you forget to source setup.sh?");
    std::string source = std::string(src);
    std::cout <<" project " << project << " source " << source << " dir " << dir << std::endl;
    project_ = source + dir;
  }
}
