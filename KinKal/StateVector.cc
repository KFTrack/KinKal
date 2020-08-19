#include "KinKal/StateVector.hh"
namespace KinKal {
  using std::string;
  using std::vector;
  vector<string> StateVector::stateTitles_ = {
    "X Position",
    "Y Position",
    "Z Position",
    "X Momentum",
    "Y Momentum",
    "Z Momentum"};
  vector<string> StateVector::stateNames_ = { "XPos", "YPos", "ZPos", "XMom", "YMom", "ZMom"};
  vector<string> StateVector::stateUnits_ = { "mm", "mm", "mm", "MeV/c", "MeV/c", "MeV/c"};
  string const& StateVector::stateName(size_t index) { return stateNames_[index];}
  string const& StateVector::stateUnit(size_t index) { return stateUnits_[index];}
  string const& StateVector::stateTitle(size_t index) { return stateTitles_[index];}
}
