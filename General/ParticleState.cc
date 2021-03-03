#include "KinKal/General/ParticleState.hh"
namespace KinKal {
  using std::string;
  using std::vector;
  const vector<string> ParticleState::stateTitles_ = {
    "X Position",
    "Y Position",
    "Z Position",
    "X Momentum",
    "Y Momentum",
    "Z Momentum"};
  const vector<string> ParticleState::stateNames_ = { "XPos", "YPos", "ZPos", "XMom", "YMom", "ZMom"};
  const vector<string> ParticleState::stateUnits_ = { "mm", "mm", "mm", "MeV/c", "MeV/c", "MeV/c"};
  string const& ParticleState::stateName(size_t index) { return stateNames_[index];}
  string const& ParticleState::stateUnit(size_t index) { return stateUnits_[index];}
  string const& ParticleState::stateTitle(size_t index) { return stateTitles_[index];}

  double ParticleStateEstimate::momentumVar() const {
    auto momdir = state_.momentum3().Unit();
    DVEC dMdm(0.0, 0.0, 0.0, momdir.X(), momdir.Y(), momdir.Z());
    return ROOT::Math::Similarity(dMdm,scovar_);
  }
}
