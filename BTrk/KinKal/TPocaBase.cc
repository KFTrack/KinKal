#include "BTrk/KinKal/TPOCABase.hh"
namespace KinKal {
  std::vector<std::string> TPOCABase::statusNames_ = { "converged", "unconverged", "outsiderange", "pocafailed", "derivfailed", "unknown"};
  std::string const& TPOCABase::statusName(TPStat status) { return statusNames_[static_cast<size_t>(status)];}
}
