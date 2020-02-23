#include "KinKal/TPocaBase.hh"
namespace KinKal {
  std::vector<std::string> TPocaBase::statusNames_ = { "converged", "unconverged", "outsiderange", "pocafailed", "derivfailed", "unknown"};
  std::string const& TPocaBase::statusName(TPStat status) { return statusNames_[static_cast<size_t>(status)];}
}
