#include "KinKal/TPocaData.hh"
namespace KinKal {
  std::vector<std::string> TPocaData::statusNames_ = { "converged", "unconverged", "oscillating", "diverged","pocafailed", "invalid"};
  std::string const& TPocaData::statusName(TPStat status) { return statusNames_[static_cast<size_t>(status)];}
}
