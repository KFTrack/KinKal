#include "Trajectory/ClosestApproachData.hh"
namespace KinKal {
  const std::vector<std::string> ClosestApproachData::statusNames_ = { "converged", "unconverged", "oscillating", "diverged","pocafailed", "invalid"};
  std::string const& ClosestApproachData::statusName(TPStat status) { return statusNames_[static_cast<size_t>(status)];}
}
