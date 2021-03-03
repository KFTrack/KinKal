#include "KinKal/Trajectory/ClosestApproachData.hh"
namespace KinKal {
  const std::vector<std::string> ClosestApproachData::statusNames_ = { "converged", "unconverged", "oscillating", "diverged","pocafailed", "invalid"};
  std::string const& ClosestApproachData::statusName(TPStat status) { return statusNames_[static_cast<size_t>(status)];}

  std::ostream& operator << (std::ostream& ost, ClosestApproachData const& cadata) {
    ost << "DOCA = " << cadata.doca() << " +- " << sqrt(cadata.docaVar()) << " sign = " << cadata.lSign()
    << " DeltaT = " << cadata.deltaT() << " +- " << sqrt(cadata.tocaVar())
    << " pdir.Dot(sdir) = " << cadata.dirDot();
    return ost;
  }
}
