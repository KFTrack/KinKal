#include "KinKal/FitStatus.hh"
namespace KinKal {

  std::vector<std::string> FitStatus::statnames{ "needsfit", "converged", "unconverged", "failed"};

  std::ostream& operator <<(std::ostream& os, FitStatus fitstatus ) {
    os << "Fit Status " << FitStatus::statusName(fitstatus.status_)
      << " iterations " << fitstatus.niter_
      << " chisq " << fitstatus.chisq_ 
      << " NDOF " << fitstatus.ndof_;
    return os;
  }
}
