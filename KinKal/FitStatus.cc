#include "KinKal/FitStatus.hh"
namespace KinKal {

  std::string FitStatus::statusName(FitStatus::status stat) {
    switch(stat) {
      case FitStatus::needsfit: default :
	return "NeedsFit";
      case FitStatus::unconverged: 
	return "Unconverged";
      case FitStatus::converged: 
	return "Converged";
      case FitStatus::oscillating: 
	return "Oscillating";
      case FitStatus::diverged: 
	return "Diverged";
      case FitStatus::lowNDOF: 
	return "LowNDOF";
      case FitStatus::failed: 
	return "Failed";
    }
  }

  std::ostream& operator <<(std::ostream& ost, FitStatus fitstatus ) {
    ost << "Fit Status " << FitStatus::statusName(fitstatus.status_)
      << fitstatus.comment_
      << " Meta-iteration " << fitstatus.miter_
      << " iteration " << fitstatus.iter_
      << " chisq " << fitstatus.chisq_ 
      << " NDOF " << fitstatus.ndof_
      << " prob " << fitstatus.prob_;
    return ost;
  }
}
