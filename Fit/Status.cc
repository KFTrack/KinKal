#include "KinKal/Fit/Status.hh"
namespace KinKal {

  std::string Status::statusName(Status::status stat) {
    switch(stat) {
      case Status::unfit: default :
	return "Unfit ";
      case Status::unconverged: 
	return "Unconverged ";
      case Status::converged: 
	return "Converged ";
      case Status::diverged: 
	return "Diverged ";
      case Status::lowNDOF: 
	return "LowNDOF ";
      case Status::failed: 
	return "Failed ";
    }
  }

  std::ostream& operator <<(std::ostream& ost, Status const& fitstatus ) {
    ost << "Fit Status " << Status::statusName(fitstatus.status_)
      << fitstatus.comment_
      << " Meta-iteration " << fitstatus.miter_
      << " iteration " << fitstatus.iter_
      <<  " " << fitstatus.chisq_;
    return ost;
  }
}
