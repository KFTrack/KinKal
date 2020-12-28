#ifndef KinKal_Status_hh
#define KinKal_Status_hh
#include "KinKal/Fit/Chisq.hh"
#include <ostream>
#include <string>
#include <vector>

namespace KinKal {
// struct to define fit status
  struct Status {
    enum status {unfit=-1,converged,unconverged,oscillating,diverged,lowNDOF,failed}; // fit status
    int miter_; // meta-iteration number;
    int iter_; // iteration number;
    status status_; // current status
    Chisq chisq_; // current chisquared
    std::string comment_; // further information about the status 
    bool usable() const { return status_ !=failed &&  status_ !=diverged; }
    bool needsFit() const { return status_ == unfit || status_ == unconverged; }
    Status(unsigned miter) : miter_(miter), iter_(-1), status_(unfit){}
    static std::string statusName(status stat);
  };
  std::ostream& operator <<(std::ostream& os, Status const& fitstatus );
}
#endif
