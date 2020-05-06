#ifndef KinKal_FitStatus_hh
#define KinKal_FitStatus_hh
#include <ostream>
#include <string>
#include <vector>
#include <limits>

namespace KinKal {
// struct to define fit status
  struct FitStatus {
    enum status {needsfit=-1,converged,unconverged,oscillating,diverged,lowNDOF,failed}; // fit status
    int miter_; // meta-iteration number;
    int iter_; // iteration number;
    status status_; // current status
    double chisq_; // current chisquared
    unsigned ndof_; // current number of degrees of freedom
    double prob_; // chisquared probability
    std::string comment_; // further information about the status 
    bool usable() const { return status_ !=failed; }
    bool needsFit() const { return status_ == needsfit || status_ == unconverged; }
    FitStatus(unsigned miter) : miter_(miter), iter_(-1), status_(needsfit), chisq_(std::numeric_limits<double>::max()), ndof_(0), prob_(-1.0){}
    static std::string statusName(status stat);
  };
  std::ostream& operator <<(std::ostream& os, FitStatus fitstatus );
}
#endif
