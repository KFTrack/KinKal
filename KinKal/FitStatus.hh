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
    int iter_; // iteration number;
    status status_; // current status
    float chisq_; // current chisquared
    int ndof_; // current number of degrees of freedom
    float prob_; // chisquared probability
    bool usable() const { return status_ == converged || status_ == unconverged; }
    FitStatus() : iter_(-1), status_(needsfit), chisq_(std::numeric_limits<float>::max()), ndof_(0), prob_(-1.0){}
    static std::string statusName(status stat);
  };
  std::ostream& operator <<(std::ostream& os, FitStatus fitstatus );
}
#endif
