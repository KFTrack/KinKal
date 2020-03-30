#ifndef KinKal_FitStatus_hh
#define KinKal_FitStatus_hh
#include <ostream>
#include <string>
#include <vector>
namespace KinKal {
// struct to define fit status
  struct FitStatus {
    enum status {needsfit=0,converged,unconverged,oscillating,diverged,failed}; // fit status
    unsigned niter_; // number of iterations executed;
    status status_; // current status
    double chisq_; // current chisquared
    int ndof_; // number of degrees of freedom
    FitStatus() : niter_(0), status_(needsfit), chisq_(0.0), ndof_(0) {}
    bool iterate() const { return status_ == needsfit; }
    static std::string statusName(status stat) { return statnames[stat]; }
    static std::vector<std::string> statnames;
  };
  std::ostream& operator <<(std::ostream& os, FitStatus fitstatus );
}
#endif
