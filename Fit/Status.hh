#ifndef KinKal_Status_hh
#define KinKal_Status_hh
#include "KinKal/General/Chisq.hh"
#include <ostream>
#include <string>
#include <vector>

namespace KinKal {
// struct to define fit status
  struct Status {
    enum status {unfit=-1,converged,unconverged,lowNDOF,gapdiverged,paramsdiverged,chisqdiverged,failed}; // fit status
    unsigned miter_; // meta-iteration number;
    unsigned iter_; // iteration number;
    status status_; // current status
    Chisq chisq_; // current chisquared
    std::string comment_; // further information about the status
    bool usable() const { return status_ > unfit && status_ < lowNDOF; }
    bool needsFit() const { return status_ == unfit || status_ == unconverged; }
    Status(unsigned miter,unsigned iter=0) : miter_(miter), iter_(iter), status_(unfit){}
    static std::string statusName(status stat);
  };
  std::ostream& operator <<(std::ostream& os, Status const& fitstatus );
}
#endif
