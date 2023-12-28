#ifndef KinKal_Domain_hh
#define KinKal_Domain_hh
//
//  Struct to define a domain used in computing magnetic field corrections: just the time range and the tolerance used to define it
//
#include "KinKal/General/TimeRange.hh"
namespace KinKal {
  struct Domain : public TimeRange {
    double tol_; // tolerance used to create this domain
    Domain(TimeRange const& trange, double tol) : TimeRange(trange), tol_(tol) {}
  }
}
#endif

