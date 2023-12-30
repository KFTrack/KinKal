#ifndef KinKal_Domain_hh
#define KinKal_Domain_hh
//
//  domain used in computing magnetic field corrections: just a time range and the tolerance used to define it
//
#include "KinKal/General/TimeRange.hh"
namespace KinKal {
  class Domain : public TimeRange {
    public:
      Domain(double lowtime, double range, double tol) : TimeRange(lowtime,lowtime+range), tol_(tol) {}
      Domain(TimeRange const& range, double tol) :TimeRange(range),  tol_(tol) {}
      bool operator < (Domain const& other) const {return begin() < other.begin(); }
      double tolerance() const { return tol_; }
    private:
      double tol_; // tolerance used to create this domain
  };
}
#endif

