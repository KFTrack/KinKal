#ifndef KinKal_Domain_hh
#define KinKal_Domain_hh
//
//  domain used in computing magnetic field corrections: just a time range and the tolerance used to define it
//
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Vectors.hh"
namespace KinKal {
  class Domain : public TimeRange {
    public:
      Domain(double lowtime, double range, VEC3 const& bnom, double tol) : TimeRange(lowtime,lowtime+range), bnom_(bnom), tol_(tol) {}
      Domain(TimeRange const& range, VEC3 const& bnom, double tol) :TimeRange(range), bnom_(bnom), tol_(tol) {}
      bool operator < (Domain const& other) const {return begin() < other.begin(); }
      double tolerance() const { return tol_; }
      void updateBNom( VEC3 const& bnom) { bnom_ = bnom; }; // used in DomainWall updating
    private:
      VEC3 bnom_; // nominal BField for this domain
      double tol_; // tolerance used to create this domain
  };
}
#endif

