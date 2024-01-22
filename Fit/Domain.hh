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
      Domain(double lowtime, double range, VEC3 const& bnom, double tol) : range_(lowtime,lowtime+range), bnom_(bnom), tol_(tol) {}
      Domain(TimeRange const& range, VEC3 const& bnom, double tol) : range_(range), bnom_(bnom), tol_(tol) {}
      bool operator < (Domain const& other) const {return begin() < other.begin(); }
      auto const& range() const { return range_; }
      // forward range functions
      double begin() const { return range_.begin();}
      double end() const { return range_.end();}
      double mid() const { return range_.mid();}
      double tolerance() const { return tol_; }
      auto const& bnom() const { return bnom_; }
      void updateBNom( VEC3 const& bnom) { bnom_ = bnom; }; // used in DomainWall updating
    private:
      TimeRange range_; // range of this domain
      VEC3 bnom_; // nominal BField for this domain
      double tol_; // tolerance used to create this domain
  };
}
#endif

