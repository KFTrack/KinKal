#ifndef KinKal_Domain_hh
#define KinKal_Domain_hh
//
//  domain used to compute magnetic field corrections. Magnetic bending not described by the intrinsic parameterization is assumed
//  to be negligible over the domain
//
#include <memory>
#include "KinKal/General/CloneContext.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Vectors.hh"
namespace KinKal {
  class Domain : public TimeRange {
    public:
      Domain(double lowtime, double range, VEC3 const& bnom) : range_(lowtime,lowtime+range), bnom_(bnom) {}
      Domain(TimeRange const& range, VEC3 const& bnom) : range_(range), bnom_(bnom) {}
      // clone op for reinstantiation
      Domain(Domain const&);
      std::shared_ptr< Domain > clone(CloneContext&) const;
      bool operator < (Domain const& other) const {return begin() < other.begin(); }
      auto const& range() const { return range_; }
      // forward range functions
      double begin() const { return range_.begin();}
      double end() const { return range_.end();}
      double mid() const { return range_.mid();}
      auto const& bnom() const { return bnom_; }
      void updateBNom( VEC3 const& bnom) { bnom_ = bnom; }; // used in DomainWall updating
    private:
      TimeRange range_; // range of this domain
      VEC3 bnom_; // nominal BField for this domain
  };

  // clone op for reinstantiation
  Domain::Domain(Domain const& rhs):
      range_(rhs.range()),
      bnom_(rhs.bnom()){
  }

  std::shared_ptr<Domain> Domain::clone(CloneContext& context) const{
    auto rv = std::make_shared<Domain>(*this);
    return rv;
  }
}
#endif

