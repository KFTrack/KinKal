//
//  Description of a rectangular planar section
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Rectangle_hh
#define KinKal_Rectangle_hh
#include "KinKal/Geometry/Plane.hh"
namespace KinKal {
  class Rectangle : public Plane {
    public:
      virtual ~Rectangle() {};
      // construct from necessary parameters
      Rectangle(VEC3 const& norm, VEC3 const& udir, VEC3 const& center, double uhalflen, double vhalflen);
      // surface interface
      bool inBounds(VEC3 const& point, double tol) const override;
      // rectangle-specific interface
      double uHalfLength() const { return uhalflen_; }
      double vHalfLength() const { return vhalflen_; }
    private:
      double uhalflen_, vhalflen_; // u and v half-lengths
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Rectangle const& rect);
#endif
