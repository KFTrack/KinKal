//
//  Description of a rectangular planar section
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Rectangle_hh
#define KinKal_Rectangle_hh
#include "KinKal/Geometry/Plane.hh"
#include <exception>
namespace KinKal {
  class Rectangle : public Plane {
    public:
      virtual ~Rectangle() {};
      // construct from necessary parameters
      Rectangle(VEC3 const& norm, VEC3 const& center, VEC3 const& uaxis, double uhalflen, double vhalflen);
      // surface interface
      bool inBounds(VEC3 const& point, double tol) const override;
      // rectangle-specific interface
      auto const& uDirection() const { return udir_; }
      auto const& vDirection() const { return vdir_; }
      double uHalfLength() const { return uhalflen_; }
      double vHalfLength() const { return vhalflen_; }
    private:
      double uhalflen_, vhalflen_; // u and v half-lengths
      VEC3 udir_; // U direction: perpendicular to the normal
      VEC3 vdir_; // UV direction: perpendicular to the normal:w
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Rectangle const& rect);
#endif
