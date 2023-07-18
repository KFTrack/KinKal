//
//  Description of an annular planar section
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Annulus_hh
#define KinKal_Annulus_hh
#include "KinKal/Geometry/Plane.hh"
namespace KinKal {
  class Annulus : public Plane {
    public:
      ~Annulus() {};
      // construct from necessary parameters
      Annulus(VEC3 const& norm, VEC3 const& center, double innerrad, double outerrad) : Plane(norm,center), irad_(innerrad), orad_(outerrad) {}
      // surface interface
      bool inBounds(VEC3 const& point, double tol) const override;
      // annulus-specific interface
      auto innerRadius() const { return irad_; }
      auto outerRadius() const { return orad_; }
    private:
      double irad_, orad_; // inner and outer radii
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Annulus const& annulus);
#endif
