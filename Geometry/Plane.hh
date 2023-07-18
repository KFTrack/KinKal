//
//  Description of an unbounded plane in space.
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Plane_hh
#define KinKal_Plane_hh
#include "KinKal/Geometry/Surface.hh"
namespace KinKal {
  class Plane : public Surface {
    public:
      virtual ~Plane() {};
      // construct from necessary parameters
      Plane(VEC3 const& norm, VEC3 const& center) : norm_(norm.Unit()), center_(center){}
      // surface interface
      bool onSurface(VEC3 const& point, double tol) const override;
      bool isInside(VEC3 const& point) const override; // defined as on the same side as the normal
      bool inBounds(VEC3 const& point, double tol) const override { return true; }
      double curvature(VEC3 const& point) const override { return 0.0; }
      IntersectFlag intersect(Ray const& ray,double& dist, bool forwards, double tol) const override;
      VEC3 normal(VEC3 const& point) const override { return norm_; }
      auto const& normal() const { return norm_; }
      auto const& center() const { return center_; }
    private:
      VEC3 norm_; // normal to the plane
      VEC3 center_; // point on the plane
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Plane const& plane);
#endif
