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
      Plane(VEC3 const& norm, VEC3 const& udir, VEC3 const& center);
      // surface interface
      bool onSurface(VEC3 const& point, double tol) const override;
      // inside is defined as on the opposite side as the normal
      // the absolute value should never matter
      bool isInside(VEC3 const& point) const override;
      bool inBounds(VEC3 const& point, double tol) const override { return true; }
      double distance(VEC3 const& point) const override;
      double curvature(VEC3 const& point) const override { return 0.0; }
      IntersectFlag intersect(Ray const& ray,double& dist, bool forwards, double tol) const override;
      VEC3 normal(VEC3 const& point) const override { return norm_; }
      auto const& uDirection() const { return udir_; }
      auto const& vDirection() const { return vdir_; }
      auto const& normal() const { return norm_; }
      Plane tangentPlane(VEC3 const& ) const override { return *this; }
      // plane interfac3
      auto const& center() const { return center_; }
    private:
      // note that UVW forms a right-handed orthonormal coordinate system
      VEC3 norm_; // normal to the plane (W direction)
      VEC3 udir_; // U direction: perpendicular to the normal
      VEC3 vdir_; // V direction: perpendicular to the normal and U
      VEC3 center_; // point on the plane
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Plane const& plane);
#endif
