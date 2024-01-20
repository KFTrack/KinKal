//
//  Description of a right circular cylindrical, for use in KinKal intersection
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Cylinder_hh
#define KinKal_Cylinder_hh
#include "KinKal/Geometry/Surface.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Rectangle.hh"
namespace KinKal {
  class Cylinder : public Surface {
    public:
      virtual ~Cylinder() {}
      // construct from necessary parameters
      Cylinder(VEC3 const& axis, VEC3 const& center, double radius, double halflen );
      // Surface interface
      bool onSurface(VEC3 const& point, double tol) const override;
      bool isInside(VEC3 const& point) const override;
      double curvature(VEC3 const& point) const override { return 1.0/radius_; }
      bool inBounds(VEC3 const& point, double tol) const override;
      double distance(VEC3 const& point) const override;
      IntersectFlag intersect(Ray const& ray,double& dist, bool forwards, double tol) const override;
      VEC3 normal(VEC3 const& point) const override; // radially outward
      Plane tangentPlane(VEC3 const& point) const override;
      // cylinder-specific interface
      auto const& axis() const { return axis_; }
      auto const& center() const { return center_; }
      double radius() const { return radius_; }
      double halfLength() const { return halflen_; }
      // front and rear boundaries of this cylinder as disks.  U direction is arbitrary
      Disk frontDisk() const;
      Disk midDisk() const;
      Disk backDisk() const;
      // inscribed rectangle with given normal direction.  U direction is long the axis
      Rectangle inscribedRectangle(VEC3 const& norm) const;
    private:
      VEC3 axis_; // symmetry axis of the cylinder
      VEC3 center_; // geometric center
      double radius_; // transverse radius
      double halflen_; // half length
      double radius2_; // squared radius (cache);
      VEC3 uDirection() const; // arbitrary direction orthogonal to the axis
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Cylinder const& cyl);
#endif
