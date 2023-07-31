//
//  Description of a right conical frustrum for use in KinKal intersection
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Frustrum_hh
#define KinKal_Frustrum_hh
#include "KinKal/Geometry/Surface.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Rectangle.hh"
namespace KinKal {
  class Frustrum : public Surface {
    public:
      virtual ~Frustrum() {}
      // construct from necessary parameters
      Frustrum(VEC3 const& axis, VEC3 const& center, double r0, double r1, double halflen );
      // Surface interface
      bool onSurface(VEC3 const& point, double tol) const override;
      bool isInside(VEC3 const& point) const override;
      double curvature(VEC3 const& point) const override;
      bool inBounds(VEC3 const& point, double tol) const override;
      double distance(VEC3 const& point) const override;
      IntersectFlag intersect(Ray const& ray,double& dist, bool forwards, double tol) const override;
      VEC3 normal(VEC3 const& point) const override; // radially outward
      Plane tangentPlane(VEC3 const& point) const override;
      // frustrum-specific interface
      auto const& axis() const { return axis_; }
      auto const& center() const { return center_; }
      auto lowRadius() const { return r0_; } // radius at the lower end
      auto highRadius() const { return r1_; } // radius at the upper end
      auto centerRadius() const { return rmid_; } // radius at the center
      double radius(VEC3 point) const;
      double radius(double alen) const { return rmid_ + alen*drda_; } // radius at a given distance along the axis, WRT the center
      auto halfLength() const { return halflen_; }
      double halfAngle() const { return atan2(sint_,cost_); }
      // front and rear boundaries of this cylinder as disks.  U direction is arbitrary
      Disk frontDisk() const;
      Disk backDisk() const;
    private:
      VEC3 axis_; // symmetry axis of the cylinder
      VEC3 center_; // geometric center
      double r0_, r1_; // lowest and highest radius
      double halflen_; // half length
      // caches used to speed calculation
      double rmid_,drda_; // mid radius and rate of radial change
      double sint_, cost_; // sine and cosine of the cone 1/2 angle
      VEC3 uDirection() const; // arbitrary direction orthogonal to the axis
  };
}
std::ostream& operator <<(std::ostream& ost, KinKal::Frustrum const& cyl);
#endif
