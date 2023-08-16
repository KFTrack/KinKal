//
//  Abstract interface to a 2-dimensional surface in 3-space
//  These are used in reconstruction to describe materials and
//  common reference locations
//  original author: David Brown (LBN) 2023
//
#ifndef KinKal_Surface_hh
#define KinKal_Surface_hh
#include "KinKal/Geometry/Ray.hh"
#include "KinKal/Geometry/IntersectFlag.hh"
#include "KinKal/General/Vectors.hh"
namespace KinKal {
  class Surface {
    public:
      virtual ~Surface() {}
      // determine if the given point is on the surface, within tolerance
      virtual bool onSurface(VEC3 const& point, double tol) const = 0;
      // determine if a point is 'inside' the surface.  For open surfaces this tests one side or the other
      virtual bool isInside(VEC3 const& point) const = 0;
      // describe the local surface curvature, defined as the maximum 2nd derivate along the surface nearest that point
      virtual double curvature(VEC3 const& point) const = 0;
      // determine if a point on the surface is in bounds
      virtual bool inBounds(VEC3 const& point, double tol) const = 0;
      // find the distance along a ray where it would intersect this surface; Returned flag describes what happened
      virtual IntersectFlag intersect(Ray const& ray,double& dist, bool forwards, double tol) const = 0;
      // find the normal to the surface at the given point.  Direction convention is surface-dependent
      virtual VEC3 normal(VEC3 const& point) const  = 0;
  };
}
#endif
