#include "KinKal/Geometry/Plane.hh"
#include <exception>
namespace KinKal {
  Plane::Plane(VEC3 const& norm, VEC3 const& udir, VEC3 const& center) : norm_(norm.Unit()), udir_(udir.Unit()), center_(center){
    // check that U is perpendicular
    if(udir_.Dot(normal()) > 1e-10) throw std::invalid_argument("U direction not perpendicular to normal");
    // V direction is implicit, to create a right-handed coordinate system
    vdir_ = normal().Cross(udir_);
  }

  bool Plane::onSurface(VEC3 const& point, double tol) const {
    return fabs(norm_.Dot(point-center_)) < tol;
  }

  bool Plane::isInside(VEC3 const& point) const {
    return norm_.Dot(point-center_) < 0.0;
  }

  double Plane::distance(VEC3 const& point) const {
    return norm_.Dot(point-center_);
  }

  IntersectFlag Plane::intersect(Ray const& ray,double& dist, bool forwards, double tol) const {
    IntersectFlag retval;
    double ddir = norm_.Dot(ray.direction());
    if(fabs(ddir)>0.0) {
      double pdist = norm_.Dot(center_ - ray.start());
      dist = pdist/ddir;
      if(dist > 0.0 || !forwards){
        retval.onsurface_ = true;
        retval.inbounds_  = inBounds(ray.position(dist),tol);
      }
    }
    return retval;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Plane const& plane) {
  ost << "Plane with center " << plane.center() << " normal " << plane.normal()
    << " U direction " << plane.uDirection() << " V direction " << plane.vDirection();
  return ost;
}
