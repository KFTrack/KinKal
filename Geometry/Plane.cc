#include "KinKal/Geometry/Plane.hh"
namespace KinKal {
  bool Plane::onSurface(VEC3 const& point, double tol) const {
    return fabs(norm_.Dot(point-center_)) < tol;
  }

  bool Plane::isInside(VEC3 const& point) const {
    return norm_.Dot(point-center_) > 0.0;
  }

  IntersectFlag Plane::intersect(Ray const& ray,double& dist, bool forwards, double tol) const {
    IntersectFlag retval;
    double ddir = norm_.Dot(ray.dir_);
    if(fabs(ddir)>0.0) {
      double pdist = norm_.Dot(center_ - ray.start_);
      dist = pdist/ddir;
      if(dist > 0.0 || !forwards){
        retval.onsurface_ = true;
        retval.inbounds_ = inBounds(ray.position(dist),tol);
      }
    }
    return retval;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Plane const& plane) {
  ost << "Plane with center " << plane.center() << " , normal " << plane.normal();
  return ost;
}
