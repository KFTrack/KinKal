#include "KinKal/Geometry/Cylinder.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace KinKal {
  bool Cylinder::onSurface(VEC3 const& point, double tol) const {
    auto rvec = point - center_;
    auto pvec = PerpVector(rvec,axis_);
    return fabs(pvec.R()-radius_) < tol;
  }

  bool Cylinder::inBounds(VEC3 const& point, double tol) const {
    auto rvec = point - center_;
    return fabs(rvec.Dot(axis_)) < halflen_ + tol;
  }

  bool Cylinder::isInside(VEC3 const& point) const {
    auto rvec = point - center_;
    return Perp2(rvec,axis_) < radius2_;
  }

  VEC3 Cylinder::normal(VEC3 const& point) const {
    // normal is perpendicular part of the difference
    auto rvec = point - center_;
    auto pvec = PerpVector(rvec,axis_);
    return pvec.Unit();
  }

  IntersectFlag Cylinder::intersect(Ray const& ray,double& dist, bool forwards, double tol) const {
    IntersectFlag retval;
    double ddot = ray.dir_.Dot(axis_);
    double alpha = (1.0 - ddot*ddot); // always positive
                                      // make sure the ray isn't co-linear on the relevant scale
    if(alpha > tol/std::max(radius_,halflen_)){
      auto rvec = ray.start_ - center_;
      double sdot = rvec.Dot(axis_);
      double beta = sdot*ddot - rvec.Dot(ray.dir_);
      double gamma = rvec.Mag2() - sdot*sdot - radius2_;
      double beta2 = beta*beta;
      double ag = alpha*gamma;
      // make sure there's a solution
      if(beta2 > ag){
        // choose the solution based on the strategy
        double delta = sqrt(beta2 - ag);
        double d1 = (beta - delta)/alpha;
        double d2 = (beta + delta)/alpha;
        if(!forwards){
          // closest solution
          dist = fabs(d1) < fabs(d2) ? d1 : d2;
        } else {
          // closest forwards solution
          if(d1 > 0.0){
            retval.onsurface_ = true;
            dist = d1;
          } else if(d2 > 0.0) {
            retval.onsurface_ = true;
            dist = d2;
          }
        }
      }
      // test that the point is on the cylinder
      if(retval.onsurface_){
        auto point = ray.position(dist);
        retval.inbounds_ = inBounds(point,tol);
      }
    }
    return retval;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Cylinder const& cyl) {
  ost << "Cylinder with center = " << cyl.center() << " , axis " << cyl.axis() << " radius " << cyl.radius() << " half-length " << cyl.halfLength();
  return ost;
}
