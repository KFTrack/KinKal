#include "KinKal/Geometry/Cylinder.hh"
#include "Math/VectorUtil.h"
#include <exception>
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

  VEC3 Cylinder::uDirection() const {
    // u direction is arbitrary; use 'x' if the axis is mostly along z
    static const VEC3 xdir(1.0,0.0,0.0);
    static const VEC3 ydir(0.0,1.0,0.0);
    static const VEC3 zdir(0.0,0.0,1.0);
    double zdot  = fabs(axis_.Dot(zdir));
    if(zdot > 0.5){
      return ydir.Cross(axis_).Unit();
    } else if(zdot > 0.01) {
      return axis_.Cross(zdir).Unit();
    } else {
      double xdot = fabs(axis_.Dot(xdir));
      if(xdot > 0.5)
        return ydir.Cross(axis_).Unit();
      else
        return xdir.Cross(axis_).Unit();
    }
  }

  Disk Cylinder::frontDisk() const {
      return Disk(axis_,uDirection(),center_-halflen_*axis_,radius_);
  }

  Disk Cylinder::backDisk() const {
    return Disk(axis_,uDirection(),center_+halflen_*axis_,radius_);
  }

  Rectangle Cylinder::inscribedRectangle(VEC3 const& norm) const {
    // make sure normal is perpendicular to the axis
    if(axis_.Dot(norm) > 1e-10) throw std::invalid_argument("normal not perpendicular to axis");
    // U points along the cylinder
    return Rectangle(norm,axis_,center_,halflen_,radius_);
  }

  Rectangle Cylinder::tangentRectangle(VEC3 const& spoint) const {
    // rectangle normal is the local cylinder normal
    auto norm = normal(spoint);
    // rectangle center is on the cylinder
    VEC3 rcenter = center_ + norm*radius_;
    return Rectangle(norm,axis_,rcenter,halflen_,radius_);
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
