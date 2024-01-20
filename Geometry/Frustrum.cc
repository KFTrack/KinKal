#include "KinKal/Geometry/Frustrum.hh"
#include "Math/VectorUtil.h"
#include <exception>
using namespace ROOT::Math::VectorUtil;
namespace KinKal {

  Frustrum::Frustrum(VEC3 const& axis, VEC3 const& center, double r0, double r1, double halflen ) : axis_(axis.Unit()), center_(center), r0_(r0), r1_(r1), halflen_(halflen), rmid_(0.5*(r0+r1)), drda_( 0.5*(r1-r0)/halflen) {
    if(r0 <= 0.0 || r1 <= 0.0 || halflen <= 0.0 ) throw std::invalid_argument("Invalid Frustrum arguments");
    double dr = r1_ - r0_;
    // compute sine and cosine of the cone 1/2 angle
    double len = sqrt( dr*dr + 4*halflen_*halflen_);
    cost_ = 2*halflen_/len;
    sint_ = dr/len;
  }

  double Frustrum::radius(VEC3 point) const {
    auto rvec = point - center_;
    double alen = rvec.Dot(axis_);
    return radius(alen);
  }

  bool Frustrum::onSurface(VEC3 const& point, double tol) const {
    auto rvec = point - center_;
    auto pvec = PerpVector(rvec,axis_);
    double alen = rvec.Dot(axis_);
    return fabs(pvec.R()-radius(alen)) < tol;
  }

  bool Frustrum::inBounds(VEC3 const& point, double tol) const {
    auto rvec = point - center_;
    return fabs(rvec.Dot(axis_)) < halflen_ + tol/cost_;
  }

  bool Frustrum::isInside(VEC3 const& point) const {
    auto rvec = point - center_;
    double alen = rvec.Dot(axis_);
    return Perp(rvec,axis_) < radius(alen);
  }

  double Frustrum::curvature(VEC3 const& point) const {
    auto rvec = point - center_;
    double alen = rvec.Dot(axis_);
    return 1.0/radius(alen);
  }

  double Frustrum::distance(VEC3 const& point) const {
    auto rvec = point - center_;
    double alen = rvec.Dot(axis_);
    return (Perp(rvec,axis_) - radius(alen))*cost_;
  }

  VEC3 Frustrum::normal(VEC3 const& point) const {
    // normal is perpendicular part of the difference
    auto rvec = point - center_;
    auto pvec = PerpVector(rvec,axis_);
    return pvec.Unit()*cost_ - axis_*sint_;
  }

  VEC3 Frustrum::uDirection() const {
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

  Disk Frustrum::frontDisk() const {
    return Disk(axis_,uDirection(),center_-halflen_*axis_,radius(-halflen_));
  }

  Disk Frustrum::midDisk() const {
    return Disk(axis_,uDirection(),center_,radius(0.0));
  }

  Disk Frustrum::backDisk() const {
    return Disk(axis_,uDirection(),center_+halflen_*axis_,radius(halflen_));
  }

  Plane Frustrum::tangentPlane(VEC3 const& spoint) const {
    auto rvec = spoint - center_;
    double alen = rvec.Dot(axis_);
    // rectangle normal is the local cone surface normal
    auto norm = normal(spoint);
    // correct for any tolerance if the point isnt exactly on the surface (up to FP accuracy)
    double rad = norm.Dot(spoint - center_);
    auto rcent = spoint + norm*(radius(alen)-rad);
    // v direction is along the cone
    auto udir = axis_.Cross(norm);
    return Plane(norm,udir,rcent);
  }

  IntersectFlag Frustrum::intersect(Ray const& ray,double& dist, bool forwards, double tol) const {
    // exact solution
    IntersectFlag retval;
    // displace the ray start to be WRT the base of the cone
    auto spos = ray.start() - center_ + halflen_*axis_;
    double z1 = ray.direction().Dot(axis_);
    double z0 = spos.Dot(axis_);
    auto pvec = ray.direction() - z1*axis_;
    auto qvec = spos - z0*axis_;
    double b2 = pvec.Mag2() - z1*z1*drda_*drda_;
    double b1 = z1*drda_*(drda_*z0 + r0_) - pvec.Dot(qvec);
    double b0 = qvec.Mag2() - pow(drda_*z0 + r0_,2);
    // check if ray is parallel to cone surface
    if(b2 == 0 && b1 != 0){
      retval.onsurface_ = true;
     dist = 0.5*b0/b1;
    } else if (b1*b1 > b0*b2){
       double t1 = b1/b2;
      double t2 = sqrt(b1*b1 - b0*b2)/b2;
      double d1 = t1 + t2;
      double d2 = t1 - t2;
      if(!forwards){
        retval.onsurface_ = true;
        dist = fabs(d1) < fabs(d2) ? d1 : d2;
      } else {
        // closest forwards solution
        if(d1 > 0.0 && d2 > 0.0 ){
          retval.onsurface_ = true;
          dist = std::min(d1,d2);
        } else if(d1 > 0.0) {
          retval.onsurface_ = true;
          dist = d1;
        } else if(d2 > 0.0) {
          retval.onsurface_ = true;
          dist = d2;
        }
      }
    }
    // test that the point is on the cone
    if(retval.onsurface_){
      auto point = ray.position(dist);
      retval.inbounds_ = inBounds(point,tol);
    }
    return retval;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Frustrum const& fru) {
  ost << "Frustrum with center = " << fru.center() << " , axis " << fru.axis()
    << " low radius " << fru.lowRadius()
    << " high radius " << fru.highRadius()
    << " half-length " << fru.halfLength()
    << " half-angle " << fru.halfAngle();
  return ost;
}
