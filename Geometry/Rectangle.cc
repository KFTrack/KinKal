#include "KinKal/Geometry/Rectangle.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace KinKal {
  Rectangle::Rectangle(VEC3 const& norm, VEC3 const& center, VEC3 const& uaxis, double uhalflen, double vhalflen) :
    Plane(norm,center) , uhalflen_(uhalflen), vhalflen_(vhalflen), udir_(uaxis.Unit()){
      // check that U is perpendicular
      if(udir_.Dot(normal()) > 1e-10) throw std::invalid_argument("U direction not perpendicular to normal");
      // V direction is implicit
      vdir_ = normal().Cross(udir_);
    }


  bool Rectangle::inBounds(VEC3 const& point, double tol) const {
    auto rvec = point - center();
    double udist = rvec.Dot(udir_);
    double vdist = rvec.Dot(vdir_);
    return fabs(udist) - uhalflen_ < tol && fabs(vdist) - vhalflen_ < tol;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Rectangle const& rect) {
  KinKal::Plane const& plane = static_cast<KinKal::Plane const&>(rect);
  ost << "Rectangle in " << plane << " U direction " << rect.uDirection() << " V direction " << rect.vDirection()
    << " U, V half-lengths " << rect.uHalfLength() << " , " << rect.vHalfLength();
  return ost;
}
