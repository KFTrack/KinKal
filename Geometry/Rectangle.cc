#include "KinKal/Geometry/Rectangle.hh"
#include "Math/VectorUtil.h"
#include <exception>
using namespace ROOT::Math::VectorUtil;
namespace KinKal {

  Rectangle::Rectangle(VEC3 const& norm, VEC3 const& udir, VEC3 const& center, double uhalflen, double vhalflen) :
    Plane(norm,udir,center),
    uhalflen_(uhalflen), vhalflen_(vhalflen) {
      if(uhalflen <= 0.0 || vhalflen <= 0.0 ) throw std::invalid_argument("Invalid Rectangle arguments");

    }
  bool Rectangle::inBounds(VEC3 const& point, double tol) const {
    auto rvec = point - center();
    double udist = rvec.Dot(uDirection());
    double vdist = rvec.Dot(vDirection());
    return fabs(udist) - uhalflen_ < tol && fabs(vdist) - vhalflen_ < tol;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Rectangle const& rect) {
  KinKal::Plane const& plane = static_cast<KinKal::Plane const&>(rect);
  ost << "Rectangle in " << plane << " U, V half-lengths " << rect.uHalfLength() << " , " << rect.vHalfLength();
  return ost;
}
