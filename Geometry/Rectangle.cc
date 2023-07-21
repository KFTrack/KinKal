#include "KinKal/Geometry/Rectangle.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace KinKal {

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
