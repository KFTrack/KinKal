#include "KinKal/Geometry/Annulus.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace KinKal {
  bool Annulus::inBounds(VEC3 const& point, double tol) const {
    auto radius = Perp(point - center(),normal());
    return radius > irad_ - tol && radius < orad_ + tol;
  }
}

std::ostream& operator <<(std::ostream& ost, KinKal::Annulus const& ann) {
  KinKal::Plane const& plane = static_cast<KinKal::Plane const&>(ann);
  ost << "Annulus in " << plane << " Inner, Outer radii " << ann.innerRadius() << " , " << ann.outerRadius();
  return ost;
}
