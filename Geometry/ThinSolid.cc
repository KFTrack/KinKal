#include "KinKal/Geometry/ThinSolid.hh"
namespace KinKal {
  bool ThinSolid::isInside(VEC3 const& point) const {
    return fabs(surf_->distance(point)) < halfthick_;
  }
}
