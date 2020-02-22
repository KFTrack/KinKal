#include "KinKal/StrawHit.hh"
#include <algorithm>

namespace KinKal {
  float StrawHit::inRange(TPOCABase const& tpoca) const {
    Vec3 tpos; tpoca.delta(tpos); // transverse to wire
    // correct for axis offset FIXME!
   float retval = std::max((double)0.0,sqrt(tpos.Mag2()) -radius_) /0.05; // error should come from hit FIXME!
   // should also check length along straw FIXME!
   return retval;
  }
  void StrawHit::update(TPOCABase const& tpoca) const {
    //FIXME!
  }
}
