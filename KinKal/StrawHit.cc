#include "KinKal/StrawHit.hh"
#include <algorithm>

namespace KinKal {
  float StrawHit::inRange(TPocaBase const& tpoca) const {
    Vec3 tpos; tpoca.delta(tpos); // transverse to wire
    // correct for axis offset FIXME!
   float retval = std::max((double)0.0,sqrt(tpos.Mag2()) -smat_.strawRadius()) /0.05; // error should come from hit FIXME!
   // should also check length along straw FIXME!
   return retval;
  }
  void StrawHit::update(TPocaBase const& tpoca) const {
    //FIXME!
  }
}
