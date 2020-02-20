#include "BTrk/KinKal/StrawHit.hh"
namespace KinKal {
  bool StrawHit::inCell(TPOCABase const& tpoca, double nsigma) const {
    Vec3 tpos; tpoca.delta(tpos); // transverse to wire
    // correct for axis offset FIXME!
   return sqrt(tpos.Mag2()) < radius_ + nsigma*0.05; // error should come from hit FIXME!
   // should also check length along straw FIXME!
  }
  void StrawHit::update(TPOCABase const& tpoca) const {
    //FIXME!
  }
}
