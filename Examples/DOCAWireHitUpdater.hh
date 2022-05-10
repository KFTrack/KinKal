#ifndef KinKal_DOCAWireHitUpdater_hh
#define KinKal_DOCAWireHitUpdater_hh
// simple WireHit updater based on DOCA
#include "KinKal/Detector/WireHitStructs.hh"
namespace KinKal {
  class DOCAWireHitUpdater {
    public:
      DOCAWireHitUpdater(double mindoca,double maxdoca ) : mindoca_(mindoca), maxdoca_(maxdoca) {}
      // define the state given the (presumably unbiased) distance of closest approach
      WireHitState wireHitState(double doca) const;
      double minDOCA() const { return mindoca_; }
      double maxDOCA() const { return maxdoca_; }
    private:
      double mindoca_; // minimum DOCA value to sign LR ambiguity
      double maxdoca_; // maximum DOCA to still use a hit
  };

  WireHitState DOCAWireHitUpdater::wireHitState(double doca) const {
    WireHitState state;
    if(fabs(doca) > maxdoca_ ) {
      state = WireHitState::inactive; // disable the hit if it's an outlier
    } else if(fabs(doca) > mindoca_ ) {
      state = doca > 0.0 ? WireHitState::right : WireHitState::left;
    } else {
      // too close to the wire: don't try to disambiguate LR sign
      state = WireHitState::null;
    }
    return state;
  }
}
#endif
