#ifndef KinKal_DOCAWireHitUpdater_hh
#define KinKal_DOCAWireHitUpdater_hh
// simple WireHit updater based on DOCA
#include "KinKal/Examples/WireHitStructs.hh"
namespace KinKal {
  class DOCAWireHitUpdater {
    public:
      DOCAWireHitUpdater(double mindoca,double maxdoca ) : mindoca_(mindoca), maxdoca_(maxdoca) {}
      // define the state given the (presumably unbiased) distance of closest approach
      void updateWireHitState(WireHitState &state, double doca) const;
      double minDOCA() const { return mindoca_; }
      double maxDOCA() const { return maxdoca_; }
    private:
      double mindoca_; // minimum DOCA value to sign LR ambiguity
      double maxdoca_; // maximum DOCA to still use a hit
  };

  void DOCAWireHitUpdater::updateWireHitState(WireHitState &state, double doca) const {
    if(fabs(doca) > maxdoca_ ) {
      state.state_ = WireHitState::inactive; // disable the hit if it's an outlier
    } else if(fabs(doca) > mindoca_ ) {
      state.state_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
    } else {
      // too close to the wire: don't try to disambiguate LR sign
      state.state_ = WireHitState::null;
    }
  }
}
#endif
