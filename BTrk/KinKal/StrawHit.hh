#ifndef KinKal_StrawHit_hh
#define KinKal_StrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/WireHit.hh"
#include "BTrk/KinKal/Types.hh"
namespace KinKal {

  class StrawHit : public WireHit {
    public:
      virtual bool inCell(TPOCABase const& tpoca, double nsigma) const override;
      virtual void update(TPOCABase const& tpoca) const;
    private:
      double radius_; // straw radius
      Pol2 offset_; // straw offset WRT wire center.  Not yet implemented FIXME!
      // add state for wire length, longitudinal resolution, transverse resolution FIXME!
  };
}
#endif
