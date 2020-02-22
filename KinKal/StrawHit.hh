#ifndef KinKal_StrawHit_hh
#define KinKal_StrawHit_hh
//
//  class representing a straw sensor measurement.  It assumes a (possibly displaced)
//  circular outer cathode locally parallel to the wire.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/WireHit.hh"
namespace KinKal {
  class StrawHit : public WireHit {
    public:
      StrawHit(TLine const& straj, Context const& context, D2T const& d2t,double radius, LRAmbig ambig=null,bool active=true) : 
	WireHit(straj,context,d2t,radius*0.289,ambig,active), radius_(radius) {} // null error set to radius/sqrt(12) FIXME!
      virtual float inRange(TPOCABase const& tpoca) const override;
      virtual void update(TPOCABase const& tpoca) const override;
      virtual ~StrawHit(){}
    private:
      double radius_; // straw radius
      Pol2 offset_; // straw offset WRT wire center.  Not yet implemented FIXME!
      // add state for longitudinal resolution, transverse resolution FIXME!
  };
}
#endif
