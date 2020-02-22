#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/THit.hh"
#include "KinKal/D2T.hh"
#include "KinKal/TLine.hh"
namespace KinKal {

  class WireHit : public THit {
    public:
    // define left-right ambiguity, = sign of angular momentum of particle trajectory WRT hit trajectory
    // Null means to ignore drift information and constrain to the wire
      enum LRAmbig { left=-1, null=0, right=1}; 
      // THit interface overrrides
      virtual void resid(TPOCABase const& tpoca, Residual& resid) const override;
      virtual unsigned nDOF() const override { return 1; }
      // construct from a D2T relationship
      virtual TLine const& sensorTraj() const override { return wire_; }
      WireHit(TLine const& wire, Context const& context, D2T const& d2t,double nullms, LRAmbig ambig=null, bool active=true) : 
	THit(context,active), wire_(wire), d2t_(d2t), nullms_(nullms), ambig_(ambig) {}
      virtual ~WireHit(){}
      LRAmbig ambig() const { return ambig_; }
      D2T const& d2T() const { return d2t_; }
    private:
      TLine wire_; // line representing the local wire position of this hit.  The range describes the wire length
      D2T const& d2t_; // distance to time relationship
      double nullms_; // mean square (not root!) of the error in space for null ambiguity
      LRAmbig ambig_; // current ambiguity assignment: can change during a fit
  };
}
#endif
