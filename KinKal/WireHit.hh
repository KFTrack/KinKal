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
      virtual void resid(TPocaBase const& tpoca, Residual& resid) const override;
      virtual unsigned nDOF() const override { return 1; }
      // construct from a D2T relationship
      virtual TLine const& sensorTraj() const override { return wire_; }
      WireHit(TLine const& wire, BField const& bfield, D2T const& d2t,double nullvar, LRAmbig ambig=null, bool active=true) : 
	THit(bfield,active), wire_(wire), d2t_(d2t), nullvar_(nullvar), ambig_(ambig) {}
      virtual ~WireHit(){}
      LRAmbig ambig() const { return ambig_; }
      D2T const& d2T() const { return d2t_; }
    private:
      TLine wire_; // local linear approximation to the wire of this hit.  The range describes the active wire length
      D2T const& d2t_; // distance to time relationship for drift in this cell
      double nullvar_; // variance of the error in space for null ambiguity
      LRAmbig ambig_; // current ambiguity assignment: can change during a fit
  };
}
#endif
