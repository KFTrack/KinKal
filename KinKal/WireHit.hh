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
      // interpet the TPOCA as a drift residual WRT the wire.  This uses the D2T function and the Cell description
      virtual bool resid(TPOCABase const& tpoca, Residual& resid, RDer const& dRdDT, double nsigma) const override;
      // construct from a D2T relationship
      virtual TLine const& sensorTraj() const override { return wire_; }
      WireHit(TLine const& wire, Context const& context, D2T const& d2t,LRAmbig ambig=null,bool active=true) : THit(context,active), wire_(wire), d2t_(d2t), ambig_(ambig) {}
      virtual ~WireHit(){}
      virtual unsigned nDOF() const override { return 1; }
      // determine if a position is inside the drift cell, within tolerance
      virtual bool inCell(TPOCABase const& tpoca, double nsigma) const = 0;
      LRAmbig ambig() const { return ambig_; }
      D2T const& d2T() const { return d2t_; }
    private:
      TLine wire_; // line representing the local wire position of this hit
      D2T const& d2t_; // distance to time relationship
      LRAmbig ambig_;
  };
}
#endif
