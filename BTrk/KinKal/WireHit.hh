#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/TrkHit.hh"
#include "BTrk/KinKal/D2T.hh"
namespace KinKal {

  template class WireHit : public TrkHit<1> {
    public:
      enum LRAmbig {
      // interpet the TPOCA as a drift residual WRT the wire.  This uses the D2T function and the Cell descriptin
      virtual bool resid(TPOCABase const& tpoca, Resid& resid, RDer const& dRdDT, double nsigma) const override;
      // construct from a D2T relationship
      WireHit(D2T const& d2t) : d2t_(dt2) {}
      // determine if a position is inside the drift cell, within tolerance
      bool inCell(TPOCABase const& tpoca, double nsigma) const = 0;
      // translate a position into the wire coordinate system.
    private:
      double wrad_; // wire radius
      D2T const& d2t_; // distance to time relationship
      LRambig ambig_; // L
  };
}
#endif
