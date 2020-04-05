#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/THit.hh"
#include "KinKal/D2T.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/LRAmbig.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class WireHit : public THit<KTRAJ> {
    public:
      using THit<KTRAJ>::bfield;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef TPoca<PKTRAJ,TLine> TPOCA;
      typedef typename KTRAJ::PDER PDER; 
      // THit interface overrrides
      virtual void resid(Residual& resid) const override;
      virtual unsigned nDOF() const override { return 1; }
      virtual void update(PKTRAJ const& pktraj) override;
      virtual TPOCA const& poca() const override { return tpoca_; }
      virtual PDER const& dDdP() const override { return tpoca_.dDdP(); }
      virtual PDER const& dTdP() const override { return tpoca_.dTdP(); }
// construct from a D2T relationship
      TLine const& wireTraj() const { return wire_; }
      WireHit(BField const& bfield, PKTRAJ const& pktraj, TLine const& wire, D2T const& d2t,double nullvar, LRAmbig ambig=LRAmbig::null, bool active=true) : 
	THit<KTRAJ>(bfield,active), wire_(wire), tpoca_(pktraj,wire_), d2t_(d2t), nullvar_(nullvar), ambig_(ambig) {}
      virtual ~WireHit(){}
      LRAmbig ambig() const { return ambig_; }
      D2T const& d2T() const { return d2t_; }
    private:
      TLine wire_; // local linear approximation to the wire of this hit.  The range describes the active wire length
      TPOCA tpoca_;
      D2T const& d2t_; // distance to time relationship for drift in this cell
      double nullvar_; // variance of the error in space for null ambiguity
      LRAmbig ambig_; // current ambiguity assignment: can change during a fit
  };

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    tpoca_ = TPOCA(pktraj,wire_);
    // could update for sag, etc. FIXME!
  }

  template <class KTRAJ> void WireHit<KTRAJ>::resid( Residual& resid) const {
    if(tpoca_.usable()){
      // translate TPOCA to residual
      if(ambig_ != LRAmbig::null){ 
	auto iambig = static_cast<std::underlying_type<LRAmbig>::type>(ambig_);
	// convert DOCA to wire-local polar coordinates.  This defines azimuth WRT the B field for ExB effects
	float rho = tpoca_.doca()*iambig; // this is allowed to go negative
	Vec3 bvec;
	bfield().fieldVect(bvec,tpoca_.particlePoca().Vect());
	auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BField
	Vec3 dvec;
	tpoca_.delta(dvec);
	float phi = asin(float(dvec.Unit().Dot(pdir)));
	Pol2 drift(rho, phi);
	float tdrift, tdvar, vdrift;
	d2T().distanceToTime(drift, tdrift, tdvar, vdrift);
	// should add intrinsic measurement effects to totvar FIXME!
	// should add annealing temperature to totvar FIXME!
	float totvar = tdvar;
	// residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
	resid = Residual(tpoca_.deltaT()-tdrift,totvar,iambig/vdrift,1.0);
      } else {
	// interpret DOCA against the wire directly as the residual.  There is no direct time dependence in this case
	// residual is in space, so unit dependendence on distance, none on time
	resid = Residual(-tpoca_.doca(),nullvar_,1.0,0.0);
      }
    } else
      throw std::runtime_error("POCA failure");
  }

}
#endif
