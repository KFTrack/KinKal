#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement.  Implemented using TPOCA between the particle traj and the wire
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/THit.hh"
#include "KinKal/D2T.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/LRAmbig.hh"
#include "KinKal/BField.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class WireHit : public THit<KTRAJ> {
    public:
      typedef THit<KTRAJ> THIT;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef TPoca<PKTRAJ,TLine> TPOCA;
      typedef Residual<KTRAJ> RESIDUAL;
      typedef DXing<KTRAJ> DXING;
      typedef std::shared_ptr<DXING> DXINGPTR;
      typedef typename KTRAJ::PDER PDER;
      // THit interface overrrides
      virtual void resid(PKTRAJ const& pktraj, RESIDUAL& resid) const override;
      virtual unsigned nDOF() const override { return 1; }
// construct from a D2T relationship; BField is needed to compute ExB effects
      TLine const& wireTraj() const { return wire_; }
      WireHit(DXINGPTR const& dxing, BField const& bfield, TLine const& wire, D2T const& d2t, float nullvar, LRAmbig ambig=LRAmbig::null, bool active=true) : 
	THIT(dxing, active), wire_(wire), d2t_(d2t), nullvar_(nullvar), ambig_(ambig), bfield_(bfield) {}
      virtual ~WireHit(){}
      LRAmbig ambig() const { return ambig_; }
      D2T const& d2T() const { return d2t_; }
    private:
      TLine wire_; // local linear approximation to the wire of this hit.  The range describes the active wire length
      D2T const& d2t_; // distance to time relationship for drift in this cell
      float nullvar_; // variance of the error in space for null ambiguity
      LRAmbig ambig_; // current ambiguity assignment: can change during a fit
      BField const& bfield_;
  };

  template <class KTRAJ> void WireHit<KTRAJ>::resid(PKTRAJ const& pktraj, RESIDUAL& resid) const {
    // compute TPOCA.  Measurement time doesn't provide a good hint
    TPOCA tpoca(pktraj,wire_);
    if(tpoca.usable()){
      // translate TPOCA to residual
      if(ambig_ != LRAmbig::null){ 
	auto iambig = static_cast<std::underlying_type<LRAmbig>::type>(ambig_);
	// convert DOCA to wire-local polar coordinates.  This defines azimuth WRT the B field for ExB effects
	float rho = tpoca.doca()*iambig; // this is allowed to go negative
	Vec3 bvec;
	bfield_.fieldVect(bvec,tpoca.particlePoca().Vect());
	auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BField
	Vec3 dvec;
	tpoca.delta(dvec);
	float phi = asin(float(dvec.Unit().Dot(pdir)));
	Pol2 drift(rho, phi);
	float tdrift, tdvar, vdrift;
	d2T().distanceToTime(drift, tdrift, tdvar, vdrift);
	// should add intrinsic measurement effects to tdvar FIXME!
	// should add annealing temperature to totvar FIXME!
	// residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
	PDER dRdP = tpoca.dDdP()*iambig/vdrift + tpoca.dTdP(); 
	resid = RESIDUAL(tpoca.particleToca(),tpoca.deltaT()-tdrift,tdvar,dRdP);
      } else {
	// interpret DOCA against the wire directly as the residual.  There is no direct time dependence in this case
	// residual is in space, so unit dependendence on distance, none on time
	resid = RESIDUAL(tpoca.particleToca(),-tpoca.doca(),nullvar_,tpoca.dDdP());
      }
    } else
      throw std::runtime_error("POCA failure");
  }
}
#endif
