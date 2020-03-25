#include "KinKal/WireHit.hh"

namespace KinKal {
  void WireHit::resid(TPocaBase const& tpoca, Residual& resid) const {
    if(ambig_ != null){ 
      // convert DOCA to wire-local polar coordinates.  This defines azimuth WRT the B field for ExB effects
      float rho = tpoca.doca()*ambig_; // this is allowed to go negative
      Vec3 bvec;
      bfield().fieldVect(bvec,tpoca.poca1().Vect());
      auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BField
      Vec3 dvec;
      tpoca.delta(dvec);
      float phi = asin(float(dvec.Unit().Dot(pdir)));
      Pol2 drift(rho, phi);
      float tdrift, tdvar, vdrift;
      d2T().distanceToTime(drift, tdrift, tdvar, vdrift);
      // should add intrinsic measurement effects to tdvar FIXME!
      float totvar = tdvar;
      // residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
      resid = Residual(tpoca.deltaT()-tdrift,totvar,ambig_/vdrift,1.0);
    } else {
      // interpret DOCA against the wire directly as the residual.  There is no direct time dependence in this case
      // residual is in space, so unit dependendence on distance, none on time
      resid = Residual(-tpoca.doca(),nullvar_,1.0,0.0);
    }
  }
}
