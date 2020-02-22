#include "KinKal/WireHit.hh"

namespace KinKal {
  void WireHit::resid(TPOCABase const& tpoca, Residual& resid) const {
    Residual::RVec rvec;
    Residual::RCov rcov;
    if(ambig_ != null){ 
      // convert DOCA to wire-local polar coordinates.  This defines azimuth WRT the B field for ExB effects
      float rho = tpoca.doca()*ambig_; // this is allowed to go negative
      Vec3 bvec;
      context().bField().fieldVect(bvec,tpoca.poca1().Vect());
      auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BField
      Vec3 dvec;
      tpoca.delta(dvec);
      float phi = asin(float(dvec.Unit().Dot(pdir)));
      Pol2 drift(rho, phi);
      float tdrift, tdriftrms, vdrift;
      d2T().distanceToTime(drift, tdrift, tdriftrms, vdrift);
      // residual is in time
      rvec(0) = tpoca.dt()-tdrift; // measurement - prediction
      rcov(0,0) = tdriftrms*tdriftrms;  // should include intrinsic measurement error and velocity dependence FIXME! 
      // Clonstruct the residual object: dRdt is just the local drift speed
      resid = Residual(rvec,rcov,vdrift);
    } else {
      // interpret DOCA against the wire directly as the residual.  There is no direct time dependence in this case
      rvec(0) = tpoca.doca();
      rcov(0,0) = nullms_;
      Residual(rvec,rcov,0.0);
    }
  }
}
