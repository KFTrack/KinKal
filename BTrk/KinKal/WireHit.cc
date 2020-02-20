#include "BTrk/KinKal/WireHit.hh"

namespace KinKal {
  bool WireHit::resid(TPOCABase const& tpoca, RESID& resid, RDer const& dRdDT, double nsigma) const {
    bool retval(true);
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
      float tdrift, tdriftrms, speed;
      d2T().distanceToTime(drift,tdrift, tdriftrms, speed);
      // residual is in time space
      RESID::RVec rvec; rvec(0) = tpoca.dt()-tdrift;
      RESID::RCov rcov; rcov(0,0) = tdriftrms*tdriftrms;  // should include intrinsic measurement error FIXME! 
      resid = RESID(rvec,rcov);
       // fill derviatives FIXME!
    } else {
      // interpret DOCA against the wire directly as the residual

    }
    return retval;
  }
}
