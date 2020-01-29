#include "BTrk/KinKal/TPOCA.hh"
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/PWire.hh"

namespace KinKal {
 template<> void FindTPOCA(LHelix const& lhelix, PWire const& pwire, TPOCA& tpoca) {
  // compute wire azimuth from direction; input time doesn't matter
  Vec3 wdir, wpos;
  pwire.direction(pwire.t0(),wdir);
  // position at t0
  pwire.pos0(wpos);
  // ignore slope; initial z value is z0 of the wire
  // compute helix azimuth angle from the hit z position
//  double tprop = lhelix.time(wpos.Z());
  double hphi = wpos.Z()/lhelix.lam() + lhelix.phi0();



  double wphi = atan2(wdir.Y(),wdir.X());
  // denomiator
  // angle difference
  double phibar = hphi-wphi;
  double sinphibar = sin(phibar);
  double cosphibar = cos(phibar);
  double sineta = sin(wphi);
  double coseta = cos(wphi);
  double l2 = lhelix.lam()*lhelix.lam();
  double r2 = lhelix.rad()*lhelix.rad();
  double s2 = sinphibar*sinphibar;
  double denom =  l2 + s2*r2;
  double Factor = lhelix.lam()/sqrt(denom);
  double dx = lhelix.cx() - wpos.X(); 
  double dy = lhelix.cy() - wpos.Y();
  double ddot = -sineta*dx + coseta*dy;
  // doca; this approximation ignores 2nd order terms.  It becomes
  // innacurate when tandip is small and z large
  double doca = -Factor*(lhelix.rad()*cosphibar - ddot );

  tpoca.p1_.SetPx(doca); // junk FIXME!

  // compute time for the different trajs FIXME!

 }

}
