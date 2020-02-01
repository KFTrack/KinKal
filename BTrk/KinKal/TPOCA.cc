#include "BTrk/KinKal/TPOCA.hh"
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/PLine.hh"

using namespace std;

namespace KinKal {
  // specialization between a helix and a line
  template<> void TPOCA<LHelix,PLine>::findPOCA() {
    auto const& lhelix = ttraj0();
    auto const& pline = ttraj1();
    // reset poca
    reset();
    // initialize the helix time using the Z position of the line
    double htime = lhelix.ztime(pline.z0());
    double ltime = pline.t0();
    // use successive linear approximation until desired precision on DOCA is met.
    double doca(1.0e5);
    static const unsigned maxiter=100; // don't allow infinite iteration
    unsigned niter(0);
    // line parameters don't change: compute them once
    Vec3 lpos,ldir;
    pline.pos0(lpos);
    pline.direction(ltime,ldir);
    // helix velocity (scalar) doesn't change
    Vec3 hvel;
    lhelix.velocity(htime,hvel);
    double hv = sqrt(hvel.Mag2());
    while(doca > precision_ && niter < maxiter) {
      // find helix local position and direction
      Vec3 hpos,hdir;
      lhelix.position(htime,hpos);
      lhelix.direction(htime,hdir);
      auto dpos = lpos-hpos;
      // dot products
      double ddot = ldir.Dot(hdir);
      double denom = 1.0 - ddot*ddot;
      // check for parallel)
      if(denom<1.0e-5){
	status_ = TPOCA::pocafailed;
	break;
      }
      double hdd = dpos.Dot(hdir);
      double ldd = dpos.Dot(ldir);
      // compute length from expansion point to POCA and convert to times
      double hlen = (hdd - ldd*ddot)/denom;
      double llen = (hdd*ddot - ldd)/denom;
      htime += hlen/hv; // helix time is iterative
      ltime = pline.t0() + llen/pline.velocity();  // line time is always WRT t0
      // compute DOCA 
      double dd2 = dpos.Mag2() - (hdd*hdd + ldd*ldd + 2.0*hdd*ldd*(1.0 - 2.0*ddot))/denom;
      if(dd2 < 0.0){
	status_ = TPOCA::pocafailed;
	break;
      }
      doca = sqrt(dd2);      
      niter++;
    }
    // finalize TPOCA
    if(status_ != TPOCA::pocafailed){
      if(niter < maxiter)
	status_ = TPOCA::converged;
      else
	status_ = TPOCA::unconverged;
      doca_ = doca;
      // set the positions
      poca_[0].SetE(htime);
      lhelix.position(poca_[0]);
      poca_[1].SetE(ltime);
      pline.position(poca_[1]);
    }
  }

  template<> void TDPOCA<LHelix,PLine,LHelix::npars_>::fillDerivatives() {
    if(!usable())return;
    auto const& lhelix = ttraj0();
    auto const& pline = ttraj1();
// pre-compute some values.
    Vec3 vdoca, ddir, ldir, hdir;
    doca(vdoca);
    ddir = vdoca.Unit();// direction vector along D(POCA) from traj 2 to 1 (line to helix)
    lhelix.direction(poca0().T(),hdir);
    pline.direction(poca1().T(),ldir);

    double hphi = lhelix.phi(poca0().T()); // local azimuth of helix
    double lphi = ldir.Phi(); // line azimuth
    double sineta = sin(lphi);
    double coseta = cos(lphi);
    double dphi = hphi - lphi;
    double sindphi = sin(dphi);
    double cosdphi = cos(dphi);
    double l2 = lhelix.lam()*lhelix.lam();
    double r2 = lhelix.rad()*lhelix.rad();
    double s2 = sindphi*sindphi;
    double denom =  l2 + s2*r2;
    double Factor = lhelix.lam()/sqrt(denom);

    double dx = lhelix.cx() - poca1().X();
    double dy = lhelix.cy() - poca1().Y();
    double ddot = -sineta*dx + coseta*dy;

    dDdP_[LHelix::t0_][0] = 0.0; // no t0 dependence, DOCA is purely geometric
    dDdP_[LHelix::cx_][0] = -Factor*sineta;
    dDdP_[LHelix::cy_][0] = Factor*coseta;
    // components that depend on phi; this ignores terms of magnitude DOCA/Radius
    dDdP_[LHelix::phi0_][0] = Factor*lhelix.rad()*sindphi;
    dDdP_[LHelix::rad_][0] = -Factor*(l2*cosdphi + lhelix.rad()*s2*ddot)/denom;
    // this becomes inaccurate when z=0, but should have no effect on the fit
    dDdP_[LHelix::lam_][0] = -Factor*lhelix.rad()*sindphi*poca0().Z()/l2; 
  }

}

