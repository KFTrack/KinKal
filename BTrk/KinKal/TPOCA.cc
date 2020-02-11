#include "BTrk/KinKal/TPOCA.hh"
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/PLine.hh"
// specializations for TPOCA
using namespace std;
namespace KinKal {

  // specialization between a helix and a line
  template<> TPOCA<LHelix,PLine>::TPOCA(LHelix const& lhelix, PLine const& pline, double precision) : TPOCABase(lhelix,pline,precision)  { 
    // reset status
    reset();
    // initialize the helix time using the Z position of the line
    double htime = lhelix.ztime(pline.z0());
    double ltime = pline.t0();
    // use successive linear approximation until desired precision on DOCA is met.
    double ddoca(1.0e5);
    double doca(0.0);
    static const unsigned maxiter=100; // don't allow infinite iteration
    unsigned niter(0);
    // line parameters don't change: compute them once
    Vec3 lpos0,ldir;
    pline.pos0(lpos0);
    pline.direction(ltime,ldir);
    // helix velocity (scalar) doesn't change
    Vec3 hvel;
    lhelix.velocity(htime,hvel);
    double hv = sqrt(hvel.Mag2());
    while(fabs(ddoca) > precision_ && niter < maxiter) {
      // find helix local position and direction
      Vec3 hpos,hdir;
      lhelix.position(htime,hpos);
      lhelix.direction(htime,hdir);
      auto dpos = lpos0-hpos;
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
      lhelix.position(htime,hpos);
      Vec3 lpos;
      pline.position(ltime,lpos);
      double dd2 = (hpos-lpos).Mag2();
      if(dd2 < 0.0){
	status_ = TPOCA::pocafailed;
	break;
      }
      ddoca = doca;
      doca = sqrt(dd2);
      ddoca -= doca;
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

  template<> TDPOCA<LHelix,PLine>::TDPOCA(TPOCA<LHelix,PLine> const& tpoca) : TPOCA<LHelix,PLine>(tpoca) {
  // check status of POCA calculation from base class
    if(usable()){
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

  template<> TDPOCA<LHelix,PLine>::TDPOCA(LHelix const& lhelix, PLine const& pline, double precision) : TDPOCA<LHelix,PLine>(TPOCA<LHelix,PLine>(lhelix,pline,precision)){}

}
