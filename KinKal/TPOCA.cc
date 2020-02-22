#include "KinKal/TPOCA.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/PKTraj.hh"
// specializations for TPOCA
using namespace std;
namespace KinKal {
  // specialization between a looping helix and a line
  template<> TPOCA<LHelix,TLine>::TPOCA(LHelix const& lhelix, TLine const& tline, double precision) : TPOCABase(lhelix,tline,precision)  { 
    // reset status
    reset();
    double htime,ltime;
    // initialize the helix time using the Z position of the line
    // this will fail if the line is long and parallel to the z axis as there can be multiple solutions FIXME!
    htime = lhelix.ztime(tline.z0());
    ltime = tline.t0();
    // use successive linear approximation until desired precision on DOCA is met.
    double ddoca(1.0e5);
    double doca(0.0);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter(0);
    // helix speed doesn't change
    double hspeed = lhelix.speed(lhelix.t0());
    Vec3 hdir;
    while(fabs(ddoca) > precision_ && niter++ < maxiter) {
      // find helix local position and direction
      Vec3 hpos;
      lhelix.position(htime,hpos);
      lhelix.direction(htime,hdir);
      auto dpos = tline.pos0()-hpos;
      // dot products
      double ddot = tline.dir().Dot(hdir);
      double denom = 1.0 - ddot*ddot;
      // check for parallel)
      if(denom<1.0e-5){
	status_ = TPOCA::pocafailed;
	break;
      }
      double hdd = dpos.Dot(hdir);
      double ldd = dpos.Dot(tline.dir());
      // compute length from expansion point to POCA and convert to times
      double hlen = (hdd - ldd*ddot)/denom;
      double llen = (hdd*ddot - ldd)/denom;
      htime += hlen/hspeed; // helix time is iterative
      ltime = tline.t0() + llen/tline.speed(ltime);  // line time is always WRT t0
      // compute DOCA
      lhelix.position(htime,hpos);
      Vec3 lpos;
      tline.position(ltime,lpos);
      double dd2 = (hpos-lpos).Mag2();
      if(dd2 < 0.0 ){
	status_ = TPOCA::pocafailed;
	break;
      }
      ddoca = doca;
      doca = sqrt(dd2);
      ddoca -= doca;
      if(isnan(ddoca)){
	status_ = TPOCA::pocafailed;
	break;
      }
    }
    // finalize TPOCA
    if(status_ != TPOCA::pocafailed){
      if(niter < maxiter)
	status_ = TPOCA::converged;
      else
	status_ = TPOCA::unconverged;
        // set the positions
      poca_[0].SetE(htime);
      lhelix.position(poca_[0]);
      poca_[1].SetE(ltime);
      tline.position(poca_[1]);
      // sign doca by angular momentum projected onto difference vector
      double lsign = tline.dir().Cross(hdir).Dot(poca_[1].Vect()-poca_[0].Vect());
      doca_ = copysign(doca,lsign);
    }
  }

  template<> TDPOCA<LHelix,TLine>::TDPOCA(TPOCA<LHelix,TLine> const& tpoca) : TPOCA<LHelix,TLine>(tpoca) {
  // check status of POCA calculation from base class
    if(usable()){
      auto const& lhelix = ttraj0();
      auto const& tline = ttraj1();
      // pre-compute some values.
      Vec3 vdoca, ddir, hdir;
      delta(vdoca);
      ddir = vdoca.Unit();// direction vector along D(POCA) from traj 2 to 1 (line to helix)
      lhelix.direction(poca0().T(),hdir);

      double hphi = lhelix.phi(poca0().T()); // local azimuth of helix
      double lphi = tline.dir().Phi(); // line azimuth
      double sineta = sin(lphi);
      double coseta = cos(lphi);
      double dphi = hphi - lphi;
      double sindphi = sin(dphi);
      double cosdphi = cos(dphi);
      double l2 = lhelix.lam()*lhelix.lam();
      double r2 = lhelix.rad()*lhelix.rad();
      double s2 = sindphi*sindphi;
      double denom =  l2 + s2*r2;
      double Factor = fabs(lhelix.lam())/sqrt(denom);

      double dx = lhelix.cx() - poca1().X();
      double dy = lhelix.cy() - poca1().Y();
      double ddot = -sineta*dx + coseta*dy;

      // no t0 dependence, DOCA is purely geometric
      dDdP_[LHelix::cx_] = -Factor*sineta;
      dDdP_[LHelix::cy_] = Factor*coseta;
      // components that depend on phi; this ignores terms of magnitude DOCA/Radius
      dDdP_[LHelix::phi0_] = Factor*lhelix.rad()*sindphi;
      dDdP_[LHelix::rad_] = -Factor*(l2*cosdphi + lhelix.rad()*s2*ddot)/denom;
      // this becomes inaccurate when z=0, but should have no effect on the fit
      dDdP_[LHelix::lam_] = -Factor*lhelix.rad()*sindphi*poca0().Z()/l2; 

      // no spatial dependence, DT is purely temporal
      dTdP_[LHelix::t0_] = 1.0; // time is 100% correlated
    }
  }

  template<> TDPOCA<LHelix,TLine>::TDPOCA(LHelix const& lhelix, TLine const& tline, double precision) : TDPOCA<LHelix,TLine>(TPOCA<LHelix,TLine>(lhelix,tline,precision)){}

  // specialization between a piecewise LHelix and a line
  typedef PKTraj<LHelix> PLHelix;
  template<> TPOCA<PLHelix,TLine>::TPOCA(PLHelix const& phelix, TLine const& tline, double precision) : TPOCABase(phelix,tline,precision)  {
// iteratively find the nearest piece, and POCA for that piece.  Start in the middle
    size_t oldindex= phelix.pieces().size();
    size_t index = size_t(rint(oldindex/2.0));
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    status_ = converged; 
    unsigned niter=0;
    while(status_ == converged && niter++ < maxiter && index != oldindex){
    // call down to LHelix TPOCA
      LHelix const& piece = phelix.pieces()[index];
      TPOCA<LHelix,TLine> tpoca(piece,tline,precision);
    // copy over state
      status_ = tpoca.status();
      poca_[0] = tpoca.poca0();
      poca_[1] = tpoca.poca1();
      doca_ = tpoca.doca();
      // prepare for the next iteration
      oldindex = index;
      index = phelix.nearestIndex(tpoca.poca0().T());
    }
    if(status_ == converged && niter >= maxiter) status_ = unconverged;
  }

  // same for TDPOCA
  template<> TDPOCA<PLHelix,TLine>::TDPOCA(TPOCA<PLHelix,TLine> const& tpoca) : TPOCA<PLHelix,TLine>(tpoca) {
    if(usable()){
      // call down to piece TDPOCA.  Unfortunately this re-computes TPOCA FIXME
      TDPOCA<LHelix,TLine> tdpoca(ttraj0().nearestPiece(poca0().T()),ttraj1(),precision());
      // copy over state
      dDdP_ = tdpoca.dDdP();
      dTdP_ = tdpoca.dTdP();
    }
  }
  template<> TDPOCA<PLHelix,TLine>::TDPOCA(PLHelix const& phelix, TLine const& tline, double precision) :
    TPOCA<PLHelix,TLine>(TPOCA<PLHelix,TLine>(phelix,tline,precision)) { }
}
