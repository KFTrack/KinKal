#include "BTrk/KinKal/TPOCA.hh"
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/PLine.hh"

using namespace std;

namespace KinKal {
  // specialization between a helix and a line
  template<> void TPOCA<LHelix,PLine>::findPOCA() {
    auto const& lhelix = t1_;
    auto const& pline = t2_;
    // reset poca
    reset();
    // initialize the helix time using the Z position of the line
    double htime = lhelix.time(pline.z0());
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
	status_ = TPOCA::failed;
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
	status_ = TPOCA::failed;
	break;
      }
      doca = sqrt(dd2);      
      niter++;
    }
    // finalize TPOCA
    if(status_ != TPOCA::failed){
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

  //    template<> void TPOCA<LHelix,PLine>::findDPOCA(LHelix::PDer& derivs) {
  //      // find POCA	
  //      findPOCA();
  //      if(usable()){
  //	Vec3 vdoca, ddir;
  //	doca(vdoca);
  //	ddir = vdoca.Unit();
  //	double Factor=1.0; //FIXME!
  //	derivs[LHelix::cx_][0] = -Factor*ddir.X();
  //      }
  //    }

}

