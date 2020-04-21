#include "KinKal/TPoca.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/PKTraj.hh"
#include <limits>
// specializations for TPoca
using namespace std;
namespace KinKal {
  // specialization between a looping helix and a line
  template<> TPoca<LHelix,TLine>::TPoca(LHelix const& lhelix, TLine const& tline, TPocaHint const& hint, float precision) : TPocaBase(precision),ktraj_(&lhelix), straj_(&tline) {
    // reset status
    reset();
    float htoca,stoca;
    // initialize the helix time using hints, if available.  If not, use the Z of the line
    if(hint.particleHint_)
      htoca = hint.particleToca_;
    else
      htoca = lhelix.ztime(tline.z0());
    // similar for line; this shouldn't matter, since the solution is linear
    if(hint.sensorHint_)
      stoca = hint.sensorToca_;
    else
      stoca = tline.t0();
    // use successive linear approximation until desired precision on DOCA is met.
    float dptoca(std::numeric_limits<float>::max()), dstoca(std::numeric_limits<float>::max());
    float doca(0.0);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter(0);
    // helix speed doesn't change
    float hspeed = lhelix.speed(lhelix.t0());
    Vec3 hdir;
    // iterate until change in TOCA is less than precision
    while((fabs(dptoca) > precision_ || fabs(dstoca) > precision_) && niter++ < maxiter) {
      // find helix local position and direction
      Vec3 hpos;
      lhelix.position(htoca,hpos);
      lhelix.direction(htoca,hdir);
      auto dpos = tline.pos0()-hpos;
      // dot products
      float ddot = tline.dir().Dot(hdir);
      float denom = 1.0 - ddot*ddot;
      // check for parallel)
      if(denom<1.0e-5){
	status_ = pocafailed;
	break;
      }
      float hdd = dpos.Dot(hdir);
      float ldd = dpos.Dot(tline.dir());
      // compute length from expansion point to POCA and convert to times
      dptoca = (hdd - ldd*ddot)/(denom*hspeed);
      dstoca = tline.t0() + (hdd*ddot - ldd)/(denom*tline.speed(stoca)) - stoca;
      htoca += dptoca; // helix time is iterative
      stoca += dstoca; // line time is always WRT t0, since it uses p0
      // compute DOCA
      lhelix.position(htoca,hpos);
      Vec3 lpos;
      tline.position(stoca,lpos);
      float dd2 = (hpos-lpos).Mag2();
      if(dd2 < 0.0 ){
	status_ = pocafailed;
	break;
      }
      doca = sqrt(dd2);
      // update convergence test
      if(isnan(doca)){
	status_ = pocafailed;
	break;
      }
    }
    // if successfull, finalize TPoca
    if(status_ != pocafailed){
      if(niter < maxiter)
        status_ = TPoca::converged;
      else
        status_ = TPoca::unconverged;
        // set the positions
      partPoca_.SetE(htoca);
      lhelix.position(partPoca_);
      sensPoca_.SetE(stoca);
      tline.position(sensPoca_);
      // sign doca by angular momentum projected onto difference vector
      float lsign = tline.dir().Cross(hdir).Dot(sensPoca_.Vect()-partPoca_.Vect());
      doca_ = copysign(doca,lsign);

      // pre-compute some values needed for the derivative calculations
      Vec3 vdoca, ddir, hdir;
      delta(vdoca);
      ddir = vdoca.Unit();// direction vector along D(POCA) from traj 2 to 1 (line to helix)
      lhelix.direction(particlePoca().T(),hdir);

      float hphi = lhelix.phi(particlePoca().T()); // local azimuth of helix
      float lphi = tline.dir().Phi(); // line azimuth
      float sineta = sin(lphi);
      float coseta = cos(lphi);
      float dphi = hphi - lphi;
      float sindphi = sin(dphi);
      float cosdphi = cos(dphi);
      float l2 = lhelix.lam()*lhelix.lam();
      float r2 = lhelix.rad()*lhelix.rad();
      float s2 = sindphi*sindphi;
      float denom =  l2 + s2*r2;
      float Factor = -copysign(1.0,lhelix.charge())*lhelix.lam()/sqrt(denom);

      float dx = lhelix.cx() - sensorPoca().X();
      float dy = lhelix.cy() - sensorPoca().Y();
      float ddot = -sineta*dx + coseta*dy;

      // no t0 dependence, DOCA is purely geometric
      dDdP_[LHelix::cx_] = -Factor*sineta;
      dDdP_[LHelix::cy_] = Factor*coseta;
      // components that depend on phi; this ignores terms of magnitude DOCA/Radius
      dDdP_[LHelix::phi0_] = Factor*lhelix.rad()*sindphi;
      dDdP_[LHelix::rad_] = -Factor*(l2*cosdphi + lhelix.rad()*s2*ddot)/denom;
      // this becomes inaccurate when z=0, but should have no effect on the fit
      dDdP_[LHelix::lam_] = -Factor*lhelix.rad()*sindphi*particlePoca().Z()/l2; 

      // no spatial dependence, DT is purely temporal
      dTdP_[LHelix::t0_] = 1.0; // time is 100% correlated
      // propagate parameter covariance to variance on doca and toca
      docavar_ = ROOT::Math::Similarity(dDdP(),lhelix.params().covariance());
      tocavar_ = ROOT::Math::Similarity(dTdP(),lhelix.params().covariance());
      // dot product between directions at POCA
      Vec3 pdir, sdir;
      lhelix.direction(particleToca(),pdir);
      tline.direction(sensorToca(),sdir);
      ddot_ = pdir.Dot(sdir);
    }
  }

  // specialization between a piecewise LHelix and a line
  typedef PKTraj<LHelix> PLHELIX;
  template<> TPoca<PLHELIX,TLine>::TPoca(PLHELIX const& phelix, TLine const& tline, TPocaHint const& hint, float precision) : TPocaBase(precision), ktraj_(&phelix), straj_(&tline)  {
    // iteratively find the nearest piece, and POCA for that piece.  Start at hints if availalble, otherwise the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter=0;
    size_t oldindex= phelix.pieces().size();
    size_t index;
    if(hint.particleHint_)
      index = phelix.nearestIndex(hint.particleToca_);
    else
      index = size_t(rint(oldindex/2.0));
    status_ = converged; 
    while(status_ == converged && niter++ < maxiter && index != oldindex){
      // call down to LHelix TPoca
      // prepare for the next iteration
      LHelix const& piece = phelix.pieces()[index];
      TPoca<LHelix,TLine> tpoca(piece,tline,hint,precision);
      status_ = tpoca.status();
      if(tpoca.usable()){
	// copy over the rest of the state
	partPoca_ = tpoca.particlePoca();
	sensPoca_ = tpoca.sensorPoca();
	doca_ = tpoca.doca();
	dDdP_ = tpoca.dDdP();
	dTdP_ = tpoca.dTdP();
	docavar_ = tpoca.docaVar();
	tocavar_ = tpoca.tocaVar();
	ddot_ = tpoca.dirDot();
      }
      oldindex = index;
      index = phelix.nearestIndex(tpoca.particlePoca().T());
    }
    if(status_ == converged && niter >= maxiter) status_ = unconverged;
  }

}
