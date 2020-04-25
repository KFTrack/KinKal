#include "KinKal/TPoca.hh"
#include "KinKal/IPHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/PKTraj.hh"
// specializations for TPoca
using namespace std;
namespace KinKal {
  // specialization between a looping helix and a line
  template<> TPoca<IPHelix,TLine>::TPoca(IPHelix const& lhelix, TLine const& tline, double precision) : TPocaBase(precision),ktraj_(&lhelix), straj_(&tline) {
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
        status_ = pocafailed;
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
        status_ = pocafailed;
        break;
      }
      ddoca = doca;

      doca = sqrt(dd2);
      ddoca -= doca;
      if(isnan(ddoca)){
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
      partPoca_.SetE(htime);
      lhelix.position(partPoca_);
      sensPoca_.SetE(ltime);
      tline.position(sensPoca_);
      // sign doca by angular momentum projected onto difference vector
      double lsign = tline.dir().Cross(hdir).Dot(sensPoca_.Vect()-partPoca_.Vect());
      doca_ = copysign(doca,lsign);

      // pre-compute some values needed for the derivative calculations
      Vec3 vdoca, ddir, hdir;
      delta(vdoca);
      ddir = vdoca.Unit();// direction vector along D(POCA) from traj 2 to 1 (line to helix)
      lhelix.direction(particlePoca().T(),hdir);

      // no t0 dependence, DOCA is purely geometric
      // I need to calculate these derivatives FIXME!
      dDdP_[IPHelix::d0_] = 0;
      dDdP_[IPHelix::phi0_] = 2 * sin(lhelix.phi0()) / (lhelix.omega() * cos(lhelix.phi0()) * cos(lhelix.phi0()) * cos(lhelix.phi0()));
      dDdP_[IPHelix::omega_] = 1 - tan(lhelix.phi0()) * tan(lhelix.phi0())/lhelix.omega();
      dDdP_[IPHelix::z0_] = 0;
      dDdP_[IPHelix::tanDip_] = 0;

      // no spatial dependence, DT is purely temporal
      dTdP_[IPHelix::t0_] = 1.0; // time is 100% correlated
      // propagate T0 errors to the error on doca
      ddoca_ = ROOT::Math::Similarity(dDdP(),lhelix.params().covariance());
      // dot product between directions at POCA
      Vec3 pdir, sdir;
      lhelix.direction(particleToca(),pdir);
      tline.direction(sensorToca(),sdir);
      ddot_ = pdir.Dot(sdir);
    }
  }

  // specialization between a piecewise IPHelix and a line
  typedef PKTraj<IPHelix> PLHELIX;
  template<> TPoca<PLHELIX,TLine>::TPoca(PLHELIX const& phelix, TLine const& tline, double precision) : TPocaBase(precision), ktraj_(&phelix), straj_(&tline)  {
    // iteratively find the nearest piece, and POCA for that piece.  Start in the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter=0;
    size_t oldindex= phelix.pieces().size();
    size_t index = size_t(rint(oldindex/2.0));
    status_ = converged; 
    while(status_ == converged && niter++ < maxiter && index != oldindex){
      // call down to IPHelix TPoca
      // prepare for the next iteration
      IPHelix const& piece = phelix.pieces()[index];
      TPoca<IPHelix,TLine> tpoca(piece,tline,precision);
      status_ = tpoca.status();
      if(tpoca.usable()){
        // copy over the rest of the state
        partPoca_ = tpoca.particlePoca();
        sensPoca_ = tpoca.sensorPoca();
        doca_ = tpoca.doca();
        dDdP_ = tpoca.dDdP();
        dTdP_ = tpoca.dTdP();
        ddoca_ = tpoca.dDoca();
        ddot_ = tpoca.dirDot();
      }
      oldindex = index;
      index = phelix.nearestIndex(tpoca.particlePoca().T());
    }
    if(status_ == converged && niter >= maxiter) status_ = unconverged;
  }

}
