#include "KinKal/TPoca.hh"
#include "KinKal/IPHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/PKTraj.hh"
#include <limits>
// specializations for TPoca
using namespace std;
namespace KinKal {
  // specialization between a looping helix and a line
  template<> TPoca<IPHelix,TLine>::TPoca(IPHelix const& iphelix, TLine const& tline, TPocaHint const& hint, double precision) : TPocaBase(precision),ktraj_(&iphelix), straj_(&tline) {
    // reset status
    reset();
    double htoca,stoca;
    // initialize the helix time using hints, if available.  If not, use the Z of the line
    if(hint.particleHint_)
      htoca = hint.particleToca_;
    else
      htoca = iphelix.ztime(tline.z0());
    // similar for line; this shouldn't matter, since the solution is linear
    if(hint.sensorHint_)
      stoca = hint.sensorToca_;
    else
      stoca = tline.t0();
    // use successive linear approximation until desired precision on DOCA is met.
    double dptoca(std::numeric_limits<double>::max()), dstoca(std::numeric_limits<double>::max());
    double doca(0.0);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter(0);
    // helix speed doesn't change
    double hspeed = iphelix.speed(iphelix.t0());
    // std::cout << "HSpeed " << hspeed <<  " " << iphelix.omega() << std::endl;
    // iterate until change in TOCA is less than precision
    while((fabs(dptoca) > precision_ || fabs(dstoca) > precision_) && niter++ < maxiter) {
      // find helix local position and direction
      Vec3 hpos = iphelix.position(htoca);
      Vec3 hdir = iphelix.direction(htoca);
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
      dptoca = (hdd - ldd*ddot)/(denom*hspeed);
      dstoca = tline.t0() + (hdd*ddot - ldd)/(denom*tline.speed(stoca)) - stoca;
      htoca += dptoca; // helix time is iterative
      stoca += dstoca; // line time is always WRT t0, since it uses p0
      // compute DOCA
      hpos = iphelix.position(htoca);
      Vec3 lpos = tline.position(stoca);
      // std::cout << "HPos " << hpos.Z() << " " << htoca << std::endl;
      // std::cout << "LPos " << lpos.X() <<  " " << lpos.Y() << " " << lpos.Z() << std::endl;
      double dd2 = (hpos-lpos).Mag2();
      if(dd2 < 0.0 ){
        status_ = pocafailed;
        break;
      }
      doca = sqrt(dd2);
      // std::cout << "DOCA " << doca << std::endl;
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
      // set the TPOCA 4-vectors
      partPoca_.SetE(htoca);
      sensPoca_.SetE(stoca);
      iphelix.position(partPoca_);
      tline.position(sensPoca_);
      // sign doca by angular momentum projected onto difference vector
      double lsign = tline.dir().Cross(iphelix.direction(partPoca_.T())).Dot(sensPoca_.Vect()-partPoca_.Vect());
      double dsign = copysign(1.0,lsign);
      doca_ = doca*dsign;

      // pre-compute some values needed for the derivative calculations
      double time = particlePoca().T();
      Vec3 ddir = delta().Vect().Unit();// direction vector along D(POCA) from traj 2 to 1 (line to helix)

      double l = iphelix.translen(CLHEP::c_light * iphelix.beta() * (time - iphelix.t0()));
      double phi = iphelix.phi(time);
      double phi0 = iphelix.phi0();
      double w = iphelix.omega();
      double d0 = iphelix.d0();

      // Mom4 mom = iphelix.momentum(time);
      // double pt = sqrt(mom.perp2());
      // double radius = fabs(pt);

      // no t0 dependence, DOCA is purely geometric
      dDdP_[IPHelix::omega_] = -dsign * (1./(w*w) * (-sin(phi) + l*w*cos(phi) + sin(phi0)) * ddir.x() +
                                         1./(w*w) * ( cos(phi) + l*w*sin(phi) - cos(phi0)) * ddir.y());
      dDdP_[IPHelix::d0_] = -dsign * (-sin(phi0) * ddir.x() + cos(phi0) * ddir.y());
      dDdP_[IPHelix::phi0_] = -dsign * (1./w * (cos(phi) - (d0*w+1)*cos(phi0)) * ddir.x() +
                                        1./w * (sin(phi) - (d0*w+1)*sin(phi0)) * ddir.y());
      dDdP_[IPHelix::tanDip_] = -dsign * l * ddir.z();
      dDdP_[IPHelix::z0_] = -dsign * ddir.z();

      // no spatial dependence, DT is purely temporal
      dTdP_[IPHelix::t0_] = -1.0; // time is 100% correlated
      // propagate parameter covariance to variance on doca and toca
      docavar_ = ROOT::Math::Similarity(dDdP(),iphelix.params().covariance());
      tocavar_ = ROOT::Math::Similarity(dTdP(),iphelix.params().covariance());
      // dot product between directions at POCA
      ddot_ = iphelix.direction(particleToca()).Dot(tline.direction(sensorToca()));
    }
  }

  // specialization between a piecewise IPHelix and a line
  typedef PKTraj<IPHelix> PIPHelix;
  template<> TPoca<PIPHelix,TLine>::TPoca(PIPHelix const& phelix, TLine const& tline, TPocaHint const& hint, double precision) : TPocaBase(precision), ktraj_(&phelix), straj_(&tline)  {
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
      // call down to IPHelix TPoca
      // prepare for the next iteration
      IPHelix const& piece = phelix.pieces()[index];
      TPoca<IPHelix,TLine> tpoca(piece,tline,hint,precision);
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
