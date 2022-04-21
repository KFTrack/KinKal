#ifndef KinKal_ClosestApproach_hh
#define KinKal_ClosestApproach_hh
///
//  This functor class finds the (spacetime) points of closest approach between a particle and sensor trajectory
//  Both trajectories must satisfy the 'TTraj' interface
//  Concrete instances are specializations and must be implemented explicity for each trajectory pair
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <iostream>
#include <ostream>

namespace KinKal {
  // Hint class for TCA calculation. TCA search will start at these TOCA values.  This allows to
  // disambiguate cases with multiple solutions (like looping trajectories), or to speed up calculations when an
  // approximate answer is already known.
  struct CAHint{
    double particleToca_, sensorToca_; // approximate values, used as starting points for cacluations
    CAHint(double ptoca,double stoca) :  particleToca_(ptoca), sensorToca_(stoca) {}
  };
  // Class to calculate DOCA and TOCA using time parameterized trajectories.
  // Templated on the types of trajectories. The actual implementations must be specializations for particular trajectory classes.
  template<class KTRAJ, class STRAJ> class ClosestApproach {
    public:
      // construct from the particle and sensor trajectories; TCA is computed on construction, given a hint as to where
      // to start looking, which disambiguates functions with multiple solutions
      ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, CAHint const& hint, double precision);
      // construct without a hint: TCA isn't calculated, state is invalid
      ClosestApproach(KTRAJ const& ptraj, STRAJ const& straj, double precision);
      // explicitly construct from all content (no calculation)
      ClosestApproach(KTRAJ const& ptraj, STRAJ const& straj, double precision,
          ClosestApproachData const& tpdata, DVEC const& dDdP, DVEC const& dTdP);
      // accessors
      ClosestApproachData const& tpData() const { return tpdata_; }
      KTRAJ const& particleTraj() const { return ktraj_; }
      STRAJ const& sensorTraj() const { return straj_; }
      // derviatives of TOCA and DOCA WRT particle trajectory parameters
      DVEC const& dDdP() const { return dDdP_; }
      DVEC const& dTdP() const { return dTdP_; }
      bool inRange() const { return particleTraj().inRange(particleToca()) && sensorTraj().inRange(sensorToca()); }
      double precision() const { return precision_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
      // forward the data payload interface.
      ClosestApproachData::TPStat status() const { return tpdata_.status_; }
      std::string const& statusName() const { return tpdata_.statusName(tpdata_.status_); }
      double doca() const { return tpdata_.doca(); }
      double docaVar() const { return tpdata_.docaVar(); }
      double tocaVar() const { return tpdata_.tocaVar(); }
      double dirDot() const { return tpdata_.dirDot(); }
      double deltaT() const { return tpdata_.deltaT(); }
      bool usable() const { return tpdata_.usable(); }
      double particleToca() const { return tpdata_.particleToca(); }
      double sensorToca() const { return tpdata_.sensorToca(); }
      double lSign() const { return tpdata_.lsign_; } // sign of angular momentum
      VEC4 const& particlePoca() const { return tpdata_.particlePoca(); }
      VEC4 const& sensorPoca() const { return tpdata_.sensorPoca(); }
      VEC4 delta() const { return tpdata_.delta(); }
      VEC3 const& particleDirection() const { return tpdata_.particleDirection(); }
      VEC3 const& sensorDirection() const { return tpdata_.sensorDirection(); }
      // calculate CA given the hint, and fill the state
      void findTCA(CAHint const& hint);
    private:
      double precision_; // precision used to define convergence
      KTRAJ const& ktraj_; // kinematic particle trajectory
      STRAJ const& straj_; // sensor trajectory
    protected:
      ClosestApproachData tpdata_; // data payload of CA calculation
      DVEC dDdP_; // derivative of DOCA WRT Parameters
      DVEC dTdP_; // derivative of TOCA WRT Parameters
  };

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, double prec) :
    precision_(prec),ktraj_(ktraj), straj_(straj) {}

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, double prec,
    ClosestApproachData const& tpdata, DVEC const& dDdP, DVEC const& dTdP) :
   precision_(prec),ktraj_(ktraj), straj_(straj), tpdata_(tpdata),dDdP_(dDdP), dTdP_(dTdP) {}

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, CAHint const& hint,
      double prec) : ClosestApproach(ktraj,straj,prec) {
    findTCA(hint);
  }

  template<class KTRAJ, class STRAJ> void ClosestApproach<KTRAJ,STRAJ>::findTCA(CAHint const& hint) {
    // reset status
    tpdata_.reset();
    // initialize TOCA using hints
    tpdata_.partCA_.SetE(hint.particleToca_);
    tpdata_.sensCA_.SetE(hint.sensorToca_);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter(0);
    // speed doesn't change
    double pspeed = ktraj_.speed(particleToca());
    double sspeed = straj_.speed(sensorToca());
    // iterate until change in TOCA is less than precision
    double dptoca(std::numeric_limits<double>::max()), dstoca(std::numeric_limits<double>::max());
    while(tpdata_.usable() && (fabs(dptoca) > precision() || fabs(dstoca) > precision()) && niter++ < maxiter) {
      // find positions and directions at the current TOCA estimate
      tpdata_.partCA_ = ktraj_.position4(tpdata_.particleToca());
      tpdata_.sensCA_ = straj_.position4(tpdata_.sensorToca());
      tpdata_.pdir_ = ktraj_.direction(particleToca());
      tpdata_.sdir_ = straj_.direction(sensorToca());
      VEC3 dpos = sensorPoca().Vect()-particlePoca().Vect();
      // dot products
      double ddot = sensorDirection().Dot(particleDirection());
      double denom = 1.0 - ddot*ddot;
      // check for parallel)
      if(denom<1.0e-5){
        tpdata_.status_ = ClosestApproachData::pocafailed;
        break;
      }
      double hdd = dpos.Dot(particleDirection());
      double ldd = dpos.Dot(sensorDirection());
      // compute the change in times
      dptoca = (hdd - ldd*ddot)/(denom*pspeed);
      dstoca = (hdd*ddot - ldd)/(denom*sspeed);
      // update the TOCA estimates
      tpdata_.partCA_.SetE(particleToca()+dptoca);
      tpdata_.sensCA_.SetE(sensorToca()+dstoca);
    }
    if(tpdata_.status_ != ClosestApproachData::pocafailed){
      if(niter < maxiter)
        tpdata_.status_ = ClosestApproachData::converged;
      else
        tpdata_.status_ = ClosestApproachData::unconverged;
      // need to add divergence and oscillation tests FIXME!
    }
    // final update
    tpdata_.partCA_ = ktraj_.position4(tpdata_.particleToca());
    tpdata_.sensCA_ = straj_.position4(tpdata_.sensorToca());
    tpdata_.pdir_ = ktraj_.direction(particleToca());
    tpdata_.sdir_ = straj_.direction(sensorToca());
    // fill the rest of the state
    if(usable()){
      // sign doca by angular momentum projected onto difference vector
      VEC3 dvec = delta().Vect();
      tpdata_.lsign_ = copysign(1.0,sensorDirection().Cross(particleDirection()).Dot(dvec));
      tpdata_.doca_ = dvec.R()*tpdata_.lsign_;
      VEC3 dvechat = dvec.Unit();
      // now variances due to the particle trajectory parameter covariance
      // for DOCA, project the spatial position derivative along the delta-CA direction
      DVDP dxdp = ktraj_.dXdPar(particleToca());
      SVEC3 dv(dvechat.X(),dvechat.Y(),dvechat.Z());
      dDdP_ = -dv*dxdp;
      dTdP_[KTRAJ::t0Index()] = -1.0;  // TOCA is 100% anti-correlated with the (mandatory) t0 component.
      // project the parameter covariance onto DOCA and TOCA
      tpdata_.docavar_ = ROOT::Math::Similarity(dDdP(),ktraj_.params().covariance());
      tpdata_.tocavar_ = ROOT::Math::Similarity(dTdP(),ktraj_.params().covariance());
    }
  }

  template<class KTRAJ, class STRAJ> void ClosestApproach<KTRAJ,STRAJ>::print(std::ostream& ost,int detail) const {
    ost << "ClosestApproach status " << statusName() << " Doca " << doca() << " +- " << sqrt(docaVar())
      << " dToca " << deltaT() << " +- " << sqrt(tocaVar()) << " cos(theta) " << dirDot() << " Precision " << precision() << std::endl;
    if(detail > 0)
      ost << "Particle Poca " << particlePoca() << " Sensor Poca " << sensorPoca() << std::endl;
    if(detail > 1)
      ost << "dDdP " << dDdP() << " dTdP " << dTdP() << std::endl;
    if(detail > 2){
      ost << "Particle ";
      particleTraj().print(ost,0);
      ost << "Sensor ";
      sensorTraj().print(ost,0);
    }
  }

}
#endif
