#ifndef KinKal_ClosestApproach_hh
#define KinKal_ClosestApproach_hh
///
//  This functor class finds the (spacetime) points of closest approach between a particle and sensor trajectory
//  Both trajectories must satisfy the 'TTraj' interface
//  Concrete instances are specializations and must be implemented explicity for each trajectory pair
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <memory>
#include <iostream>
#include <ostream>
#include <cmath>

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
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      // construct from the particle and sensor trajectories; TCA is computed on construction, given a hint as to where
      // to start looking, which disambiguates functions with multiple solutions
      ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, CAHint const& hint, double precision);
      // same, using Ptrs
      ClosestApproach(KTRAJPTR const& ktrajptr, STRAJ const& straj, CAHint const& hint, double precision);
      // construct without a hint: TCA isn't calculated, state is invalid
      ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, double precision);
      ClosestApproach(KTRAJPTR const& ktrajptr, STRAJ const& straj, double precision);
      // explicitly construct from all content (no calculation)
      ClosestApproach(KTRAJPTR const& ktrajptr, STRAJ const& straj, double precision,
          ClosestApproachData const& tpdata, DVEC const& dDdP, DVEC const& dTdP);
      // accessors
      ClosestApproachData const& tpData() const { return tpdata_; }
      KTRAJ const& particleTraj() const { return *ktrajptr_; }
      KTRAJPTR const& particleTrajPtr() const { return ktrajptr_; }
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
      // return the hint from the current state
      CAHint hint() const { return CAHint(particleToca(),sensorToca()); }
      // equivalence
      ClosestApproach& operator = (ClosestApproach const& other);
      // other accessors
      void setTrajectory(KTRAJPTR ktrajptr){ ktrajptr_ = ktrajptr; }
    private:
      double precision_; // precision used to define convergence
      KTRAJPTR ktrajptr_; // kinematic particle trajectory
      STRAJ const& straj_; // sensor trajectory This should be a shared ptr TODO
    protected:
      // calculate CA given the hint, and fill the state
      void findTCA(CAHint const& hint);
      ClosestApproachData tpdata_; // data payload of CA calculation
      DVEC dDdP_; // derivative of DOCA WRT Parameters
      DVEC dTdP_; // derivative of TOCA WRT Parameters
      void setParticleTOCA( double ptoca);
      void setSensorTOCA( double stoca);
  };

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, double prec) :
    precision_(prec),ktrajptr_(new KTRAJ(ktraj)), straj_(straj) {}

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJPTR const& ktrajptr, STRAJ const& straj, double prec) :
    precision_(prec),ktrajptr_(ktrajptr), straj_(straj) {}

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJPTR const& ktrajptr, STRAJ const& straj, double prec,
    ClosestApproachData const& tpdata, DVEC const& dDdP, DVEC const& dTdP) :
   precision_(prec),ktrajptr_(ktrajptr), straj_(straj), tpdata_(tpdata),dDdP_(dDdP), dTdP_(dTdP) {}

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJ const& ktraj, STRAJ const& straj, CAHint const& hint,
      double prec) : ClosestApproach(ktraj,straj,prec) {
    findTCA(hint);
  }

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>::ClosestApproach(KTRAJPTR const& ktrajptr, STRAJ const& straj, CAHint const& hint,
      double prec) : ClosestApproach(ktrajptr,straj,prec) {
    findTCA(hint);
  }

  template<class KTRAJ, class STRAJ> ClosestApproach<KTRAJ,STRAJ>& ClosestApproach<KTRAJ,STRAJ>::operator = (ClosestApproach const& other)
  {
    tpdata_ = other. tpData();
    dDdP_ = other.dDdP();
    dTdP_ = other.dTdP();
    ktrajptr_ = other.ktrajptr_;
    // make sure the sensor traj is the same
    if(&straj_ != &other.sensorTraj()) throw std::invalid_argument("Inconsistent ClosestApproach SensorTraj");
    return *this;
  }

  template<class KTRAJ, class STRAJ> void ClosestApproach<KTRAJ,STRAJ>::findTCA(CAHint const& hint) {
    // reset status
    tpdata_.reset();
    // initialize using hints
    setParticleTOCA(hint.particleToca_);
    setSensorTOCA(hint.sensorToca_);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter TODO
    unsigned niter(0);
    // speed doesn't change
    double pspeed = ktrajptr_->speed(particleToca());
    double sspeed = straj_.speed(sensorToca());
    // iterate until change in TOCA is less than precision
    double dptoca(std::numeric_limits<double>::max()), dstoca(std::numeric_limits<double>::max());
    while(tpdata_.usable() && (fabs(dptoca) > precision() || fabs(dstoca) > precision()) && niter++ < maxiter) {
      // find positions and directions at the current TOCA estimate
      VEC3 delta = sensorPoca().Vect()-particlePoca().Vect();
      double ddot = sensorDirection().Dot(particleDirection());
      double denom = 1.0 - ddot*ddot;
      // check for error conditions
      if(denom<1.0e-5){
        tpdata_.status_ = ClosestApproachData::parallel;
        break;
      }
      if(std::isnan(denom)){
        tpdata_.status_ = ClosestApproachData::failed;
        break;
      }
      double pdd = delta.Dot(particleDirection());
      double sdd = delta.Dot(sensorDirection());
      // compute the change in TOCA
      dptoca = (pdd - sdd*ddot)/(denom*pspeed);
      dstoca = (pdd*ddot - sdd)/(denom*sspeed);
      // update the TOCA estimates
      setParticleTOCA(particleToca()+dptoca);
      setSensorTOCA(sensorToca()+dstoca);
      // need to add divergence and oscillation tests TODO
    }
    // fill the rest of the state
    if(usable()){
      if(fabs(dptoca) < precision() && fabs(dstoca) < precision() )tpdata_.status_ = ClosestApproachData::converged;
     // sign doca by angular momentum projected onto difference vector.  This is just a convention
      VEC3 dvec = delta().Vect();
      VEC3 pdir = particleDirection();
      tpdata_.lsign_ = copysign(1.0,sensorDirection().Cross(pdir).Dot(dvec));
      tpdata_.doca_ = dvec.R()*tpdata_.lsign_;
      // now variances due to the particle trajectory parameter covariance
      // for DOCA, project the spatial position derivative along the delta-CA direction
      DVDP dxdp = ktrajptr_->dXdPar(particleToca()); // position change WRT parameters
      // change in DOCA WRT parameters; this is straight forwards
      VEC3 dvechat = dvec.Unit();
      SVEC3 dvh(dvechat.X(),dvechat.Y(),dvechat.Z());
      dDdP_ = tpdata_.lsign_*dvh*dxdp;
      // change in particle DeltaT WRT parameters; both PTOCA and STOCA terms are important.
      // The dependendence on the momentum direction change is an order of magnitude smaller but
      // might be included someday TODO

      double dd = sensorDirection().Dot(particleDirection());
      double denom = (1.0-dd*dd);
      double pfactor = 1.0/(pspeed*denom);
      VEC3 beta = pfactor*(sensorDirection()*dd - particleDirection());
      SVEC3 sbeta(beta.X(),beta.Y(),beta.Z());
      double sfactor = 1.0/(sspeed*denom);
      VEC3 gamma = sfactor*(sensorDirection() - dd*particleDirection());
      SVEC3 sgamma(gamma.X(),gamma.Y(),gamma.Z());
      dTdP_ = (sbeta - sgamma)*dxdp;
      tpdata_.docavar_ = ROOT::Math::Similarity(dDdP(),ktrajptr_->params().covariance());
      tpdata_.tocavar_ = ROOT::Math::Similarity(dTdP(),ktrajptr_->params().covariance());
    }
  }

  template<class KTRAJ, class STRAJ> void ClosestApproach<KTRAJ,STRAJ>::setParticleTOCA(double ptoca) {
    tpdata_.partCA_ = ktrajptr_->position4(ptoca);
    tpdata_.pdir_ = ktrajptr_->direction(ptoca);
  }

  template<class KTRAJ, class STRAJ> void ClosestApproach<KTRAJ,STRAJ>::setSensorTOCA(double stoca) {
    tpdata_.sensCA_ = straj_.position4(stoca);
    tpdata_.sdir_ = straj_.direction(stoca);
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
