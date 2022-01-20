#ifndef KinKal_PointClosestApproach_hh
#define KinKal_PointClosestApproach_hh
///
//  This functor class finds the (spacetime) point of closest approach between a particle and
//  a point in space-time
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include <iostream>
#include <ostream>

namespace KinKal {
  // Hint class for TCA calculation. TCA search will start at this TOCA values.  This allows to
  // disambiguate cases with multiple solutions (like looping trajectories), or to speed up calculations when an
  // approximate answer is already known.
  struct PCAHint{
    double particleToca_; // approximate values, used as starting points for cacluations
    PCAHint(double ptoca) :  particleToca_(ptoca) {}
  };
  template<class KTRAJ> class PointClosestApproach {
    public:
      // construct from the particle and sensor trajectories; TCA is computed on construction, given a hint as to where
      // to start looking, which disambiguates functions with multiple solutions
      PointClosestApproach(KTRAJ const& ktraj, VEC4 const& point, PCAHint const& hint, double precision);
      // construct without a hint: TCA isn't calculated, state is invalid
      PointClosestApproach(KTRAJ const& ptraj, VEC4 const& point, double precision);
      // accessors
      ClosestApproachData const& tpData() const { return tpdata_; }
      KTRAJ const& particleTraj() const { return ktraj_; }
      // derviatives of TOCA and DOCA WRT particle trajectory parameters
      DVEC const& dDdP() const { return dDdP_; }
      DVEC const& dTdP() const { return dTdP_; }
      bool inRange() const { return particleTraj().inRange(particleToca()); }
      double precision() const { return precision_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
      // calculate CA given the hint, and fill the state
      void findTCA(PCAHint const& hint);
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
      VEC4 const& point() const { return tpdata_.sensorPoca(); }
      VEC4 delta() const { return tpdata_.delta(); }
      VEC3 const& particleDirection() const { return tpdata_.particleDirection(); }
      VEC3 const& pointDirection() const { return tpdata_.sensorDirection(); }
      // calculate CA given the hint, and fill the state
  private:
      double precision_; // precision used to define convergence
      ClosestApproachData tpdata_; // data payload of CA calculation
      KTRAJ const& ktraj_; // kinematic particle trajectory
      DVEC dDdP_; // derivative of DOCA WRT Parameters
      DVEC dTdP_; // derivative of TOCA WRT Parameters
  };

  template<class KTRAJ> PointClosestApproach<KTRAJ>::PointClosestApproach(KTRAJ const& ktraj, VEC4 const& point, double prec) :
      PointClosestApproach(ktraj,point,PCAHint(point.T()), prec) {}

  template<class KTRAJ> PointClosestApproach<KTRAJ>::PointClosestApproach(KTRAJ const& ktraj, VEC4 const& point, PCAHint const& hint,
      double prec) : precision_(prec), ktraj_(ktraj) {
    // sensor CA is fixed to the point
    tpdata_.sensCA_ = point;
    findTCA(hint);
  }

   template<class KTRAJ> void PointClosestApproach<KTRAJ>::findTCA(PCAHint const& hint) {
    // reset status
    tpdata_.reset();
    // initialize TOCA using hints
    tpdata_.partCA_.SetE(hint.particleToca_);
    static const unsigned maxiter=100; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter(0);
    // speed doesn't change
    double pspeed = ktraj_.speed(particleToca());
    // iterate until change in TOCA is less than precision
    double dptoca(std::numeric_limits<double>::max());
    while(tpdata_.usable() && fabs(dptoca) > precision() && niter++ < maxiter) {
      // find positions and directions at the current TOCA estimate
      tpdata_.partCA_ = ktraj_.position4(tpdata_.particleToca());
      tpdata_.pdir_ = ktraj_.direction(particleToca());
      auto dpos = point().Vect()-particlePoca().Vect();
      // compute the change in times
      dptoca = dpos.Dot(tpdata_.pdir_)/pspeed;
      // update the TOCA estimates
      tpdata_.partCA_.SetE(particleToca()+dptoca);
    }
    if(niter < maxiter)
      tpdata_.status_ = ClosestApproachData::converged;
    else
      tpdata_.status_ = ClosestApproachData::unconverged;
    // final update
    tpdata_.partCA_ = ktraj_.position4(tpdata_.particleToca());
    tpdata_.pdir_ = ktraj_.direction(particleToca());
    tpdata_.sdir_ = tpdata_.delta().Vect();
    // fill the rest of the state
    if(usable()){
      // sign doca by angular momentum projected onto difference vector
      VEC3 dvec = delta().Vect();
      tpdata_.lsign_ = copysign(1.0,pointDirection().Cross(particleDirection()).Dot(dvec));
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


  template<class KTRAJ> void PointClosestApproach<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "PointClosestApproach status " << statusName() << " Doca " << doca() << " +- " << sqrt(docaVar())
      << " dToca " << deltaT() << " +- " << sqrt(tocaVar()) << " cos(theta) " << dirDot() << " Precision " << precision() << std::endl;
    if(detail > 0)
      ost << "Particle Poca " << particlePoca() << " Point " << point() << std::endl;
    if(detail > 1)
      ost << "dDdP " << dDdP() << " dTdP " << dTdP() << std::endl;
    if(detail > 2){
      ost << "Particle ";
      particleTraj().print(ost,0);
    }
  }

}
#endif
