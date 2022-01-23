#ifndef KinKal_PointPiecewiseClosestApproach_hh
#define KinKal_PointPiecewiseClosestApproach_hh
//
// ClosestApproach class dealing with piecewise trajectories
//
#include "KinKal/Trajectory/PointClosestApproach.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class PointPiecewiseClosestApproach {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using KTCA = PointClosestApproach<KTRAJ>;
      // the constructor is the only non-inherited function
      PointPiecewiseClosestApproach(PKTRAJ const& pktraj, VEC4 const& point, PCAHint const& hint, double precision);
      PointPiecewiseClosestApproach(PKTRAJ const& pktraj, VEC4 const& point, double precision);
      // copy the TCA interface.  This is ugly and a maintenance burden, but avoids inheritance problems
      ClosestApproachData::TPStat status() const { return tpdata_.status(); }
      std::string const& statusName() const { return tpdata_.statusName(); }
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
      ClosestApproachData const& tpData() const { return tpdata_; }
      PKTRAJ const& particleTraj() const { return pktraj_; }
      size_t particleTrajIndex() const { return pindex_; }
      DVEC const& dDdP() const { return dDdP_; }
      DVEC const& dTdP() const { return dTdP_; }
      bool inRange() const { return particleTraj().inRange(particleToca()); }
      double precision() const { return precision_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
    private:
      double precision_; // precision used to define convergence
      ClosestApproachData tpdata_; // data payload of CA calculation
      PKTRAJ const& pktraj_;
      size_t pindex_; // indices to the local traj used in TCA calculation
      DVEC dDdP_;
      DVEC dTdP_;
  };

  template<class KTRAJ> PointPiecewiseClosestApproach<KTRAJ>::PointPiecewiseClosestApproach(PKTRAJ const& pktraj, VEC4 const& point, double prec) :
  PointPiecewiseClosestApproach(pktraj,point, PCAHint(point.T()), prec) {}

  // iteratively find the nearest piece, and CA for that piece.  Start at hints if availalble, otherwise the middle
  template<class KTRAJ> PointPiecewiseClosestApproach<KTRAJ>::PointPiecewiseClosestApproach(PKTRAJ const& pktraj, VEC4 const& point, PCAHint const& hint, double prec) : precision_(prec), pktraj_(pktraj){
    // iteratively find the nearest piece, and CA for that piece.  Start at hints if availalble, otherwise the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter=0;
    size_t oldindex= pktraj_.pieces().size();
    pindex_ = pktraj_.nearestIndex(hint.particleToca_);
    // copy over the hint: it needs to evolve
    PCAHint phint = hint;
    // iterate until TCA is on the same piece
    do{
      KTCA tpoca(pktraj_.piece(pindex_),point,phint,prec);
      // copy the state
      tpdata_ = tpoca.tpData();
      dDdP_ = tpoca.dDdP();
      dTdP_ = tpoca.dTdP();
      //      inrange = tpoca.inRange();
      // update the hint
      phint.particleToca_ = tpoca.particleToca();
      // update the piece (if needed)
      oldindex = pindex_;
      pindex_ = pktraj_.nearestIndex(tpoca.particlePoca().T());
    } while( pindex_ != oldindex && usable() && niter++ < maxiter);
    // overwrite the status if we oscillated on the piece
    if(tpdata_.status() == ClosestApproachData::converged && niter >= maxiter)
      tpdata_.status_ = ClosestApproachData::unconverged;
    // should test explicitly for piece oscillation FIXME!
    // test if the solution is on a cusp and if so, chose the one with the smallest DOCA TODO
  }


  template<class KTRAJ> void PointPiecewiseClosestApproach<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "PointPiecewiseClosestApproach status " << statusName() << " Doca " << doca() << " +- " << sqrt(docaVar())
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

