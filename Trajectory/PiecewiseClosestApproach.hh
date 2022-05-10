#ifndef KinKal_PiecewiseClosestApproach_hh
#define KinKal_PiecewiseClosestApproach_hh
//
// ClosestApproach class dealing with piecewise trajectories
//
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <ostream>

namespace KinKal {
  template<class KTRAJ, class STRAJ> class PiecewiseClosestApproach : public ClosestApproach<ParticleTrajectory<KTRAJ>,STRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using KTCA = ClosestApproach<KTRAJ,STRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      PiecewiseClosestApproach(PKTRAJ const& pktraj, STRAJ const& straj, CAHint const& hint, double precision);
      // provide access to the local (non-piecewise) information implicit in this class
      size_t particleTrajIndex() const { return pindex_; }
      KTRAJ const& localParticleTraj() const { return this->particleTraj().piece(pindex_); }
      KTRAJPTR const& localTraj() const { return this->particleTraj().indexTraj(pindex_); }
      KTCA localClosestApproach() const { return KTCA(localTraj(),this->sensorTraj(),this->precision(),this->tpData(),this->dDdP(),this->dTdP()); }
    private:
      size_t pindex_; // indices to the local traj used in TCA calculation
  };

  template<class KTRAJ, class STRAJ> PiecewiseClosestApproach<KTRAJ,STRAJ>::PiecewiseClosestApproach(ParticleTrajectory<KTRAJ> const& pktraj, STRAJ const& straj, CAHint const& hint, double prec) : ClosestApproach<ParticleTrajectory<KTRAJ>,STRAJ>(pktraj,straj,prec) {
    // iteratively find the nearest piece, and CA for that piece.  Start at hints if availalble, otherwise the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter TODO
    unsigned niter=0;
    size_t oldindex= this->particleTraj().pieces().size();
    this->pindex_ = this->particleTraj().nearestIndex(hint.particleToca_);
    // copy over the hint: it needs to evolve
    CAHint phint = hint;
    // iterate until TCA is on the same piece
    do{
      KTCA tpoca(this->particleTraj().piece(pindex_),this->sensorTraj(),phint,prec);
      // copy the state
      this->tpdata_ = tpoca.tpData();
      this->dDdP_ = tpoca.dDdP();
      this->dTdP_ = tpoca.dTdP();
      //      inrange = tpoca.inRange();
      // update the hint
      phint.particleToca_ = tpoca.particleToca();
      phint.sensorToca_ = tpoca.sensorToca();
      // update the piece (if needed)
      oldindex = pindex_;
      pindex_ = this->particleTraj().nearestIndex(tpoca.particlePoca().T());
    } while( pindex_ != oldindex && this->usable() && niter++ < maxiter);
    // overwrite the status if we oscillated on the piece
    if(this->tpdata_.status() == ClosestApproachData::converged && niter >= maxiter)
      this->tpdata_.status_ = ClosestApproachData::unconverged;
    // should test explicitly for piece oscillation FIXME!
    // test if the solution is on a cusp and if so, chose the one with the smallest DOCA TODO
  }
}
#endif

