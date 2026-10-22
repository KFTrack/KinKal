#ifndef KinKal_PiecewiseClosestApproach_hh
#define KinKal_PiecewiseClosestApproach_hh
//
// ClosestApproach class dealing with piecewise trajectories
//
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <ostream>

namespace KinKal {
  template<class KTRAJ, class STRAJ> class PiecewiseClosestApproach : public ClosestApproach<KTRAJ,STRAJ> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTCA = ClosestApproach<KTRAJ,STRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using STRAJPTR = std::shared_ptr<STRAJ>;
      // maxiter is the maximum number of moves to the piece holding the current solution. With maxiter = 0 the hint
      // time defines the reference piece, not just a starting value: the solution is computed on the piece holding the
      // hint time and its TOCA may lie outside that piece's range.
      PiecewiseClosestApproach(PTRAJ const& ptraj, STRAJPTR strajptr, CAHint const& hint, double precision,
          unsigned maxiter=10);
  };

  template<class KTRAJ, class STRAJ> PiecewiseClosestApproach<KTRAJ,STRAJ>::PiecewiseClosestApproach(ParticleTrajectory<KTRAJ> const& ptraj, STRAJPTR strajptr, CAHint const& hint, double prec,
      unsigned maxiter) : KTCA(ptraj.nearestTraj(hint.particleToca_),strajptr,prec) {
    // find the CA on the piece holding the hint time
    auto pindex = ptraj.nearestIndex(hint.particleToca_);
    this->findTCA(hint);
    if(maxiter == 0) return;
    // iterate until the local TCA is on the same piece, moving at most maxiter times
    unsigned niter=0;
    while(true){
      auto oldindex = pindex;
      pindex = ptraj.nearestIndex(this->particleToca());
      this->ktrajptr_ = ptraj.indexTraj(pindex);
      if(pindex == oldindex || !this->usable()) break;
      if(niter++ >= maxiter){
        // overwrite the status if we didn't converge on the piece
        if(this->status() == ClosestApproachData::converged) this->tpdata_.status_ = ClosestApproachData::unconverged;
        break;
      }
      this->findTCA(this->hint());
    }
    // should test explicitly for piece oscillation FIXME!
    // test if the solution is on a cusp and if so, chose the one with the smallest DOCA TODO
  }
}
#endif

