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
      PiecewiseClosestApproach(PTRAJ const& ptraj, STRAJPTR strajptr, CAHint const& hint, double precision);
  };

  template<class KTRAJ, class STRAJ> PiecewiseClosestApproach<KTRAJ,STRAJ>::PiecewiseClosestApproach(ParticleTrajectory<KTRAJ> const& ptraj, STRAJPTR strajptr, CAHint const& hint, double prec) : KTCA(ptraj.nearestTraj(hint.particleToca_),strajptr,prec) {
    // iteratively find the nearest piece, and CA for that piece.  Start at hint
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter TODO
    unsigned niter=0;
    size_t oldindex= ptraj.pieces().size();
    auto pindex = ptraj.nearestIndex(hint.particleToca_);
    // copy over the hint: it needs to evolve
    CAHint phint = hint;
    // iterate until the local TCA is on the same piece
    do{
      this->findTCA(phint);
      phint = this->hint();
      oldindex = pindex;
      pindex = ptraj.nearestIndex(hint.particleToca_);
      this->ktrajptr_ = ptraj.indexTraj(pindex);
    } while( pindex != oldindex && this->usable() && niter++ < maxiter);
    // overwrite the status if we didn't converge on the piece
    if(this->status() == ClosestApproachData::converged && pindex != oldindex)
      this->tpdata_.status_ = ClosestApproachData::unconverged;
    // should test explicitly for piece oscillation FIXME!
    // test if the solution is on a cusp and if so, chose the one with the smallest DOCA TODO
  }
}
#endif

