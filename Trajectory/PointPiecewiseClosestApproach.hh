#ifndef KinKal_PointPiecewiseClosestApproach_hh
#define KinKal_PointPiecewiseClosestApproach_hh
//
// ClosestApproach class dealing with piecewise trajectories
//
#include "KinKal/Trajectory/PointClosestApproach.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class PointPiecewiseClosestApproach : public PointClosestApproach<KTRAJ> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTCA = PointClosestApproach<KTRAJ>;
      // the constructor is the only non-inherited function
      PointPiecewiseClosestApproach(PTRAJ const& ptraj, VEC4 const& point, PCAHint const& hint, double precision);
      PointPiecewiseClosestApproach(PTRAJ const& ptraj, VEC4 const& point, double precision);
  };

  template<class KTRAJ> PointPiecewiseClosestApproach<KTRAJ>::PointPiecewiseClosestApproach(PTRAJ const& ptraj, VEC4 const& point, double prec) :
  PointPiecewiseClosestApproach(ptraj,point, PCAHint(point.T()), prec) {}

  // iteratively find the nearest piece, and CA for that piece.  Start at hints if availalble, otherwise the middle
  template<class KTRAJ> PointPiecewiseClosestApproach<KTRAJ>::PointPiecewiseClosestApproach(PTRAJ const& ptraj, VEC4 const& point, PCAHint const& hint, double prec) :
    KTCA(ptraj.nearestTraj(hint.particleToca_),point,prec) {
    // iteratively find the nearest piece, and CA for that piece
    static const unsigned maxiter=10; // don't allow infinite iteration.
    unsigned niter=0;
    size_t oldindex= ptraj.pieces().size();
    auto pindex = ptraj.nearestIndex(hint.particleToca_);
    // copy over the hint: it needs to evolve
    PCAHint phint = hint;
    // iterate until TCA is on the same piece
    do{
      this->findTCA(phint);
      phint.particleToca_ = tpoca.particleToca();
      oldindex = pindex;
      pindex = ptraj.nearestIndex(tpoca.particlePoca().T());
      this->ktrajptr_ = ptraj.indexTraj(pindex);
    } while( pindex != oldindex && usable() && niter++ < maxiter);
    // overwrite the status if we oscillated on the piece
    if(this->tpdata_.status() == ClosestApproachData::converged && pindex != oldindex)
      this->tpdata_.status_ = ClosestApproachData::unconverged;
    // should test explicitly for piece oscillation FIXME!
    // test if the solution is on a cusp and if so, chose the one with the smallest DOCA TODO

  }
}
#endif

