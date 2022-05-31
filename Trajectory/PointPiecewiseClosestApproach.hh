#ifndef KinKal_PointPiecewiseClosestApproach_hh
#define KinKal_PointPiecewiseClosestApproach_hh
//
// ClosestApproach class dealing with piecewise trajectories
//
#include "KinKal/Trajectory/PointClosestApproach.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class PointPiecewiseClosestApproach : public PointClosestApproach<ParticleTrajectory<KTRAJ>> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTCA = PointClosestApproach<KTRAJ>;
      // the constructor is the only non-inherited function
      PointPiecewiseClosestApproach(PTRAJ const& pktraj, VEC4 const& point, PCAHint const& hint, double precision);
      PointPiecewiseClosestApproach(PTRAJ const& pktraj, VEC4 const& point, double precision);
         // provide access to the local (non-piecewise) information implicit in this class
      size_t particleTrajIndex() const { return pindex_; }
      KTRAJ const& localParticleTraj() const { return this->particleTraj().piece(pindex_); }
      KTCA localClosestApproach() const { return KTCA(localParticleTraj(),this->point(),this->precision(),this->tpData(),this->dDdP(),this->dTdP()); }
 private:
      size_t pindex_; // indices to the local traj used in TCA calculation
  };

  template<class KTRAJ> PointPiecewiseClosestApproach<KTRAJ>::PointPiecewiseClosestApproach(PTRAJ const& pktraj, VEC4 const& point, double prec) :
  PointPiecewiseClosestApproach(pktraj,point, PCAHint(point.T()), prec) {}

  // iteratively find the nearest piece, and CA for that piece.  Start at hints if availalble, otherwise the middle
  template<class KTRAJ> PointPiecewiseClosestApproach<KTRAJ>::PointPiecewiseClosestApproach(PTRAJ const& pktraj, VEC4 const& point, PCAHint const& hint, double prec) :
    KTCA(pktraj,point,prec) {
    // iteratively find the nearest piece, and CA for that piece.  Start at hints if availalble, otherwise the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter=0;
    size_t oldindex= this->particleTraj().pieces().size();
    pindex_ = this->particleTraj().nearestIndex(hint.particleToca_);
    // copy over the hint: it needs to evolve
    PCAHint phint = hint;
    // iterate until TCA is on the same piece
    do{
      KTCA tpoca(this->particleTraj().piece(pindex_),this->point(),phint,this->precision());
      // copy the state
      this->tpdata_ = tpoca.tpData();
      this->dDdP_ = tpoca.dDdP();
      this->dTdP_ = tpoca.dTdP();
      //      inrange = tpoca.inRange();
      // update the hint
      phint.particleToca_ = tpoca.particleToca();
      // update the piece (if needed)
      oldindex = pindex_;
      pindex_ = this->particleTraj().nearestIndex(tpoca.particlePoca().T());
    } while( pindex_ != oldindex && usable() && niter++ < maxiter);
    // overwrite the status if we oscillated on the piece
    if(this->tpdata_.status() == ClosestApproachData::converged && niter >= maxiter)
      this->tpdata_.status_ = ClosestApproachData::unconverged;
    // should test explicitly for piece oscillation FIXME!
    // test if the solution is on a cusp and if so, chose the one with the smallest DOCA TODO

  }
}
#endif

