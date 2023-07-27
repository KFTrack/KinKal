// particle state with estimated errors (covariance), used to seed a fit,
// or to express a fit result as a particle state with errors.  Transforming between ParticleStateEstimate
// and a kinematic trajectory is lossless.
#ifndef KinKal_ParticleStateEstimate_hh
#define KinKal_ParticleStateEstimate_hh
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/MomBasis.hh"

namespace KinKal {
  class ParticleStateEstimate : public ParticleState {
    public:
      // construct from from raw information
      ParticleStateEstimate(ParticleState const& state, DMAT const& scovar) : ParticleState(state), scovar_(scovar) {}
      ParticleStateEstimate(SVEC6 const& state, DMAT const& scovar,double time, double mass, int charge) : ParticleState(state, time, mass, charge), scovar_(scovar) {}
      ParticleStateEstimate() {}
      DMAT const& stateCovariance() const { return scovar_; }
      // project the variance onto the scalar momentum
      double momentumVariance() const;
      // project the position variance in given direction.  Note this will throw if given the momentum direction, as that variance is infinite
      double positionVariance(MomBasis::Direction dir) const;
      // project the position variance onto the Plane defined by the perp direction (u) and the phi direction (v)
      PMAT planeCovariance() const;
    private:
      DMAT scovar_; // covariance of state vector
  };
}
#endif
