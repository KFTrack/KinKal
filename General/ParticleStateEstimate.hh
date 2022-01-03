// particle state with estimated errors (covariance), used to seed a fit,
// or to express a fit result as a particle state with errors.  Transforming between ParticleStateEstimate
// and a kinematic trajectory is lossless.
#ifndef KinKal_ParticleStateEstimate_hh
#define KinKal_ParticleStateEstimate_hh
#include "KinKal/General/ParticleState.hh"

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
    private:
      DMAT scovar_; // covariance of state vector
  };
}
#endif
