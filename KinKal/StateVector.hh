#ifndef KinKal_StateVector_hh
#define KinKal_StateVector_hh
// Describe the state vector (position and momentum) of a particle.  This (with mass and charge) unambgiuously defines the particle and
// can be used as input or output of a track fit
#include "KinKal/Vectors.hh"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
namespace KinKal {
  typedef ROOT::Math::SVector<double,6> SVEC; // type for state vector payload
  typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> > SMat;  // matrix type for state vector covariance
  typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepStd<double,6,6> > DSDP; // matrix for state derivatives WRT parameters and vice-versa
  typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepStd<double,6,6> > DPDS; // matrix for state derivatives WRT parameters and vice-versa

  class StateVector {
    public:
      // construct from position and momentum 3-vectors
      StateVector(Vec3 const& pos, Vec3 const& mom) : state_(pos.X(),pos.Y(),pos.Z(),mom.X(),mom.Y(),mom.Z()) {}
      // construct from raw information
      StateVector(SVEC const& state) : state_(state) {}
      StateVector() {}
      // direct accessor
      SVEC const& state() const { return state_; }
      // explicit component accessors.  Note these return by value.  Unfortunately Root doesn't provide a more elegant conversion operator
      Vec3 position() const { return Vec3(state_[0],state_[1],state_[2]); }
      Vec3 momentum() const { return Vec3(state_[3],state_[4],state_[5]); }
    private:
      SVEC state_; // state vector payload of this particle
  };

  // the same, including covariance ( when used as a measurement with estimated errors)
  class StateVectorMeasurement  {
    public:
    // construct from from raw information
      StateVectorMeasurement(StateVector const& state, SMat const& scovar) : state_(state.state()), scovar_(scovar) {}
      StateVectorMeasurement(SVEC const& state, SMat const& scovar) : state_(state), scovar_(scovar) {}
      StateVectorMeasurement() {}
      StateVector const& stateVector() const { return state_; }
      SMat const& stateCovariance() const { return scovar_; }
    private:
      StateVector state_; // state
      SMat scovar_; // covariance of state vector
  };
}
#endif
