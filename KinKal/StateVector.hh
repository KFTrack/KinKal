#ifndef KinKal_StateVector_hh
#define KinKal_StateVector_hh
// Describe the state vector (position and momentum) of a particle.  This (with mass and charge) unambgiuously defines the particle and
// can be used as input or output of a track fit
#include "KinKal/Vectors.hh"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <vector>
#include <string>

namespace KinKal {
  using SVEC = ROOT::Math::SVector<double,6>; // type for state vector payload
  using SMat = ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6>>;  // matrix type for state vector covariance
  using DSDP = ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepStd<double,6,6>>; // matrix for state derivatives WRT parameters and vice-versa
  using DPDS = ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepStd<double,6,6>>; // matrix for state derivatives WRT parameters and vice-versa

  class StateVector {
    public:
      static size_t constexpr dimension() { return 6; }
      static std::string const& stateName(size_t index);
      static std::string const& stateUnit(size_t index);
      static std::string const& stateTitle(size_t index);

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
      static std::vector<std::string> stateTitles_;
      static std::vector<std::string> stateNames_;
      static std::vector<std::string> stateUnits_;
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
