#ifndef KinKal_ParticleState_hh
#define KinKal_ParticleState_hh
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

  class ParticleState {
    public:
      static size_t constexpr dimension() { return 6; }
      static std::string const& stateName(size_t index);
      static std::string const& stateUnit(size_t index);
      static std::string const& stateTitle(size_t index);

      // construct from position and momentum vectors
      ParticleState(VEC3 const& pos, VEC3 const& mom, double time, double mass) : state_(pos.X(),pos.Y(),pos.Z(),mom.X(),mom.Y(),mom.Z()), time_(time), mass_(mass) {}
      ParticleState(VEC4 const& pos4, MOM4 const& mom4) : state_(pos4.Vect().X(),pos4.Vect().Y(),pos4.Vect().Z(),mom4.Vect().X(),mom4.Vect().Y(),mom4.Vect().Z()), time_(pos4.E()), mass_(mom4.M()) {}
      // construct from raw information
      ParticleState(SVEC const& state, double time, double mass) : state_(state), time_(time), mass_(mass) {}
      ParticleState() {}
      // direct accessor
      SVEC const& state() const { return state_; }
      // explicit component accessors.  Note these return by value.  Unfortunately Root doesn't provide a more elegant conversion operator
      VEC3 position3() const { return VEC3(state_[0],state_[1],state_[2]); }
      VEC4 position4() const { return VEC4(state_[0],state_[1],state_[2],time_); }
      VEC3 momentum3() const { return VEC3(state_[3],state_[4],state_[5]); }
      MOM4 momentum4() const { return MOM4(state_[3],state_[4],state_[5],mass_); }

      double mass() const { return mass_; } 
      double time() const { return time_; } 
    private:
      SVEC state_; // state vector payload of this particle, with position and momentum information
      double time_; // time particle had this state
      double mass_; // particle mass
      static std::vector<std::string> stateTitles_;
      static std::vector<std::string> stateNames_;
      static std::vector<std::string> stateUnits_;
  };

  // the same, including covariance ( when used as a measurement with estimated errors)
  class ParticleStateMeasurement  {
    public:
    // construct from from raw information
      ParticleStateMeasurement(ParticleState const& state, SMat const& scovar) : state_(state), scovar_(scovar) {}
      ParticleStateMeasurement(SVEC const& state, SMat const& scovar,double time, double mass) : state_(state, time, mass), scovar_(scovar) {}
      ParticleStateMeasurement() {}
      ParticleState const& stateVector() const { return state_; }
      SMat const& stateCovariance() const { return scovar_; }
    private:
      ParticleState state_; // state
      SMat scovar_; // covariance of state vector
  };
}
#endif
