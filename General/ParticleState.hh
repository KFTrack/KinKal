#ifndef KinKal_ParticleState_hh
#define KinKal_ParticleState_hh
// Describe the state vector (position and momentum) of a particle.  This (with mass and charge) unambgiuously defines the particle kinematics
// so this class can provide a full kinematic interface.  The algebraic parameterization can be used as both input and output to a fit.
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <vector>
#include <string>

namespace KinKal {

  class ParticleState {
    public:
      static size_t constexpr dimension() { return 6; }
      static std::string const& stateName(size_t index);
      static std::string const& stateUnit(size_t index);
      static std::string const& stateTitle(size_t index);

      // construct from position and momentum vectors
      ParticleState(VEC3 const& pos, VEC3 const& mom, double time, double mass, int charge) : state_(pos.X(),pos.Y(),pos.Z(),mom.X(),mom.Y(),mom.Z()), time_(time), mass_(mass), charge_(charge) {}
      ParticleState(VEC4 const& pos4, MOM4 const& mom4, int charge) : state_(pos4.Vect().X(),pos4.Vect().Y(),pos4.Vect().Z(),mom4.Vect().X(),mom4.Vect().Y(),mom4.Vect().Z()), time_(pos4.E()), mass_(mom4.M()) , charge_(charge){}
      // construct from raw information
      ParticleState(SVEC6 const& state, double time, double mass, int charge) : state_(state), time_(time), mass_(mass), charge_(charge) {}
      ParticleState() : time_(0.0), mass_(0.0), charge_(0) {}
      // direct accessor to the state content in its algebraic form.
      SVEC6 const& state() const { return state_; }
      // explicit component accessors.  Note these return by value.  Unfortunately Root doesn't provide a more elegant conversion operator
      VEC3 position3() const { return VEC3(state_[0],state_[1],state_[2]); }
      VEC4 position4() const { return VEC4(state_[0],state_[1],state_[2],time_); }
      VEC3 momentum3() const { return VEC3(state_[3],state_[4],state_[5]); }
      MOM4 momentum4() const { return MOM4(state_[3],state_[4],state_[5],mass_); }
      double momentum() const { return momentum3().R(); }
      double mass() const { return mass_; }
      double time() const { return time_; }
      int charge() const { return charge_; }
      double energy() const { auto mom = momentum(); return sqrt(mom*mom + mass_*mass_); }
      double beta() const { auto mom = momentum(); return mom/sqrt(mom*mom + mass_*mass_); }
      double gamma() const { return energy()/mass_; }
      double speed() const { return CLHEP::c_light*beta(); }
      VEC3 velocity() const { return speed()*momentum3()/momentum(); }
    private:
      SVEC6 state_; // state vector payload of this particle, with position and momentum information
      double time_; // time particle had this state
      double mass_; // particle mass
      int charge_; // charge in units of the proton charge
      const static std::vector<std::string> stateTitles_;
      const static std::vector<std::string> stateNames_;
      const static std::vector<std::string> stateUnits_;
  };
}
#endif
