#ifndef KinKal_Effect_hh
#define KinKal_Effect_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Fit/FitState.hh"
#include "KinKal/Fit/Chisq.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/TimeRange.hh"
#include <array>
#include <memory>
#include <ostream>

namespace KinKal {
 
  template<class KTRAJ> class Effect {
    public:
      enum State{unprocessed=-1,processed,updated,failed};
      static std::string const& stateName(State state);
      // type of the data payload used for processing the fit
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      // properties common to all effects
      virtual double time() const = 0; // time of this effect
      virtual bool active() const = 0; // whether this effect is/was used in the fit
      // Add this effect to the ongoing fit in the given direction.
      virtual void process(FitState& kkdata,TimeDir tdir) = 0;
      // update this effect for a new reference trajectory within the existing algebraic iteration sequence
      virtual void update(PKTRAJ const& ref) = 0;
      // update this effect to start a new algebraic iteration squence using the new reference trajectory and configuration
      virtual void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) = 0;
      // diagnostic printout
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      // the following only has a non-trivial implementation for effects which (potentially) add information content to the fit
      // chisquared (quality) associated with the most recent processing and current reference.  This is used to determine
      // fit convergence
      virtual Chisq chisq() const { return Chisq();} 
      // chisquared WRT a given local parameter set.  This is a purely diagnostic function
      virtual Chisq chisq(Parameters const& pdata) const { return Chisq();} // chisq contribution WRT parameters 
      // The following only has a non-trivial implemetation for effects which (potentially) alter the physical particle trajectory
      virtual void append(PKTRAJ& fit) {};
      // disallow copy and equivalence
      Effect(Effect const& ) = delete; 
      Effect& operator =(Effect const& ) = delete; 
      State state(TimeDir tdir) const { return state_[static_cast<std::underlying_type<TimeDir>::type>(tdir)]; }
      void setState(TimeDir tdir, State state) { state_[static_cast<std::underlying_type<TimeDir>::type>(tdir)] = state; }
      bool wasProcessed(TimeDir tdir) const { return state(tdir) == processed; }
      void updateState() { state_[0] = state_[1] = updated; }
      Effect() : state_{{unprocessed,unprocessed}} {}
      virtual ~Effect(){} 
    private:
      std::array<State,2> state_; // state of processing in each direction
  };
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Effect<KTRAJ> const& eff) {
    ost << (eff.active() ? "Active " : "Inactive ") << "time " << eff.time() << " state " <<
    TimeDir::forwards << " " << eff.stateName(eff.state(TimeDir::forwards))  << " : " <<
    TimeDir::backwards << " " << eff.stateName(eff.state(TimeDir::backwards));
    return ost;
  }

  template <class KTRAJ> std::string const& Effect<KTRAJ>::stateName(Effect::State state) {
    const static std::vector<std::string> stateNames_ = { "Unprocessed", "Processed", "Updated", "Failed" };
    switch (state) {
      case unprocessed: default:
	return stateNames_[0];
      case processed:
	return stateNames_[1];
      case updated:
	return stateNames_[2];
      case failed:
	return stateNames_[3];
    }
  }
}

#endif
