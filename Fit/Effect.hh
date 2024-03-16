#ifndef KinKal_Effect_hh
#define KinKal_Effect_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/Chisq.hh"
#include "KinKal/Fit/FitState.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/General/TimeRange.hh"
#include <array>
#include <memory>
#include <ostream>

namespace KinKal {

  template<class KTRAJ> class Effect {
    public:
      // type of the data payload used for processing the fit
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      Effect() {}
      virtual ~Effect(){}
      // Effect interface
      virtual double time() const = 0; // time of this effect
      virtual bool active() const = 0; // whether this effect is/was used in the fit
      // Add this effect to the ongoing fit in the given direction.
      virtual void process(FitState& kkdata,TimeDir tdir) = 0;
      // update this effect for a new algebraic iteration, perhaps with a new config
      virtual void updateState(MetaIterConfig const& miconfig,bool first) = 0;
      // update this effect for a new configuration
      virtual void updateConfig(Config const& config) =0;
      // add this effect to a trajectory in the given direction
      virtual void append(PTRAJ& fit,TimeDir tdir) =0;
      // update the reference trajectory for this effect
      virtual void updateReference(PTRAJ const& ptraj) =0;
      // chisquared WRT a given local parameter set, assumed uncorrelatedd  This is used for convergence testing
      virtual Chisq chisq(Parameters const& pdata) const  = 0;
      // diagnostic printout
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      // disallow copy and equivalence
      Effect(Effect const& ) = delete;
      Effect& operator =(Effect const& ) = delete;
  };

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Effect<KTRAJ> const& eff) {
    ost << (eff.active() ? "Active " : "Inactive ") << "time " << eff.time();
    return ost;
  }
}

#endif
