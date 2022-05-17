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
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      Effect() {}
      virtual ~Effect(){}
      // Effect interface
      virtual double time() const = 0; // time of this effect
      virtual bool active() const = 0; // whether this effect is/was used in the fit
      // Add this effect to the ongoing fit in the given direction.
      virtual void process(FitState& kkdata,TimeDir tdir) = 0;
      // update this effect for a new reference trajectory within the existing algebraic iteration sequence
      virtual void update(PKTRAJ const& ref) = 0;
      // update this effect to start a new algebraic iteration squence using the new reference trajectory and configuration
      virtual void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) = 0;
      // update this effect for a new configuration
      virtual void update(Config const& config) =0;
      // update the particle trajectory for this effect
      virtual void append(PKTRAJ& fit) =0;
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
