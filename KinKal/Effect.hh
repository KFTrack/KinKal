#ifndef KinKal_Effect_hh
#define KinKal_Effect_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include "KinKal/ParticleTrajectory.hh"
#include "KinKal/FitState.hh"
#include "KinKal/Config.hh"
#include <array>
#include <memory>
#include <ostream>

namespace KinKal {
 
  template<class KTRAJ> class Effect {
    public:
      enum Status{unprocessed=-1,processed,updated,failed};
      static std::string const& statusName(Status status);
      // type of the data payload used for processing the fit
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      virtual double time() const = 0; // time of this effect
      virtual unsigned nDOF() const {return 0; }; // how/if this effect contributes to the measurement NDOF
      virtual bool isActive() const = 0; // whether this effect is/was used in the fit
       // Add this effect to the ongoing fit in a give direction.
      virtual void process(FitState& kkdata,TimeDir tdir) = 0;
      virtual double fitChi() const { return 0.0;} // unbiased chi contribution of this effect after fitting
      virtual double chisq(Parameters const& pdata) const { return 0.0;} // chisq contribution WRT parameters 
      // update this effect for a new reference trajectory; this is used in the algebraic iteration step
      virtual void update(PKTRAJ const& ref) = 0;
      // update this effect for a new configuration and reference trajectory; this is used to start a meta-iteration step
      virtual void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) = 0;
      // append this effects trajectory change (if appropriate)
      virtual void append(PKTRAJ& fit) {};
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      virtual ~Effect(){} 
      // disallow copy and equivalence
      Effect(Effect const& ) = delete; 
      Effect& operator =(Effect const& ) = delete; 
      Status status(TimeDir tdir) const { return status_[static_cast<std::underlying_type<TimeDir>::type>(tdir)]; }
      void setStatus(TimeDir tdir, Status status) { status_[static_cast<std::underlying_type<TimeDir>::type>(tdir)] = status; }
      bool wasProcessed(TimeDir tdir) const { return status(tdir) == processed; }
      void updateStatus() { status_[0] = status_[1] = updated; }
      Effect() : status_{{unprocessed,unprocessed}} {}
    private:
      std::array<Status,2> status_; // status of processing in each direction
  };
  
  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Effect<KTRAJ> const& eff) {
    ost << (eff.isActive() ? "Active " : "Inactive ") << "time " << eff.time() << " status " <<
    TimeDir::forwards << " " << eff.statusName(eff.status(TimeDir::forwards))  << " : " <<
    TimeDir::backwards << " " << eff.statusName(eff.status(TimeDir::backwards));
    return ost;
  }

  template <class KTRAJ> std::string const& Effect<KTRAJ>::statusName(Effect::Status status) {
    static std::vector<std::string> statusNames_ = { "Unprocessed", "Processed", "Updated", "Failed" };
    switch (status) {
      case unprocessed: default:
	return statusNames_[0];
      case processed:
	return statusNames_[1];
      case updated:
	return statusNames_[2];
      case failed:
	return statusNames_[3];
    }
  }
}

#endif
