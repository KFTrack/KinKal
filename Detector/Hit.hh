#ifndef KinKal_Hit_hh
#define KinKal_Hit_hh
//
//  Base class to describe a measurement that constrains some aspect of the track fit
//  The hit may be associated with a piece of detector material as well
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/General/Weights.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Chisq.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Fit/Config.hh"
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class Hit {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      // default
      Hit(){}
      virtual ~Hit(){}
      // disallow copy and equivalence
      Hit(Hit const& ) = delete; 
      Hit& operator =(Hit const& ) = delete;
      // the constraint this hit implies WRT the current reference, expressed as a weight
      virtual Weights weight() const =0;
      // hits may be active (used in the fit) or inactive; this is a pattern recognition feature
      virtual bool active() const =0;
      virtual Chisq chisq() const =0; // least-squares distance to reference parameters
      virtual Chisq chisq(Parameters const& params) const =0;  // least-squares distance to given parameters
      virtual double time() const = 0;  // time of this hit: this is WRT the reference trajectory
      // update to a new reference, without changing state
      virtual void update(PKTRAJ const& pktraj) = 0;
      // update the internals of the hit, specific to this meta-iteraion
      virtual void updateState(PKTRAJ const& pktraj, MetaIterConfig const& config) = 0;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const = 0;
  };

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Hit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
  }

}
#endif

