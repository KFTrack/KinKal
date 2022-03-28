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
#include "KinKal/Fit/MetaIterConfig.hh"
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class Hit {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      // default
      Hit() : wscale_(1.0){}
      virtual ~Hit(){}
      // disallow copy and equivalence
      Hit(Hit const& ) = delete;
      Hit& operator =(Hit const& ) = delete;
      Weights const& weight() const { return weight_; }
      // hits may be active (used in the fit) or inactive; this is a pattern recognition feature
      virtual bool active() const =0;
      virtual Chisq chisq(Parameters const& params) const =0;  // least-squares distance to given parameters
      virtual Chisq chisq() const =0;  // least-squares distance to reference parameters
      virtual double time() const = 0;  // time of this hit: this is WRT the reference trajectory
      // update to a new reference, without changing state
      virtual void update(PKTRAJ const& pktraj);
      // update the internals of the hit, specific to this meta-iteraion
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& config);
      virtual void print(std::ostream& ost=std::cout,int detail=0) const = 0;
      // accessors
      double weightScale() const { return wscale_; }
      Parameters const& referenceParameters() const { return refparams_; }
      // the constraint this hit implies WRT the current reference, expressed as a weight
      Weights scaledweight() const { auto wt =  weight_; wt *= wscale_; return wt; }
      // unbiased least-squares distance to reference parameters
      Chisq chisquared() const;
    protected:
      Weights weight_; // weight representation of the hits constraint.  Subclasses must set this in update
      double wscale_; // current annealing weight scaling
      Parameters refparams_; // reference parameters, used to compute reference residuals
  };

  template <class KTRAJ> void Hit<KTRAJ>::update(PKTRAJ const& pktraj) {
  // update the reference parameters
    refparams_ = pktraj.nearestPiece(time()).params();
  }

  template <class KTRAJ> void Hit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    wscale_ = 1.0/miconfig.varianceScale();
    update(pktraj);
  }

  template<class KTRAJ> Chisq Hit<KTRAJ>::chisquared() const {
    if(active()){
      // subtract out the effect of this hit's weight from the reference parameters
      Weights wt(refparams_);
      wt -= weight();
      Parameters uparams(wt);
      return chisq(uparams);
    } else
      return Chisq();
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Hit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
  }

}
#endif

