#ifndef KinKal_Hit_hh
#define KinKal_Hit_hh
//
//  Base class to describe a measurement that constrains some aspect of the track fit
//  The constraint is expressed as a weight WRT a set of reference parameters.
//  The base class is a purely algebraic object.
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
      Hit(PKTRAJ const& pktraj,double time) : reftraj_(pktraj.nearestPiece(time)),wscale_(1.0){}
      Hit(KTRAJ const& ktraj) : reftraj_(ktraj),wscale_(1.0){}
      virtual ~Hit(){}
      // disallow copy and equivalence
      Hit(Hit const& ) = delete;
      Hit& operator =(Hit const& ) = delete;
     // hits may be active (used in the fit) or inactive; this is a pattern recognition feature
      virtual bool active() const =0;
      virtual Chisq chisq(Parameters const& params) const =0;  // least-squares distance to given parameters
      virtual double time() const = 0;  // time of this hit: this is WRT the reference trajectory
      // update to a new reference, without changing internal state
      virtual void update(PKTRAJ const& pktraj) = 0;
      // update the internals of the hit, specific to this meta-iteraion
      virtual void update(PKTRAJ const& pktraj,MetaIterConfig const& config) = 0;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const = 0;
      // accessors
      // the constraint this hit implies WRT the current reference, expressed as a weight.  This will be used in the next fit iteration
      Weights const& weight() const { return weight_; }
      // the constraint used in making the current reference
      Weights const& referenceWeight() const { return refweight_; }
      // the same, scaled for annealing
      double weightScale() const { return wscale_; }
      // parameters WRT which this hit's residual and weights are set.  These are generally biased
      // in that they contain the information of this hit
      Parameters const& referenceParameters() const { return reftraj_.params(); }
      // Unbiased parameters, taking out the hits effect
      Parameters unbiasedParameters() const;
      // unbiased least-squares distance to reference parameters
      Chisq chisquared() const;
    protected:
      // update the weight
      void setWeight(Weights const& weight){
        refweight_ = weight_;
        weight_ = weight;
        weight_ *= wscale_;
      }
      KTRAJ reftraj_; // reference parameters used for this hit's weight  Should be private FIXME
      double wscale_; // current annealing weight scaling Should be private FIXME
    private:
      Weights weight_; // weight representation of the hit's constraint
      Weights refweight_; // weight used in the previous iteration
  };

  template<class KTRAJ> Parameters Hit<KTRAJ>::unbiasedParameters() const {
    if(active()){
      // convert the parameters to a weight, and subtract this hit's weight
      Weights wt(referenceParameters());
      // subtract out the effect of this hit's reference weight from the reference parameters
      wt -= referenceWeight();
      return Parameters(wt);
    } else
      return reftraj_.params();
  }

  template<class KTRAJ> Chisq Hit<KTRAJ>::chisquared() const {
    if(active()){
      return chisq(unbiasedParameters());
    } else
      return Chisq();
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Hit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
  }

}
#endif

