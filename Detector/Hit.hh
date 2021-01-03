#ifndef KinKal_Hit_hh
#define KinKal_Hit_hh
//
//  Base class to describe a measurement that constrains some parameters of the fit
//  The hit constraint is described as a set of Residuals.
//  Residuals may constrain any physical aspect of the fit (time, position, time+position, momentum, ...)
//  Each residual esimates the 1-dimenaional difference between the measurement and the value predicted by a reference trajectory
//  The hit may be associated with a piece of detector material as well
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/General/Weights.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Detector/Residual.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Fit/Config.hh"
#include <memory>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class Hit {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
     // default
      Hit(){}
      virtual ~Hit(){}
      // disallow copy and equivalence
      Hit(Hit const& ) = delete; 
      Hit& operator =(Hit const& ) = delete;
      // access the cache used by every hit
      KTRAJ const& refTraj() const { return reftraj_; }
      // the constraint this hit implies WRT the current reference, expressed as a weight
      Weights const& weight() const { return weight_; }
      // number of degrees of freedom constrained by this hit; this is just the number of active residuals
      virtual unsigned nDOF() const;
      // convert active residuals into chisq, WRT either reference or other parameters
      double chisq() const; // chisq WRT reference parameters
      double chisq(Parameters const& params) const;
      // hits are active if any residual is active
      bool active() const { return nDOF() > 0; }
      // residuals corrected to refer to the given set of parameters (1st-order)
      Residual residual(Parameters const& params, unsigned ires) const;
      // virtual interface, provided by subclasses
      // number of residuals in this hit
      virtual unsigned nResid() const = 0;
      // individual residuals may be active or inactive
      virtual bool active(unsigned ires) const = 0;
      // reference residuals for this hit.  iDOF indexs the measurement and is hit-specific, outside the range will throw
      virtual Residual const& residual(unsigned ires) const = 0;
      // time of this hit: this is WRT the reference trajectory
      virtual double time() const = 0;
      // update to a new reference, without changing state
      virtual void update(PKTRAJ const& pktraj) = 0;
      // update the internals of the hit, specific to this meta-iteraion
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& config) = 0;
      // associated material information; null means no material
      virtual EXINGPTR const& detXingPtr() const = 0;
      bool hasMaterial() const { return (bool)detXingPtr(); }
      virtual void print(std::ostream& ost=std::cout,int detail=0) const = 0;
    protected:
      KTRAJ reftraj_; // reference trajectory, used to compute reference residuals
      void setWeight(); // set the weight based on current active residuals
    private:
      Weights weight_; // measurement weight WRT most recent parameters
  };

  template <class KTRAJ> unsigned HIT<ktraj>::nDOF() const {
    unsigned retval(0);
    for(unsigned ires=0; ires< nResid(); ires++) 
      if(active(ires)) retval++;
    return retval;
  }

  template <class KTRAJ> double HIT<ktraj>::chisq() const {
    double retval(0.0);
    for(unsigned ires=0; ires< nResid(); ires++) {
      if(active(ires)) {
	// find the reference residual
	auto const& res = residual(ires);
	// project the parameter covariance into a residual space variance
	double rvar = ROOT::Math::Similarity(res.dRdP(),pdata.covariance());
	// add the measurement variance
	rvar +=  res.variance();
	// add chisq for this DOF
	retval += (res.value()*res.value())/rvar;
      }
    }
    return retval;
  }

  template <class KTRAJ> double HIT<ktraj>::chisq(Parameters const& params) const {
    double retval(0.0);
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = pdata.parameters() - reftraj_.params().parameters(); 
    for(unsigned ires=0; ires< nResid(); ires++) {
      if(active(ires)) {
	// compute the residual WRT the given parameters
	auto res = residual(params,ires);
	// use the differnce to 'correct' the reference residual value to be WRT these parameters
	double uresid = res.value() - ROOT::Math::Dot(dpvec,res.dRdP());
	// project the parameter covariance into a residual space variance
	double rvar = ROOT::Math::Similarity(res.dRdP(),pdata.covariance());
	// add the measurement variance
	rvar +=  res.variance();
	// add chisq for this DOF
	retval += (uresid*uresid)/rvar;
      }
    }
    return retval;
  }

  template <class KTRAJ> Residual Hit<KTRAJ>::residual(Parameters const& pdata,unsigned ires) const {
    auto const& resid = residual(ires);
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = pdata.parameters() - reftraj_.params().parameters(); 
    // use the differnce to 'correct' the reference residual to be WRT these parameters
    double uresid = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
    return Residual(uresid,resid.variance(),resid.dRdP());
  }

  template <class KTRAJ> void HIT<ktraj>::setWeight() {
    // start with a null weight
    weight_ = Weight();
    for(unsigned ires=0; ires< nResid(); ires++) {
      if(active(ires)) {
 	auto const& res = residual(ires);
	// convert derivatives vector to a Nx1 matrix
	ROOT::Math::SMatrix<double,NParams(),1> dRdPM;
	dRdPM.Place_in_col(res.dRdP(),0,0);
	// convert the variance into a 1X1 matrix
	ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> RVarM;
	// weight by inverse variance
	double tvar = res.variance();
	RVarM(0,0) = 1.0/tvar;
	// expand these into the weight matrix
	DMAT wmat = ROOT::Math::Similarity(dRdPM,RVarM);
	// translate residual value into weight vector WRT the reference parameters
	// sign convention reflects resid = measurement - prediction
	DVEC wvec = wmat*reftraj_.params().parameters() + res.dRdP()*res.value()/tvar;
	// weights are linearly additive
	weight_ += Weights(wvec,wmat);
      }
    }
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, Hit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
 }

}
#endif

