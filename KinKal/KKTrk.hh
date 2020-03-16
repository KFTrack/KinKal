#ifndef KinKal_KKTrk_hh
#define KinKal_KKTrk_hh
//
// Primary class of the Kinematic Kalman fit.  This class owns the state describing
// the fit (hits, material interactions, BField corrections) and coordinates the
// fit processing.  It uses
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the kinematic trajectory used in this fit
//
#include "KinKal/PKTraj.hh"
#include "KinKal/KKEff.hh"
#include "KinKal/KKEnd.hh"
#include "KinKal/KKMHit.hh"
#include "KinKal/KKHit.hh"
#include "KinKal/KKMat.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/THit.hh"
#include "KinKal/FitStatus.hh"
//#include "KinKal/DetElem.hh"
#include <set>
#include <vector>
#include <memory>
#include <stdexcept>

namespace KinKal {
  // struct to define fit parameters
  struct Config {
    unsigned maxniter_; // maximum number of iterations
    double dwt_; // dweighting of initial seed covariance
    double mindchisq_; // minimum change in chisquared for convergence
    bool addmat_; // add material
    double tbuff_; // time buffer for final fit
    Config() : maxniter_(10), dwt_(1.0e6), mindchisq_(1.0e-3), addmat_(true), tbuff_(10.0) {} 
  };

  template<class KTRAJ> class KKTrk {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      struct KKEFFComp { // comparator, used in sorting
	bool operator()(std::unique_ptr<KKEFF> const& a, std::unique_ptr<KKEFF> const&  b) const {
	  if(a.get() != b.get())
	    return a->time() < b->time();
	  else
	    return false;
	}
      };
      typedef std::set<std::unique_ptr<KKEFF>,KKEFFComp > KKCON; // container type for effects
      // construct from a set of hits (and materials FIXME!)
      KKTrk(PKTRAJ const& reftraj, std::vector<const THit*> hits, Config const& config ); 
      void fit(bool force=false); // process the effects.
      FitStatus const& status() const { return status_; }
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKCON const& effects() const { return effects_; }
    private:
      // helper functions
      bool update();
      bool fitIteration();
      // payload
      Config config_; // configuration
      FitStatus status_; // current fit status
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated and the current fit performed
      PKTRAJ fittraj_; // result of the current fit, becomes the reference when the fit is iterated
      KKCON effects_; // effects used in this fit, sorted by time
  };

  template <class KTRAJ> KKTrk<KTRAJ>::KKTrk(PKTRAJ const& reftraj, std::vector<const THit*> hits,Config const& config) : 
    config_(config), reftraj_(reftraj), fittraj_(reftraj.range(),reftraj.mass(),reftraj.charge()){
      // loop over the hits
      for(auto hit : hits ) {
	// create the hit effects and insert them in the set
	// if there's associated material, create a combined material and hit effect, otherwise just a hit effect
	if(config_.addmat_ && hit->material() != 0){
	  auto iemp = effects_.emplace(std::make_unique<KKMHit<KTRAJ> >(*hit,reftraj));
	  if(!iemp.second)throw std::runtime_error("Insertion failure");
	} else{ 
	  auto iemp = effects_.emplace(std::make_unique<KKHit<KTRAJ> >(*hit,reftraj));
	  if(!iemp.second)throw std::runtime_error("Insertion failure");
	}
      }
      //add pure material and BField inhomogeneity effects FIXME!
      // reset the range if necessary
      if( reftraj_.range().infinite()
	  || reftraj_.range().low() > effects_.begin()->get()->time()
	  || reftraj_.range().high() < effects_.rbegin()->get()->time()){
	reftraj_.setRange(TRange(std::min(reftraj_.range().low(),effects_.begin()->get()->time() - config_.tbuff_),
	    std::max(reftraj_.range().high(),effects_.rbegin()->get()->time() + config_.tbuff_)));
      }
      // create end effects
      auto iemp = effects_.emplace(std::make_unique<KKEnd<KTRAJ>>(reftraj,TDir::forwards,config_.dwt_));
      if(!iemp.second)throw std::runtime_error("Insertion failure");
      iemp = effects_.emplace(std::make_unique<KKEnd<KTRAJ>>(reftraj,TDir::backwards,config_.dwt_));
      if(!iemp.second)throw std::runtime_error("Insertion failure");
    }

  // fit iteration management.  This should cover simulated annealing too FIXME!
  template <class KTRAJ> void KKTrk<KTRAJ>::fit(bool force) {
    // don't fit unless necessary or forced
    if(status_.status_ == FitStatus::needsfit || force ){
      // iterate to convergence
      status_.niter_ = 0;
      bool fitOK(true);
      while(fitOK && status_.niter_ < config_.maxniter_ && status_.status_ != FitStatus::converged) {
	if(status_.niter_>0) fitOK = update();
	double oldchisq = status_.chisq_;
	if(fitOK) fitOK = fitIteration(); // this call updates chisq
	// test for convergence/divergence 
	if(fabs(oldchisq-status_.chisq_) < config_.mindchisq_) status_.status_ = FitStatus::converged;
	status_.niter_++;
      }
      if(!fitOK){
	status_.status_ = FitStatus::failed;
      } else if ( status_.status_ != FitStatus::converged) 
	status_.status_ = FitStatus::unconverged;
    } else {
      status_.status_ = FitStatus::converged;
    }
  }

  // single fit iteration 
  template <class KTRAJ> bool KKTrk<KTRAJ>::fitIteration() {
    bool retval(true);
    status_.ndof_ = -(int)KTRAJ::NParams();
    status_.chisq_ = 0.0;
    // fit in both directions
    for(TDir tdir=TDir::forwards; tdir != TDir::end; ++tdir){
    // start with empty fit information; each effect will modify this as necessary, and cache what it needs for later processing
      KKData<KTRAJ::NParams()> fitdata;
      if(tdir==TDir::forwards){
	auto ieff = effects_.begin();
	while(ieff != effects_.end()){
	  // update chisq, NDOF only forwards to save time
	  if(ieff->get()->nDOF() > 0){
	    status_.ndof_ += ieff->get()->nDOF();
	    // I'm not sure this is the most efficient way to get chisq FIXME!
	    status_.chisq_ += ieff->get()->chisq(fitdata.pData());
	  }
	  if(!ieff->get()->process(fitdata,tdir)){
	    retval = false;
	    break;
	  }
	  ieff++;
	}
      } else {
	auto ieff = effects_.rbegin();
	while(ieff != effects_.rend()){
	  if(!ieff->get()->process(fitdata,tdir)){
	    retval = false;
	    break;
	  }
	  ieff++;
	}
      }
    }
    // convert the fit result into a new trajectory
    fittraj_ = PKTRAJ(reftraj_.range(),reftraj_.mass(),reftraj_.charge());
    // process forwards by convention
    for(auto& ieff : effects_) {
      ieff->append(fittraj_);
    }
    return retval;
  }

  // update fit internal state between iterations
  template <class KTRAJ> bool KKTrk<KTRAJ>::update() {
    //swap the fit trajectory to the reference
    reftraj_ = fittraj_;
    // update the effects to use the new reference
    bool retval = true;
    for(auto& ieff : effects_) {
      retval &= ieff->update(reftraj_);
      if(!retval)break;
    }
    return retval;
  }
}
#endif
