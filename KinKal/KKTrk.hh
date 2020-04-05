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
#include "KinKal/BField.hh"
#include "TMath.h"
#include <set>
#include <vector>
#include <memory>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <ostream>

namespace KinKal {
  // struct to define iteration configuration
  struct Config {
    int maxniter_; // maximum number of iterations for this configurations
    double dwt_; // dweighting of initial seed covariance
    double convdchisq_; // maximum change in chisquared for convergence
    double divdchisq_; // minimum change in chisquared for divergence
    double oscdchisq_; // maximum change in chisquared for oscillation
    bool addmat_; // add material
    bool addfield_; // add field inhomogeneity effects
    double tbuff_; // time buffer for final fit (ns)
    double dtol_; // tolerance on direction change in BField integration (dimensionless)
    double ptol_; // tolerance on position change in BField integration (mm)
    int minndof_; // minimum number of DOFs to continue fit
    Config() : maxniter_(10), dwt_(1.0e6), convdchisq_(0.1), divdchisq_(100.0), oscdchisq_(1.0), addmat_(true), addfield_(true), tbuff_(0.5), dtol_(0.1), ptol_(0.1), minndof_(5) {} 
  };

  template<class KTRAJ> class KKTrk {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef THit<KTRAJ> THIT;
      typedef DXing<KTRAJ> DXING;
      typedef std::shared_ptr<THIT> THITPTR;
      typedef std::shared_ptr<DXING> DXINGPTR;
      typedef std::vector<THITPTR> THITCOL;
      typedef std::vector<DXINGPTR> DXINGCOL;
      struct KKEFFComp { // comparator, used in sorting
	bool operator()(std::unique_ptr<KKEFF> const& a, std::unique_ptr<KKEFF> const&  b) const {
	  if(a.get() != b.get())
	    return a->time() < b->time();
	  else
	    return false;
	}
      };
      typedef std::set<std::unique_ptr<KKEFF>,KKEFFComp > KKEFFCON; // container type for effects
      // construct from a set of hits and passive material crossings
      KKTrk(PKTRAJ const& reftraj, BField const& bfield, THITCOL& thits, DXINGCOL& dxings, Config const& config ); 
      void fit(); // process the effects.  This creates the fit
      void reset() { history_.push_back(FitStatus()); } // reset to allow renewed fitting
      // accessors
      std::vector<FitStatus> const& statusHistory() const { return history_; }
      FitStatus const& fitStatus() const { return history_.back(); } // most recent status
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKEFFCON const& effects() const { return effects_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;

    private:
      // helper functions
      bool update();
      bool fitIteration(FitStatus& status);
      bool canIterate() const;
      bool oscillating(FitStatus const& status) const;
      void createFieldDomains();
      // payload
      Config config_; // configuration
      std::vector<FitStatus> history_; // fit status history
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated and the current fit performed
      PKTRAJ fittraj_; // result of the current fit, becomes the reference when the fit is iterated
      KKEFFCON effects_; // effects used in this fit, sorted by time
      BField const& bfield_; // magnetic field map
  };

  template <class KTRAJ> KKTrk<KTRAJ>::KKTrk(PKTRAJ const& reftraj, BField const& bfield, THITCOL& thits, DXINGCOL& dxings, Config const& config) : 
    config_(config), history_(1,FitStatus()), reftraj_(reftraj), fittraj_(reftraj.range(),reftraj.mass(),reftraj.charge()), bfield_(bfield) {
      // loop over the hits
      for(auto& thit : thits ) {
	// create the hit effects and insert them in the set
	// if there's associated material, create a combined material and hit effect, otherwise just a hit effect
	if(config_.addmat_ && thit->detCrossing() != 0){
	  auto iemp = effects_.emplace(std::make_unique<KKMHit<KTRAJ> >(thit,reftraj));
	  if(!iemp.second)throw std::runtime_error("Insertion failure");
	} else{ 
	  auto iemp = effects_.emplace(std::make_unique<KKHit<KTRAJ> >(thit,reftraj));
	  if(!iemp.second)throw std::runtime_error("Insertion failure");
	}
      }
      //add pure material effects
      for(auto& dxing : dxings) {
	auto iemp = effects_.emplace(std::make_unique<KKMat<KTRAJ> >(dxing,reftraj));
	if(!iemp.second)throw std::runtime_error("Insertion failure");
      }
      // reset the range 
      reftraj_.setRange(TRange(std::min(reftraj_.range().low(),effects_.begin()->get()->time() - config_.tbuff_),
	    std::max(reftraj_.range().high(),effects_.rbegin()->get()->time() + config_.tbuff_)));
      // add BField inhomogeneity effects; not yet implemented FIXME!
      if(config_.addfield_){
	createFieldDomains();
      }
      // create end effects; this should be last to avoid confusing the BField correction
      auto iemp = effects_.emplace(std::make_unique<KKEnd<KTRAJ>>(reftraj,TDir::forwards,config_.dwt_));
      if(!iemp.second)throw std::runtime_error("Insertion failure");
      iemp = effects_.emplace(std::make_unique<KKEnd<KTRAJ>>(reftraj,TDir::backwards,config_.dwt_));
      if(!iemp.second)throw std::runtime_error("Insertion failure");
      // now fit the track
      fit();
    }

  // fit iteration management.  This should include an update/iteration schedule FIXME!
  template <class KTRAJ> void KKTrk<KTRAJ>::fit() {
    while(canIterate()) {
      // create a new status for this iteration
      FitStatus fstat = fitStatus();
      fstat.status_ = FitStatus::unconverged; // by default, assume unconverged
      fstat.ndof_ = -(int)KTRAJ::NParams();
      fstat.chisq_ = 0.0;
      fstat.iter_++;

      bool fitOK(true);
      if(fstat.iter_ > 0) fitOK = update();
      if(fitOK)fitIteration(fstat);
      if(!fitOK)fstat.status_ = FitStatus::failed;
      // record this status in the history
      history_.push_back(fstat);
    }
  }

  // single fit iteration 
  template <class KTRAJ> bool KKTrk<KTRAJ>::fitIteration(FitStatus& fstat) {
    bool fitOK(true);
    // fit forwards then backwards
    auto ieff = effects_.begin();
    // start with empty fit information; each effect will modify this as necessary, and cache what it needs for later processing
    KKData<KTRAJ::NParams()> fitdata;
    while(ieff != effects_.end()){
      // update chisq, NDOF only forwards to save time
      if(ieff->get()->nDOF() > 0){
	fstat.ndof_ += ieff->get()->nDOF();
	// I'm not sure this is the most efficient way to get chisq FIXME!
	fstat.chisq_ += ieff->get()->chisq(fitdata.pData());
	if(isnan(fstat.chisq_)){
	  fitOK = false;
	  break;
	}
      }
      if(!ieff->get()->process(fitdata,TDir::forwards)){
	fitOK = false;
	break;
      }
      ieff++;
    }
    if(fitOK){
    // reset the fit information and process backwards
      fitdata = KKData<KTRAJ::NParams()>();
      auto ieff = effects_.rbegin();
      while(ieff != effects_.rend()){
	if(!ieff->get()->process(fitdata,TDir::backwards)){
	  fitOK = false;
	  break;
	}
	ieff++;
      }
    }
    if(fitOK){
      // convert the fit result into a new trajectory; start with an empty ptraj
      fittraj_ = PKTRAJ(reftraj_.range(),reftraj_.mass(),reftraj_.charge());
      // process forwards, adding pieces as necessary
      for(auto& ieff : effects_) {
	fitOK &= ieff->append(fittraj_);
	if(!fitOK)break;
      }
      if(fitOK) {
	// trim the range to the physical elements (past the end sites
	auto feff = effects_.begin(); feff++;
	auto leff = effects_.rbegin(); leff++;
	fittraj_.range().low() = (*feff)->time() - config_.tbuff_;
	fittraj_.range().high() = (*leff)->time() + config_.tbuff_;
      }
    }
    // update the status
    if(!fitOK || isnan(fstat.chisq_)){
      fstat.status_ = FitStatus::failed;
    } else if(fabs(fstat.chisq_ -fitStatus().chisq_) < config_.convdchisq_) {
      fstat.status_ = FitStatus::converged;
    } else if (fstat.chisq_-fitStatus().chisq_ > config_.divdchisq_) {
      fstat.status_ = FitStatus::diverged;
    } else if (fstat.ndof_ < config_.minndof_){
      fstat.status_ = FitStatus::lowNDOF;
    } else if(oscillating(fstat)){
      fstat.status_ = FitStatus::oscillating;
    }
    fstat.prob_ = TMath::Prob(fstat.chisq_,fstat.ndof_);
    return fitOK;
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

  template<class KTRAJ> bool KKTrk<KTRAJ>::canIterate() const {
    auto const& fstat = fitStatus();
    if( (fstat.status_ == FitStatus::needsfit ||
	  fstat.status_ == FitStatus::unconverged) &&
	fstat.iter_ < config_.maxniter_) {
      return true;
    } else
      return false;
  }

  template<class KTRAJ> bool KKTrk<KTRAJ>::oscillating(FitStatus const& fstat) const {
    // check for oscillation in the history; changing sign and similar values between iterations
    return false; // FIXME!
  }

  struct BFieldDomain {
    TRange range_; // time range of this domain
    Vec3 pdir_; // (average) momentum direction of this domain
  };

  template <class KTRAJ> void KKTrk<KTRAJ>::createFieldDomains() {
    // use the reference trajectory to define the magnetic 'domains' where
    // the trajectory direction isn't changing rapidly.
    //
    //// Maximum time step keeping momentum direction within tolerance
    //  double tstep = reftraj_.timeStep(reftraj_.range().mid(),config_.dtol_);
    //// re-align to a fixed set of covering steps
    //  unsigned nsteps = ciel(reftraj_.range().range()/tstep);
    //  tstep = reftraj_.range().range()/(nsteps-1);
    //  // sample the BField at these points
    //  std::vector<BFieldDomain> domains;
    //  domains.reserve(nsteps);
    //  for(unsigned istep=0;istep<nsteps;istep++){
    //    BFieldDomain domain;
    //    double tlow = reftraj_.range().low() + tstep*istep;
    //    domain.range_ = TRange(tlow,tlow+tstep);
    //    double time = domain.range_.mid();
    //    Vec3 tpos,fvec;
    //    reftraj_.position(time,tpos);
    //    bfield_.fieldVect(fvec,tpos);
    //
    //    // integrate 
    //  }
  }

  template <class KTRAJ> void KKTrk<KTRAJ>::print(std::ostream& ost, int detail) const {
    using std::endl;
    ost <<  "Fit Status:" << fitStatus() << endl;
    ost << "Fit Result";
    fitTraj().print(ost,detail);
    if(detail > 0){
      ost <<  "Fit History:" << endl;
      for(auto const& stat : statusHistory()) ost << stat << endl;
    }
    if(detail > 1) {
      ost << "Reference ";
      refTraj().print(ost,detail);
      ost << "Effects:" << endl;
      for(auto const& eff : effects()) eff.get()->print(ost,detail);
    }
  }
}
#endif
