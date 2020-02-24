#ifndef KinKal_KalSite_hh
#define KinKal_KalSite_hh
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
#include "KinKal/KKHit.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/THit.hh"
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
  };
// struct to define fit status
  struct FitStatus {
    enum status {current=0,needsfit,unconverged,failed}; // fit status
    unsigned niter_; // number of iterations executed;
    status status_; // current status
    FitStatus() : niter_(0), status_(needsfit) {}
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
      // construct from a set of hits and material elements.  Must deweight initial covariance matrix
      KKTrk(PKTRAJ const& reftraj, std::vector<const THit*> hits, Config const& config ); 
      void fit(bool force=false); // process the effects.
      FitStatus status() const { return status_; }
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKCON const& effects() const { return effects_; }
    private:
      // helper functions
      bool buildTraj();
      bool update();
      bool fitIteration();
      // payload
      Config config_; // configuration
      FitStatus status_; // current fit status
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated
      PKTRAJ fittraj_; // result of the current fit
      KKCON effects_; // effects used in this fit, sorted by time
  };

  template <class KTRAJ> KKTrk<KTRAJ>::KKTrk(PKTRAJ const& reftraj, std::vector<const THit*> hits,Config const& config) : 
    config_(config), reftraj_(reftraj), fittraj_(reftraj){
      // create end sites
      auto iemp = effects_.emplace(new KKEnd(reftraj,TDir::forwards,config_.dwt_));
      if(!iemp.second)throw std::runtime_error("Insertion failure");
      iemp = effects_.emplace(new KKEnd(reftraj,TDir::backwards,config_.dwt_));
      if(!iemp.second)throw std::runtime_error("Insertion failure");
      // loop over the hits
      for(auto hit : hits ) {
	// use POCA to find the local piece associated with this hit on the reference
	TPoca<PKTRAJ,TLine> tp(reftraj,hit->sensorTraj());
	if(!tp.usable()) throw std::runtime_error("TPoca failure");
	auto ltraj = reftraj_.nearestPiece(tp.poca0().T());
	// create the hit effects and insert them in the set
	iemp = effects_.emplace(new KKHit(*hit,ltraj));
	if(!iemp.second)throw std::runtime_error("Insertion failure");
      }
    }

  // fit iteration management.  This should cover simulated annealing too FIXME!
  template <class KTRAJ> void KKTrk<KTRAJ>::fit(bool force) {
    // don't fit unless necessary or forced
    if(!force && status_.status_ == FitStatus::current)return;
    // iterate to convergence
    status_.niter_ = 0;
    bool fitOK(true);
    while(fitOK && status_.niter_ < config_.maxniter_ && status_.status_ != FitStatus::current) {
      if(status_.niter_>0) fitOK = update();
      if(fitOK) fitOK = fitIteration();
      // test for convergence/divergence FIXME!
      status_.niter_++;
    }
    if(!fitOK){
      status_.status_ = FitStatus::failed;
    } else if(status_.niter_ == config_.maxniter_) {
      status_.status_ = FitStatus::unconverged;
    } else {
      status_.status_ = FitStatus::current;
    }
  }

  // single fit iteration 
  template <class KTRAJ> bool KKTrk<KTRAJ>::fitIteration() {
  // fit in both directions
    for(TDir tdir=TDir::forwards; tdir != TDir::end; ++tdir){
      if(tdir==TDir::forwards){
	auto ieff = effects_.begin(); ieff++;
	while(ieff != effects_.end()){
	  auto iprev = ieff; iprev--;
	  if(!ieff->get()->process(**iprev,tdir)) return false;
	  ieff++;
	}
      } else {
	auto ieff = effects_.rbegin(); ieff++;
	while(ieff != effects_.rend()){
	  auto iprev = ieff; iprev--;
	  if(!ieff->get()->process(**iprev,tdir)) return false;
	  ieff++;
	}
      }
    }
    // build the fit trajectory
    return buildTraj();
  }

  // convert a final fit into a trajectory
  template <class KTRAJ> bool KKTrk<KTRAJ>::buildTraj() {
  // convert the fit result into a new trajectory.  Start at one end
  // for now, take the fit at one end, this won't work with materials FIXME!
    auto newtraj = reftraj_.front();
    // overwrite the parameters
    auto const& front = *effects_.begin();
    newtraj.params() = front->params(TDir::backwards);
    // create the piecetraj from that
    fittraj_ = PKTRAJ(newtraj);
    return true;
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
