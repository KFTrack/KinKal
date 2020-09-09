#ifndef KinKal_Track_hh
#define KinKal_Track_hh
//
//  Primary class of the Kinematic Kalman fit.  This class owns the state describing
//  the fit inputs (hits, material interactions, BFieldMap corrections, etc), the result of the fit,
//  and the methods for computing it.  The fit result is expressed as a piecewise kinematic covariant
//  particle trajectory, providing position, momentum etc information about the particle with covariance
//  as a function of physical time.
//
//  Track is templated on a simple kinematic trajectory class representing the 1-dimensional path and
//  momentum of a particle traveling through empty space in a constant magnetic field, as a function of time.
//  The piecewise kinematic trajectory fit result is expressed as a time sequence of these simple trajectory objects.
//  Material effects and spatial variation of magnetic fields are modeled through changes between adjacent simple
//  trajectories.
//  To instantiate Track the particle trajectory class must satisfy a geometric, kinematic, and parametric interface.
//  The geometric interface includes functions for position, direction, etc.
//  The kinematic interface includes functions for velocity, momentum, etc.
//  The parametric interface includes functions for parameter values, covariance, derivatives, etc.
//  An example is the LoopHelix.hh or CentralHelix.hh classes.
//
//  The Parameters object provides a minimal basis from which the geometric and kinematic properties of the particle as a function
//  of time can be computed.  For instance, a kinematic helix in space requires a Parameters instance with 6 parameters.  The physical
//  interpretation of the Parameters payload is made in the kinematic trajectory class.
//
//  Track uses the root SVector and SMatrix classes for algebraic manipulation, and GenVector classes for geometric and
//  kinematic particle descriptions, both part of the root Math package.  These are described on the root website https://root.cern.ch/root/html608/namespaceROOT_1_1Math.html
//
//  The underlying processing model is a progressive BLUE fit first used in the geometric track fit implementation used by the BaBar
//  collaboration, described in "D.N. Brown, E.A. Charles, D.A. Roberts, The BABAR track fitting algorithm, Proceedings of CHEP 2000, Padova, Italy, 2000"
//
//  Track is constructed from a configuration object which can be shared between many instances, and a unique set of hit and
//  material interactions.  The configuration object controls the fit iteration convergence testing, including simulated
//  annealing and interactions with the external environment such as the material model and the magnetic field map.
//  The fit is performed on construction.
//
//  The KinKal package is licensed under Adobe v2, and is hosted at https://github.com/KFTrack/KinKal.git
//  David N. Brown, Lawrence Berkeley National Lab
//
#include "KinKal/ParticleTrajectory.hh"
#include "KinKal/FitData.hh"
#include "KinKal/Effect.hh"
#include "KinKal/TrackEnd.hh"
#include "KinKal/MaterialHit.hh"
#include "KinKal/Hit.hh"
#include "KinKal/Material.hh"
#include "KinKal/BField.hh"
#include "KinKal/ClosestApproach.hh"
#include "KinKal/Hit.hh"
#include "KinKal/Config.hh"
#include "KinKal/FitStatus.hh"
#include "KinKal/BFieldMap.hh"
#include "KinKal/BFieldUtils.hh"
#include "TMath.h"
#include <set>
#include <vector>
#include <iterator>
#include <memory>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class Track {
    public:
      using KKEFF = Effect<KTRAJ>;
      using KKHIT = Hit<KTRAJ>;
      using KKMAT = Material<KTRAJ>;
      using KKEND = TrackEnd<KTRAJ>;
      using KKMHIT = MaterialHit<KTRAJ>;
      using KKBFIELD = BField<KTRAJ>;
      using KKCONFIGPTR = std::shared_ptr<Config>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using THIT = DetectorHit<KTRAJ>;
      using THITPTR = std::shared_ptr<THIT>;
      using THITCOL = std::vector<THITPTR>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
      using DXINGCOL = std::vector<DXINGPTR>;
      struct KKEFFComp { // comparator to sort effects by time
	bool operator()(std::unique_ptr<KKEFF> const& a, std::unique_ptr<KKEFF> const&  b) const {
	  if(a.get() != b.get())
	    return a->time() < b->time();
	  else
	    return false;
	}
      };
      typedef std::vector<std::unique_ptr<KKEFF>> KKEFFCOL; // container type for effects
      // construct from a set of hits and passive material crossings
      Track(KKCONFIGPTR kkconfig, KTRAJ const& seedtraj, THITCOL& thits, DXINGCOL& dxings ); 
      void fit(); // process the effects.  This creates the fit
      // accessors
      std::vector<FitStatus> const& history() const { return history_; }
      FitStatus const& fitStatus() const { return history_.back(); } // most recent status
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKEFFCOL const& effects() const { return effects_; }
      Config const& config() const { return *kkconfig_; }
      THITCOL const& timeHits() const { return thits_; } 
      void print(std::ostream& ost=std::cout,int detail=0) const;
    private:
      // helper functions
      void update(FitStatus const& fstat, MetaIterConfig const& miconfig);
      void fitIteration(FitStatus& status, MetaIterConfig const& miconfig);
      bool canIterate() const;
      bool oscillating(FitStatus const& status, MetaIterConfig const& miconfig) const;
      void createRefTraj(KTRAJ const& seedtraj);
      // payload
      KKCONFIGPTR kkconfig_; // configuration
      std::vector<FitStatus> history_; // fit status history; records the current iteration
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated and the current fit performed
      PKTRAJ fittraj_; // result of the current fit, becomes the reference when the fit is algebraically iterated
      KKEFFCOL effects_; // effects used in this fit, sorted by time
      THITCOL thits_; // shared collection of hits
  };

// construct from configuration, reference (seed) fit, hits,and materials specific to this fit.  Note that hits
// can contain associated materials.
  template <class KTRAJ> Track<KTRAJ>::Track(KKCONFIGPTR kkconfig, KTRAJ const& seedtraj,  THITCOL& thits, DXINGCOL& dxings) : 
    kkconfig_(kkconfig), thits_(thits) {
      // Create the initial reference traj.  This also divides the range into domains of ~constant BFieldMap and creates correction effects for inhomogeneity
      createRefTraj(seedtraj);
      // create the effects.  First, loop over the hits
      for(auto& thit : thits_ ) {
	// create the hit effects and insert them in the set
	// if there's associated material, create a combined material and hit effect, otherwise just a hit effect
	if(kkconfig_->addmat_ && thit->hasMaterial()){
	  effects_.emplace_back(std::make_unique<KKMHIT>(thit,reftraj_));
	} else{ 
	  effects_.emplace_back(std::make_unique<KKHIT>(thit,reftraj_));
	}
      }
      //add pure material effects
      if(kkconfig_->addmat_){
	for(auto& dxing : dxings) {
	  effects_.emplace_back(std::make_unique<KKMAT>(dxing,reftraj_));
	}
      }
      // preliminary sort; this makes sure the range is accurate when computing BFieldMap corrections
      std::sort(effects_.begin(),effects_.end(),KKEFFComp ());
      // reset the range 
      reftraj_.setRange(TimeRange(std::min(reftraj_.range().begin(),effects_.begin()->get()->time() - config().tbuff_),
	    std::max(reftraj_.range().end(),effects_.rbegin()->get()->time() + config().tbuff_)));
      // create the end effects: these help manage the fit
      effects_.emplace_back(std::make_unique<KKEND>(reftraj_,TimeDir::forwards,config().dwt_));
      effects_.emplace_back(std::make_unique<KKEND>(reftraj_,TimeDir::backwards,config().dwt_));
      // now fit the track
      fit();
      if(kkconfig_->plevel_ > Config::none)print(std::cout, kkconfig_->plevel_);
    }

  // fit iteration management 
  template <class KTRAJ> void Track<KTRAJ>::fit() {
    // execute the schedule of meta-iterations
    for(auto imiconfig=config().schedule().begin(); imiconfig != config().schedule().end(); imiconfig++){
      auto miconfig  = *imiconfig;
      miconfig.miter_  = std::distance(config().schedule().begin(),imiconfig);
      // algebraic convergence iteration
      FitStatus fstat(miconfig.miter_);
      history_.push_back(fstat);
      if(kkconfig_->plevel_ >= Config::basic)std::cout << "Processing fit meta-iteration " << miconfig << std::endl;
      while(canIterate()) {
	// catch exceptions and record them in the status
	try {
	  update(fstat,miconfig);
	  fitIteration(fstat,miconfig);
	} catch (std::exception const& error) {
	  fstat.status_ = FitStatus::failed;
	  fstat.comment_ = error.what();
	}
	// record this status in the history
	history_.push_back(fstat);
      }
      if(!fstat.usable())break;
    }
  }

  // single algebraic iteration 
  template <class KTRAJ> void Track<KTRAJ>::fitIteration(FitStatus& fstat, MetaIterConfig const& miconfig) {
    if(kkconfig_->plevel_ >= Config::complete)std::cout << "Processing fit iteration " << fstat.iter_ << std::endl;
    // reset counters
    fstat.chisq_ = 0.0;
    fstat.ndof_ = -(int)NParams();
    fstat.iter_++;
    // fit in both directions (order doesn't matter)
    auto feff = effects_.begin();
    // start with empty fit information; each effect will modify this as necessary, and cache what it needs for later processing
    FitData ffitdata;
    while(feff != effects_.end()){
      auto ieff = feff->get();
      // update chisquared; only needed forwards
      fstat.ndof_ += ieff->nDOF();
      double dchisq = ieff->chisq(ffitdata.pData());
      fstat.chisq_ += dchisq;
      // process
      ieff->process(ffitdata,TimeDir::forwards);
      if(kkconfig_->plevel_ >= Config::detailed){
	std::cout << "Chisq total " << fstat.chisq_ << " increment " << dchisq << " ";
	ieff->print(std::cout,kkconfig_->plevel_);
      }
      feff++;
    }
    fstat.prob_ = TMath::Prob(fstat.chisq_,fstat.ndof_);
    // reset the fit information and process backwards (the order does not matter)
    FitData bfitdata;
    auto beff = effects_.rbegin();
    while(beff != effects_.rend()){
      auto ieff = beff->get();
      ieff->process(bfitdata,TimeDir::backwards);
      beff++;
    }
    // convert the fit result into a new trajectory; start with an empty ptraj
    fittraj_ = PKTRAJ();
    // process forwards, adding pieces as necessary
    for(auto& ieff : effects_) {
      ieff->append(fittraj_);
    }
    // trim the range to the physical elements (past the end sites)
    feff = effects_.begin(); feff++;
    beff = effects_.rbegin(); beff++;
    fittraj_.front().range().begin() = (*feff)->time() - config().tbuff_;
    fittraj_.back().range().end() = (*beff)->time() + config().tbuff_;
    // update status.  Convergence criteria is iteration-dependent
    double dchisq = (fstat.chisq_ -fitStatus().chisq_)/fstat.ndof_;
    if (fstat.ndof_ < config().minndof_){
      fstat.status_ = FitStatus::lowNDOF;
    } else if(fabs(dchisq) < miconfig.convdchisq_) {
      fstat.status_ = FitStatus::converged;
    } else if (dchisq > miconfig.divdchisq_) {
      fstat.status_ = FitStatus::diverged;
    } else if(oscillating(fstat,miconfig)){
      fstat.status_ = FitStatus::oscillating;
    } else
      fstat.status_ = FitStatus::unconverged;
  }

  // update between iterations 
  template <class KTRAJ> void Track<KTRAJ>::update(FitStatus const& fstat, MetaIterConfig const& miconfig) {
    if(fstat.iter_ < 0) { // 1st iteration of a meta-iteration: update the state
      if(miconfig.miter_ > 0)// if this isn't the 1st meta-iteration, swap the fit trajectory to the reference
	reftraj_ = fittraj_;
      for(auto& ieff : effects_ ) ieff->update(reftraj_,miconfig);
    } else {
      //swap the fit trajectory to the reference
      reftraj_ = fittraj_;
      // update the effects to use the new reference
      for(auto& ieff : effects_) ieff->update(reftraj_);
    }
    // sort the effects by time
    std::sort(effects_.begin(),effects_.end(),KKEFFComp ());
  }

  template<class KTRAJ> bool Track<KTRAJ>::canIterate() const {
    return fitStatus().needsFit() && fitStatus().iter_ < config().maxniter_;
  }

  template<class KTRAJ> bool Track<KTRAJ>::oscillating(FitStatus const& fstat, MetaIterConfig const& miconfig) const {
    if(history_.size()>=3 &&history_[history_.size()-3].miter_ == fstat.miter_ ){
      double d1 = fstat.chisq_ - history_.back().chisq_;
      double d2 = fstat.chisq_ - history_[history_.size()-2].chisq_;
      double d3 = history_.back().chisq_ - history_[history_.size()-3].chisq_;
      if(d1*d2 < 0.0 && d1*d3 > 0.0 && fabs(d1) - fabs(d2) < miconfig.oscdchisq_ && fabs(fabs(d2) - fabs(d3)) < miconfig.oscdchisq_) return true;
    }
    return false;
  }

  template <class KTRAJ> void Track<KTRAJ>::createRefTraj(KTRAJ const& seedtraj ) {
  // initialize the reftraj.  Range is taken from the seed
    double tstart = seedtraj.range().begin();
    VEC3 bf;
    if(kkconfig_->bfcorr_ == Config::variable) {
      // initialize BNom at the start of the range. it will change with each piece
      bf = kkconfig_->bfield_.fieldVect(seedtraj.position(tstart)); 
      // recast the seed parameters so they give the same state vector with the field at the starting point
      KTRAJ piece(seedtraj,bf,tstart);
      reftraj_ = PKTRAJ(piece);
    } else {
      // use the seed BFieldMap, fixed for the whole fit
      reftraj_ = PKTRAJ(seedtraj); // the initial ref traj is just the seed.
      bf = reftraj_.bnom(tstart); // freeze the nominal BFieldMap to be the one from the seed
    }
    if(kkconfig_->bfcorr_ != Config::nocorr) { 
      // divide up this traj into domains.  start the field at the start of the seed trajectory
      TimeRange drange(tstart,reftraj_.range().end());
      while(drange.begin() < reftraj_.range().end()){
	// see how far we can go before the BFieldMap change cause the traj to go out of tolerance
	drange.end() = BFieldUtils::rangeInTolerance(drange.begin(),kkconfig_->bfield_, reftraj_, kkconfig_->tol_);
	if(kkconfig_->bfcorr_ == Config::variable) {
	  // create the next piece and append.  The domain transition is set to the middle of the integration range, so the effects coincide
	  double tdomain = drange.mid();
	  bf = kkconfig_->bfield_.fieldVect(reftraj_.position(tdomain));
	  // update the parameters to correspond to the same state but referencing the local field.
	  // this allows the effects built on this traj to reference the correct parameterization
	  KTRAJ newpiece(reftraj_.back(),bf,tdomain);
	  newpiece.range() = TimeRange(tdomain,std::max(drange.end(),reftraj_.range().end()));
	  reftraj_.append(newpiece);
	}
	// create the BFieldMap effect for integrated differences over this range
	effects_.emplace_back(std::make_unique<KKBFIELD>(kkconfig_->bfield_,reftraj_,drange,kkconfig_->bfcorr_));
	drange.begin() = drange.end(); // reset for next domain
      }
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::print(std::ostream& ost, int detail) const {
    using std::endl;
    if(detail == Config::minimal) 
      ost <<  fitStatus();
    else {
      ost <<  "Fit History " << endl;
      for(auto const& stat : history_) ost << stat << endl;
    }
    ost << " Fit Result ";
    fitTraj().print(ost,detail);
    if(detail > Config::basic) {
      ost << " Reference ";
      refTraj().print(ost,detail-2);
    }
    if(detail > Config::complete) {
      ost << " Effects " << endl;
      for(auto const& eff : effects()) eff.get()->print(ost,detail-3);
    }
  }

}
#endif
