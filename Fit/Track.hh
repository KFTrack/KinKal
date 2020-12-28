#ifndef KinKal_Track_hh
#define KinKal_Track_hh
//
//  Primary class of the Kinematic Kalman fit.  This class owns the state describing
//  the fit inputs (measurements, material interactions, BField corrections, etc), the result of the fit,
//  and the methods for computing it.  The fit result is expressed as a piecewise kinematic covariant
//  particle trajectory, providing position, momentum etc information about the particle with covariance
//  as a function of physical time.
//
//  Track is templated on a simple kinematic trajectory class representing the 1-dimensional path and
//  momentum of a particle traveling through empty space in a constant magnetic field, as a function of time.
//  Material effects and spatial variation of magnetic fields are modeled through changes between adjacent simple  trajectories.
//  The particle trajectory is expressed as a piecewise sequence of these simple trajectory objects, joined
//  at specific times, providing a continuous (in time) description of particle position and momentum.
//  To instantiate Track the kinematic trajectory class must satisfy a geometric, kinematic, and parametric interface.
//  The geometric interface includes functions for position, direction, etc.
//  The kinematic interface includes functions for velocity, momentum, etc.
//  The parametric interface includes functions for parameter values, covariance, derivatives, etc.
//  Fully functional examples are provided, including LoopHelix.hh, CentralHelix, and KinematicLine classes.
//
//  The Parameters object provides a minimal basis from which the geometric and kinematic properties of the particle as a function
//  of time can be computed.  The physical interpretation of the Parameters payload is made in the kinematic trajectory class.
//
//  Track uses the root SVector and SMatrix classes for algebraic manipulation, and GenVector classes for geometric and
//  kinematic particle descriptions, both part of the root Math package.  These are described on the root website https://root.cern.ch/root/html608/namespaceROOT_1_1Math.html
//
//  The underlying processing model is a progressive BLUE fit first used in the geometric track fit implementation used by the BaBar
//  collaboration, described in "D.N. Brown, E.A. Charles, D.A. Roberts, The BABAR track fitting algorithm, Proceedings of CHEP 2000, Padova, Italy, 2000"
//
//  Track is constructed from a configuration object which can be shared between many instances, and a unique set of measurements and
//  material interactions.  The configuration object controls the fit iteration convergence testing, including simulated
//  annealing and interactions with the external environment such as the material model and the magnetic field map.
//  The fit is performed on construction.
//
//  The KinKal package is licensed under Adobe v2, and is hosted at https://github.com/KFTrack/KinKal.git
//  David N. Brown, Lawrence Berkeley National Lab
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Fit/FitState.hh"
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Fit/TrackEnd.hh"
#include "KinKal/Fit/Constraint.hh"
#include "KinKal/Fit/Material.hh"
#include "KinKal/Fit/BFieldEffect.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/Status.hh"
#include "KinKal/Detector/BFieldUtils.hh"
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
      using KKHIT = Constraint<KTRAJ>;
      using KKMAT = Material<KTRAJ>;
      using KKEND = TrackEnd<KTRAJ>;
      using KKBFIELD = BFieldEffect<KTRAJ>;
      using KKCONFIGPTR = std::shared_ptr<Config>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using MEAS = Hit<KTRAJ>;
      using MEASPTR = std::shared_ptr<MEAS>;
      using MEASCOL = std::vector<MEASPTR>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
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
      Track(KKCONFIGPTR config, KTRAJ const& seedtraj, MEASCOL& thits, EXINGCOL& dxings ); 
      void fit(); // process the effects.  This creates the fit
      // accessors
      std::vector<Status> const& history() const { return history_; }
      Status const& fitStatus() const { return history_.back(); } // most recent status
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKEFFCOL const& effects() const { return effects_; }
      Config const& config() const { return *config_; }
      MEASCOL const& timeHits() const { return thits_; } 
      void print(std::ostream& ost=std::cout,int detail=0) const;
    private:
      // helper functions
      void update(Status const& fstat, MetaIterConfig const& miconfig);
      void fitIteration(Status& status, MetaIterConfig const& miconfig);
      bool canIterate() const;
      bool oscillating(Status const& status, MetaIterConfig const& miconfig) const;
      void createRefTraj(KTRAJ const& seedtraj);
      // payload
      KKCONFIGPTR config_; // configuration
      std::vector<Status> history_; // fit status history; records the current iteration
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated and the current fit performed
      PKTRAJ fittraj_; // result of the current fit, becomes the reference when the fit is algebraically iterated
      KKEFFCOL effects_; // effects used in this fit, sorted by time
      MEASCOL thits_; // shared collection of hits
  };

// construct from configuration, reference (seed) fit, hits,and materials specific to this fit.  Note that hits
// can contain associated materials.
  template <class KTRAJ> Track<KTRAJ>::Track(KKCONFIGPTR cfg, KTRAJ const& seedtraj,  MEASCOL& thits, EXINGCOL& dxings) : 
    config_(cfg), thits_(thits) {
      // Create the initial reference traj.  This also divides the range into domains of ~constant BField and creates correction effects for inhomogeneity
      createRefTraj(seedtraj);
      // create the effects.  First, loop over the hits
      for(auto& thit : thits_ ) {
	// create the hit effects and insert them in the set
	// if there's associated material, create an effect for that too
	effects_.emplace_back(std::make_unique<KKHIT>(thit,reftraj_));
	if(config_->addmat_ && thit->hasMaterial())
	  effects_.emplace_back(std::make_unique<KKMAT>(thit->detXingPtr(),reftraj_));
      }
      //add pure material effects
      if(config_->addmat_){
	for(auto& dxing : dxings) {
	  effects_.emplace_back(std::make_unique<KKMAT>(dxing,reftraj_));
	}
      }
      // preliminary sort; this makes sure the range is accurate when computing BField corrections
      std::sort(effects_.begin(),effects_.end(),KKEFFComp ());
      // reset the range 
      reftraj_.setRange(TimeRange(std::min(reftraj_.range().begin(),effects_.begin()->get()->time() - config().tbuff_),
	    std::max(reftraj_.range().end(),effects_.rbegin()->get()->time() + config().tbuff_)));
      // create the end effects: these help manage the fit
      effects_.emplace_back(std::make_unique<KKEND>(config(), reftraj_,TimeDir::forwards));
      effects_.emplace_back(std::make_unique<KKEND>(config(), reftraj_,TimeDir::backwards));
      // now fit the track
      fit();
      if(config_->plevel_ > Config::none)print(std::cout, config_->plevel_);
    }

  // fit iteration management 
  template <class KTRAJ> void Track<KTRAJ>::fit() {
    // execute the schedule of meta-iterations
    for(auto imiconfig=config().schedule().begin(); imiconfig != config().schedule().end(); imiconfig++){
      auto miconfig  = *imiconfig;
      miconfig.miter_  = std::distance(config().schedule().begin(),imiconfig);
      // algebraic convergence iteration
      Status fstat(miconfig.miter_);
      history_.push_back(fstat);
      if(config_->plevel_ >= Config::basic)std::cout << "Processing fit meta-iteration " << miconfig << std::endl;
      while(canIterate()) {
	// catch exceptions and record them in the status
	try {
	  update(fstat,miconfig);
	  fitIteration(fstat,miconfig);
	} catch (std::exception const& error) {
	  fstat.status_ = Status::failed;
	  fstat.comment_ = error.what();
	}
	// record this status in the history
	history_.push_back(fstat);
      }
      if(!fstat.usable())break;
    }
  }

  // single algebraic iteration 
  template <class KTRAJ> void Track<KTRAJ>::fitIteration(Status& fstat, MetaIterConfig const& miconfig) {
    if(config_->plevel_ >= Config::complete)std::cout << "Processing fit iteration " << fstat.iter_ << std::endl;
    // reset counters
    fstat.chisq_ = Chisq(0.0, -(int)NParams());
    fstat.iter_++;
    // fit in both directions (order doesn't matter)
    auto feff = effects_.begin();
    // start with empty fit information; each effect will modify this as necessary, and cache what it needs for later processing
    FitState forwardstate;
    while(feff != effects_.end()){
      auto ieff = feff->get();
      // update chisquared increment WRT the current state: only needed forwards
      Chisq dchisq = ieff->chisq(forwardstate.pData());
      fstat.chisq_ += dchisq;
      // process
      ieff->process(forwardstate,TimeDir::forwards);
      if(config_->plevel_ >= Config::detailed){
	std::cout << "Chisq total " << fstat.chisq_ << " increment " << dchisq << " ";
	ieff->print(std::cout,config_->plevel_);
      }
      feff++;
    }
    // reset the fit information and process backwards (the order does not matter)
    FitState backwardstate;
    auto beff = effects_.rbegin();
    while(beff != effects_.rend()){
      auto ieff = beff->get();
      ieff->process(backwardstate,TimeDir::backwards);
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
    // update status.  Convergence criteria is iteration-dependent.
    double dchisq = fstat.chisq_.chisqPerNDOF() - fitStatus().chisq_.chisqPerNDOF();
    if (fstat.chisq_.nDOF() < config().minndof_){
      fstat.status_ = Status::lowNDOF;
    } else if(fabs(dchisq) < miconfig.convdchisq_) {
      fstat.status_ = Status::converged;
    } else if (dchisq > miconfig.divdchisq_) {
      fstat.status_ = Status::diverged;
    } else if(oscillating(fstat,miconfig)){
      fstat.status_ = Status::oscillating;
    } else
      fstat.status_ = Status::unconverged;
  }

  // update between iterations 
  template <class KTRAJ> void Track<KTRAJ>::update(Status const& fstat, MetaIterConfig const& miconfig) {
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

  template<class KTRAJ> bool Track<KTRAJ>::oscillating(Status const& fstat, MetaIterConfig const& miconfig) const {
    if(history_.size()>=3 &&history_[history_.size()-3].miter_ == fstat.miter_ ){
      double d1 = fstat.chisq_.chisqPerNDOF() - history_.back().chisq_.chisqPerNDOF();
      double d2 = fstat.chisq_.chisqPerNDOF() - history_[history_.size()-2].chisq_.chisqPerNDOF();
      double d3 = history_.back().chisq_.chisqPerNDOF() - history_[history_.size()-3].chisq_.chisqPerNDOF();
      if(d1*d2 < 0.0 && d1*d3 > 0.0 && fabs(d1) - fabs(d2) < miconfig.oscdchisq_ && fabs(fabs(d2) - fabs(d3)) < miconfig.oscdchisq_) return true;
    }
    return false;
  }

  template <class KTRAJ> void Track<KTRAJ>::createRefTraj(KTRAJ const& seedtraj ) {
    if(config_->bfcorr_ != Config::nocorr) {
    // find the nominal BField.  This can be fixed or variable
      VEC3 bf;
      double tstart = seedtraj.range().begin();
      if(config_->bfcorr_ == Config::fixed) // fixed field: take the middle of the range
	bf = config_->bfield_.fieldVect(seedtraj.position3(seedtraj.range().mid()));
      else // this will change with piece: start with the begining
	bf = config_->bfield_.fieldVect(seedtraj.position3(tstart));
	// create the first piece
      KTRAJ newpiece(seedtraj,bf,tstart);
      reftraj_ = PKTRAJ(newpiece);
      // divide the range up into magnetic 'domains'.  start with the full range
      double tend = tstart;
      do {
	// see how far we can go on the current traj before the BField change causes it to go out of tolerance
	// that defines the end of this domain
	tend = BFieldUtils::rangeInTolerance(tstart,config_->bfield_, reftraj_, config_->tol_);
	// for local correction there is also tolerance coming from 2nd order terms in the rotation of the BField: this is proportional
	// to the lever arm.
	if(config_->localBFieldCorr()){
	  double dx;
	  do{
	    auto epos = reftraj_.position3(tend);
	    auto ebf = config_->bfield_.fieldVect(epos);
	    dx = epos.R()*(1.0-bf.Dot(ebf)/(bf.R()*ebf.R())); // there may be magnitude-based 2nd order terms too TODO
	    if(dx > config_->tol_){
	      double factor = std::min(0.9,0.9*config_->tol_/dx);
	      // decrease the time step
	      tend = tstart + factor*(tend-tstart);
	    }
	  } while(dx > config_->tol_);
	}
	// create the BField effect for integrated differences over this range
	effects_.emplace_back(std::make_unique<KKBFIELD>(config(),reftraj_,TimeRange(tstart,tend)));
	// if we're using a local BField correction, create a new piece that uses the local BField
	if(tend < reftraj_.range().end() && config_->localBFieldCorr()) {
	  // update the BF for the next piece: it is at the end of this one
	  bf = config_->bfield_.fieldVect(reftraj_.position3(tend));
	  // update the trajectory parameters to correspond to the same particle state but referencing the local field.
	  // this allows the effects built on this traj to reference the correct parameterization
	  KTRAJ newpiece(reftraj_.back(),bf,tend);
	  newpiece.range() = TimeRange(tend,reftraj_.range().end());
	  reftraj_.append(newpiece);
	}
	// prepare for the next domain
	tstart = tend;
      } while(tstart < reftraj_.range().end());
    } else {
      // use the seed BField, fixed for the whole fit
      reftraj_ = PKTRAJ(seedtraj); // the initial ref traj is just the seed.  The nominal BField is taken from the seed
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
