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
#include "KinKal/Fit/HitConstraint.hh"
#include "KinKal/Fit/Material.hh"
#include "KinKal/Fit/BFieldEffect.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/Status.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/TimeDir.hh"
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
      using KKHIT = HitConstraint<KTRAJ>;
      using KKMAT = Material<KTRAJ>;
      using KKEND = TrackEnd<KTRAJ>;
      using KKBFIELD = BFieldEffect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;
      using HITCOL = std::vector<HITPTR>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      using DOMAINCOL = std::vector<TimeRange>;
      using CONFIGCOL = std::vector<Config>;
      struct KKEFFComp { // comparator to sort effects by time
        bool operator()(std::unique_ptr<KKEFF> const& a, std::unique_ptr<KKEFF> const&  b) const {
          if(a.get() != b.get())
            return a->time() < b->time();
          else
            return false;
        }
      };
      using KKEFFCOL = std::vector<std::unique_ptr<KKEFF>> ; // container type for effects
      // construct from a set of hits and passive material crossings
      Track(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, HITCOL& hits, EXINGCOL& exings );
      // extend an existing track with either new configuration, new hits, and/or new material xings
      void extend(Config const& config, HITCOL& hits, EXINGCOL& exings );
      // accessors
      std::vector<Status> const& history() const { return history_; }
      Status const& fitStatus() const { return history_.back(); } // most recent status
      KTRAJ const& seedTraj() const { return seedtraj_; }
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKEFFCOL const& effects() const { return effects_; }
      Config const& config() const { return config_.back(); }
      CONFIGCOL const* configs() const { return config_; }
      BFieldMap const& bfield() const { return bfield_; }
      HITCOL const& hits() const { return hits_; }
      EXINGCOL const& exings() const { return exings_; }
      DOMAINCOL const& domains() const { return domains_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
    protected:
      Track(Config const& cfg, BFieldMap const& bfield, KTRAJ const& seedtraj );
      void fit(HITCOL& hits, EXINGCOL& exings );
    private:
      // helper functions
      TimeRange getRange(HITCOL& hits, EXINGCOL& exings) const;
      void fit(); // process the effects and create the trajectory.  This executes the current schedule
      void update(Status const& fstat, MetaIterConfig const& miconfig);
      void fitIteration(Status& status, MetaIterConfig const& miconfig);
      bool canIterate() const;
      void createEffects( HITCOL& hits, EXINGCOL& exings, DOMAINCOL const& domains);
      void createRefTraj(KTRAJ const& seedtraj,TimeRange const& refrange, DOMAINCOL const& domains);
      void replaceRefTraj(DOMAINCOL const& domains);
      void extendRefTraj(DOMAINCOL const& domains);
      // divide a kinematic trajectory range into magnetic 'domains' within which the BField inhomogeneity effects are within tolerance
      void createDomains(PKTRAJ const& pktraj, TimeRange const& range, std::vector<TimeRange>& ranges, TimeDir tdir=TimeDir::forwards) const;
      // payload
      CONFIGCOL config_; // configuration
      BFieldMap const& bfield_; // magnetic field map
      std::vector<Status> history_; // fit status history; records the current iteration
      KTRAJ seedtraj_; // seed for the fit
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated and the current fit performed
      PKTRAJ fittraj_; // result of the current fit, becomes the reference when the fit is algebraically iterated
      KKEFFCOL effects_; // effects used in this fit, sorted by time
      HITCOL hits_; // hits used in this fit
      EXINGCOL exings_; // material xings used in this fit
      DOMAINCOL domains_; // BField domains used in this fit
  };
  // sub-class constructor, based just on the seed.  It requires added hits to
  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield, KTRAJ const& seedtraj ) :
    bfield_(bfield), seedtraj_(seedtraj) {
      config_.push_back(cfg);
      if(config().schedule().size() ==0)throw std::invalid_argument("Invalid configuration: no schedule");
    }

  // construct from configuration, reference (seed) fit, hits,and materials specific to this fit.
  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield, KTRAJ const& seedtraj,  HITCOL& hits, EXINGCOL& exings) : Track(cfg,bfield,seedtraj) {
    fit(hits,exings);
  }
  template <class KTRAJ> void Track<KTRAJ>::fit(HITCOL& hits, EXINGCOL& exings) {
    // set the seed time based on the min and max time from the inputs
    TimeRange refrange = getRange(hits,exings);
    seedtraj_.setRange(refrange);
    // if correcting for BField effects, define the domains
    DOMAINCOL domains;
    PKTRAJ pkseed(seedtraj_);
    if(config().bfcorr_ ) createDomains(pkseed, refrange, domains);
    // Create the initial reference trajectory
    createRefTraj(seedtraj_,refrange,domains);
    // create the end effects: these help manage the fit
    effects_.reserve(hits.size()+exings.size()+domains.size()+2);
    effects_.emplace_back(std::make_unique<KKEND>(config(), bfield_, reftraj_,TimeDir::forwards));
    effects_.emplace_back(std::make_unique<KKEND>(config(), bfield_, reftraj_,TimeDir::backwards));
    // create all the other effects
    createEffects(hits,exings,domains);
    // now fit the track
    fit();
  }

  // extend an existing track
  template <class KTRAJ> void Track<KTRAJ>::extend(Config const& cfg, HITCOL& hits, EXINGCOL& exings) {
    // update the configuration
    auto const& oldconfig(config());
    config_.push_back(cfg);
    // configuation check
    if(config().schedule().size() ==0)throw std::invalid_argument("Invalid configuration: no schedule");
    // require the existing fit to be usable
    if(!fitStatus().usable())throw std::invalid_argument("Cannot extend unusable fit");
    // find the range of the added information, and extend as needed
    TimeRange exrange = getRange(hits,exings);
    exrange.begin() = std::min(exrange.begin(),reftraj_.range().begin());
    exrange.end() = std::max(exrange.end(),reftraj_.range().end());
    DOMAINCOL domains;
    if(config().bfcorr_ ) {
      // if the previous configuration didn't have domains, then create them for the full reference range
      if(!oldconfig.bfcorr_){
        // create domains for the whole range
        createDomains(reftraj_, exrange,domains);
        // replace the reftraj with one with BField rotations
        replaceRefTraj(domains);
      } else {
        // create domains just for the extensions, and extend the reftraj as needed
        TimeRange exlow(exrange.begin(),reftraj_.range().begin());
        if(exlow.range()>0.0) {
          DOMAINCOL lowdomains;
          createDomains(reftraj_, exlow, lowdomains, TimeDir::backwards);
          extendRefTraj(domains);
          domains.insert(domains.begin(),lowdomains.begin(),lowdomains.end());
        }
        TimeRange exhigh(reftraj_.range().end(),exrange.end());
        if(exhigh.range()>0.0){
          DOMAINCOL highdomains;
          createDomains(reftraj_, exhigh, highdomains,TimeDir::forwards);
          extendRefTraj(highdomains);
          domains.insert(domains.end(),highdomains.begin(),highdomains.end());
        }
      }
    }
    // create the effects for the new info and the new domains
    createEffects(hits,exings,domains);
    // now refit the track
    fit();
  }

  // replace the reference traj with one describing the 'same' trajectory in space, but using the local BField as reference
  template <class KTRAJ> void Track<KTRAJ>::replaceRefTraj(DOMAINCOL const& domains) {
  // create new traj
    PKTRAJ newref;
    // loop over domains
    for(auto const& domain : domains) {
      double dtime = domain.begin();
      // Set the BField to the start of this domain
      auto bf = bfield_.fieldVect(reftraj_.position3(dtime));
      // loop until we're either out of this domain or the piece is out of this domain
      while(dtime < domain.end()){
        // find the nearest piece of the current reftraj
        auto const& oldpiece = reftraj_.nearestPiece(dtime);
        // create a new piece
        KTRAJ newpiece(oldpiece,bf,dtime);
        // set the range as needed
        newpiece.range() = TimeRange(dtime,std::min(domain.end(),oldpiece.range().end()));
        reftraj_.append(newpiece);
        // update the time
        static double epsilon(1e-10);
        dtime = newpiece.range().end()+epsilon; // to avoid boundary
      }
    }
    // actually replace the reftraj
    reftraj_ = newref;
  }

  template <class KTRAJ> void Track<KTRAJ>::extendRefTraj(DOMAINCOL const& domains ) {
    // dummy implementation: need to put in parameter rotation at each domain boundary FIXME!
    if(domains.size() > 0){
      // extend the reftraj range
      TimeRange newrange(std::min(reftraj_.range().begin(),domains.front().begin()),
          std::max(reftraj_.range().end(),domains.back().end()));
      reftraj_.setRange(newrange);
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::createRefTraj(KTRAJ const& seedtraj , TimeRange const& range, DOMAINCOL const& domains ) {
    // if we're making local BField corrections, divide the trajectory into domain pieces.  Each will have equivalent parameters, but relative
    // to the local field
    if(config().bfcorr_ ) {
      if(reftraj_.pieces().size() != 0)throw std::invalid_argument("Initial reference trajectory must be empty");
      if(domains.size() == 0)throw std::invalid_argument("Empty domain");
      for(auto const& domain : domains) {
        // Set the BField to the start of this domain
        auto bf = bfield_.fieldVect(seedtraj.position3(domain.begin()));
        KTRAJ newpiece(seedtraj,bf,domain.begin());
        newpiece.range() = domain;
        reftraj_.append(newpiece);
      }
    } else {
      // use the middle of the range as the nominal BField for this fit:
      double tref = range.mid();
      VEC3 bf = bfield_.fieldVect(seedtraj.position3(tref));
      // create the first piece.  Note this constructor adjusts the parameters according to the local field
      KTRAJ firstpiece(seedtraj,bf,tref);
      firstpiece.range() = range;
      // create the piecewise trajectory from this
      reftraj_ = PKTRAJ(firstpiece);
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::createEffects( HITCOL& hits, EXINGCOL& exings,DOMAINCOL const& domains ) {
    // pre-allocate as needed
    effects_.reserve(effects_.size()+hits.size()+exings.size()+domains.size());
    // append the effects.  First, loop over the hits
    for(auto& hit : hits ) {
      // create the hit effects and insert them in the set
      effects_.emplace_back(std::make_unique<KKHIT>(hit,reftraj_));
    }
    //add material effects
    for(auto& exing : exings) {
      effects_.emplace_back(std::make_unique<KKMAT>(exing,reftraj_));
    }
    // add BField effects
    for( auto const& domain : domains) {
      // create the BField effect for integrated differences over this range
      effects_.emplace_back(std::make_unique<KKBFIELD>(config(),bfield_,domain));
    }
    // sort
    std::sort(effects_.begin(),effects_.end(),KKEFFComp ());
    // store the inputs
    hits_.insert(hits_.end(),hits.begin(),hits.end());
    exings_.insert(exings_.end(),exings.begin(),exings.end());
    domains_.insert(domains_.end(),domains.begin(),domains.end());
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
    if(config().plevel_ > Config::none)print(std::cout, config().plevel_);
  }

  // single algebraic iteration
  template <class KTRAJ> void Track<KTRAJ>::fitIteration(Status& fstat, MetaIterConfig const& miconfig) {
    if(config().plevel_ >= Config::complete)std::cout << "Processing fit iteration " << fstat.iter_ << std::endl;
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
      if(config().plevel_ >= Config::detailed){
        std::cout << "Chisq total " << fstat.chisq_ << " increment " << dchisq << " ";
        ieff->print(std::cout,config().plevel_);
      }
      feff++;
    }
    // reset the fit information and process backwards
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
    // compute parameter difference WRT reference.  Compare in the middle
    auto const& mtraj = fittraj_.nearestPiece(fittraj_.range().mid());
    auto const& rtraj = reftraj_.nearestPiece(fittraj_.range().mid());
    DVEC dpar = mtraj.params().parameters() - rtraj.params().parameters();
    DMAT refwt = rtraj.params().covariance();
    if(!refwt.Invert())throw std::runtime_error("Reference covariance uninvertible");
    double delta = ROOT::Math::Similarity(dpar,refwt);
    double dchisq = fstat.chisq_.chisqPerNDOF() - fitStatus().chisq_.chisqPerNDOF();
    // update status.  Convergence criteria is iteration-dependent.
    if (delta > config().pdchi2_) {
      fstat.status_ = Status::paramsdiverged;
      // skip divergence comparsion in first iteration after a meta-iteration, as that
      // is affected by the change in temperature
    } else if (fstat.iter_ > 0 && dchisq > config().divdchisq_) {
      fstat.status_ = Status::diverged;
    } else if (fstat.chisq_.nDOF() < (int)config().minndof_){
      fstat.status_ = Status::lowNDOF;
    } else if(fabs(dchisq) < config().convdchisq_) {
      fstat.status_ = Status::converged;
    } else
      fstat.status_ = Status::unconverged;
  }

  // update between iterations
  template <class KTRAJ> void Track<KTRAJ>::update(Status const& fstat, MetaIterConfig const& miconfig) {
    if(fstat.iter_ < 0) { // 1st iteration of a meta-iteration: update the state
      if(miconfig.miter_ > 0){ // if this isn't the 1st meta-iteration, swap the fit trajectory to the reference
        reftraj_ = fittraj_;
        for(auto& ieff : effects_ ) ieff->update(reftraj_,miconfig);
      }
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
  // divide a trajectory into magnetic 'domains' used to apply the BField corrections
  template<class KTRAJ> void Track<KTRAJ>::createDomains(PKTRAJ const& pktraj, TimeRange const& range, std::vector<TimeRange>& ranges,
  TimeDir tdir) const {
    double tstart;
    tstart = range.begin();
    do {
      // see how far we can go on the current traj before the BField change causes the momentum estimate to go out of tolerance
      auto const& ktraj = pktraj.nearestPiece(tstart);
      double tend = bfield_.rangeInTolerance(ktraj,tstart,config().tol_);
      ranges.emplace_back(tstart,tend);
      // start the next domain and the end of this one
      tstart = tend;
    } while(tstart < range.end());
  }

  template<class KTRAJ> TimeRange Track<KTRAJ>::getRange(HITCOL& hits, EXINGCOL& exings) const {
    double tmin = std::numeric_limits<double>::max();
    double tmax = -std::numeric_limits<double>::max();
    for(auto const& hit : hits){
      tmin = std::min(tmin,hit->time());
      tmax = std::max(tmax,hit->time());
    }
    for(auto const& exing : exings){
      tmin = std::min(tmin,exing->crossingTime());
      tmax = std::max(tmax,exing->crossingTime());
    }
    // add a buffer to the time.  This must cover the uncertainty on t0 as the fit iterates
    return TimeRange(tmin-config().tbuff_,tmax+config().tbuff_);
  }
}
#endif
