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
#include "KinKal/Fit/Measurement.hh"
#include "KinKal/Fit/Material.hh"
#include "KinKal/Fit/BField.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/ExtraConfig.hh"
#include "KinKal/Fit/Status.hh"
#include "KinKal/Fit/Domain.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/TimeDir.hh"
#include "TMath.h"
#include <set>
#include <vector>
#include <deque>
#include <list>
#include <array>
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
      struct KKEFFComp { // comparator to sort effects by time
        bool operator()(std::unique_ptr<KKEFF> const& a, std::unique_ptr<KKEFF> const&  b) const {
          if(a.get() != b.get())
            return a->time() < b->time();
          else
            return false;
        }
      };
      using KKEFFCOL = std::list<std::unique_ptr<KKEFF>>;
      using KKEFFFWD = typename KKEFFCOL::iterator;
      using KKEFFREV = typename KKEFFCOL::reverse_iterator;
      using KKEFFFWDBND = std::array<KKEFFFWD,2>;
      using KKEFFREVBND = std::array<KKEFFREV,2>;
      using KKMEAS = Measurement<KTRAJ>;
      using KKMAT = Material<KTRAJ>;
      using KKBFIELD = BField<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using PTRAJPTR = std::unique_ptr<PTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;
      using HITCOL = std::vector<HITPTR>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      using DOMAINCOL = std::set<Domain>;
      using CONFIGCOL = std::vector<Config>;
      using FitStateArray = std::array<FitState,2>;
      // construct from a set of hits and passive material crossings
      Track(Config const& config, BFieldMap const& bfield, PTRAJ const& seedtraj, HITCOL& hits, EXINGCOL& exings );
      // extend an existing track with either new configuration, new hits, and/or new material xings
      void extend(Config const& config, HITCOL& hits, EXINGCOL& exings );
      // extrapolate the fit with the given config until the given predicate is satisfied. This function requires
      // the fit be valid, otherwise the return code is false.  If successful the status, domains, and trajectory of the fit are updated
      // Note that the actual fit itself is unchanged
      template <class XTEST> bool extrapolate(ExtraConfig const& xconfig, XTEST const& XTest);
      // accessors
      std::vector<Status> const& history() const { return history_; }
      Status const& fitStatus() const { return history_.back(); } // most recent status
      PTRAJ const& seedTraj() const { return seedtraj_; }
      PTRAJ const& fitTraj() const { return *fittraj_; }
      KKEFFCOL const& effects() const { return effects_; }
      Config const& config() const { return config_.back(); }
      CONFIGCOL const& configs() const { return config_; }
      BFieldMap const& bfield() const { return bfield_; }
      HITCOL const& hits() const { return hits_; }
      EXINGCOL const& exings() const { return exings_; }
      DOMAINCOL const& domains() const { return domains_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
    protected:
      Track(Config const& cfg, BFieldMap const& bfield, PTRAJ const& seedtraj );
      void fit(HITCOL& hits, EXINGCOL& exings );
    private:
      // helper functions
      TimeRange getRange(HITCOL& hits, EXINGCOL& exings) const;
      void fit(); // process the effects and create the trajectory.  This executes the current schedule
      void setBounds(KKEFFFWDBND& fwdbnds, KKEFFREVBND& revbnds);
      void iterate(MetaIterConfig const& miconfig);
      void setStatus(PTRAJPTR& ptrajptr);
      void initFitState(FitStateArray& states, TimeRange const& fitrange, double dwt=1.0);
      bool canIterate() const;
      void createEffects( HITCOL& hits, EXINGCOL& exings, DOMAINCOL const& domains);
      void createTraj(PTRAJ const& seedtraj,TimeRange const& refrange, DOMAINCOL const& domains);
      void replaceDomains(DOMAINCOL const& domains);
      void extendTraj(DOMAINCOL const& domains);
      void processEnds();
      // add a single domain within the tolerance and extend the fit in the specified direction.
      void addDomain(Domain const& domain,TimeDir const& tdir);
      auto& status() { return history_.back(); } // most recent status
                                                 // divide a kinematic trajectory range into magnetic 'domains' within which the BField inhomogeneity effects are within tolerance
      void createDomains(PTRAJ const& ptraj, TimeRange const& range, DOMAINCOL& domains, double tol) const;
      // payload
      CONFIGCOL config_; // configuration
      BFieldMap const& bfield_; // magnetic field map
      std::vector<Status> history_; // fit status history; records the current iteration
      PTRAJ seedtraj_; // seed for the fit
      PTRAJPTR fittraj_; // result of the current fit
      KKEFFCOL effects_; // effects used in this fit, sorted by time
      HITCOL hits_; // hits used in this fit
      EXINGCOL exings_; // material xings used in this fit
      DOMAINCOL domains_; // BField domains used in this fit, contiguous and sorted by time
  };
  // sub-class constructor, based just on the seed.  It requires added hits to create a functional track
  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield, PTRAJ const& seedtraj ) :
    bfield_(bfield), seedtraj_(seedtraj)
  {
    config_.push_back(cfg);
    if(config().schedule().size() ==0)throw std::invalid_argument("Invalid configuration: no schedule");
  }

  // construct from configuration, reference (seed) fit, hits,and materials specific to this fit.
  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield, PTRAJ const& seedtraj,  HITCOL& hits, EXINGCOL& exings) : Track(cfg,bfield,seedtraj) {
    fit(hits,exings);
  }
  template <class KTRAJ> void Track<KTRAJ>::fit(HITCOL& hits, EXINGCOL& exings) {
    // set the seed time based on the min and max time from the inputs
    TimeRange refrange = getRange(hits,exings);
    seedtraj_.setRange(refrange);
    // if correcting for BField effects, define the domains
    DOMAINCOL domains;
    if(config().bfcorr_ ) createDomains(seedtraj_, refrange, domains, config().tol_);
    // Create the initial reference trajectory from the seed trajectory
    createTraj(seedtraj_,refrange,domains);
    // create all the other effects
//    effects_.reserve(hits.size()+exings.size()+domains.size());
    createEffects(hits,exings,domains);
    // now fit the track
    fit();
  }

  // extend an existing track
  template <class KTRAJ> void Track<KTRAJ>::extend(Config const& cfg, HITCOL& hits, EXINGCOL& exings) {
    // remember the previous config
    auto const& oldconfig = config_.back();
    // update the configuration
    config_.push_back(cfg);
    // configuation check
    if(config().schedule().size() ==0)throw std::invalid_argument("Invalid configuration: no schedule");
    // require the existing fit to be usable
    if(!fitStatus().usable())throw std::invalid_argument("Cannot extend unusable fit");
    // find the range of the added information, and extend as needed
    TimeRange exrange = fittraj_->range();
    if(hits.size() >0 || exings.size() > 0){
      TimeRange newrange = getRange(hits,exings);
      exrange.combine(newrange);
    }
    DOMAINCOL domains;
    if(config().bfcorr_ ) {
      // if the previous configuration didn't have domains, then create them for the full reference range
      if(!oldconfig.bfcorr_ || oldconfig.tol_ != config().tol_){
        // create domains for the whole range
        createDomains(*fittraj_, exrange, domains, config().tol_);
        // replace the domains.  This also replace the trajectory, as that must reference the new domains
        replaceDomains(domains);
        // tolerance has changed: remove the old domains and start again
      } else {
        // create domains just for the extensions
        TimeRange exlow(exrange.begin(),fittraj_->range().begin());
        if(exlow.range()>0.0) {
          DOMAINCOL lowdomains;
          createDomains(*fittraj_, exlow, lowdomains, config().tol_);
          domains.insert(lowdomains.begin(),lowdomains.end());
        }
        TimeRange exhigh(fittraj_->range().end(),exrange.end());
        if(exhigh.range()>0.0){
          DOMAINCOL highdomains;
          createDomains(*fittraj_, exhigh, highdomains, config().tol_);
          domains.insert(highdomains.begin(),highdomains.end());
        }
      }
    }
    // Extend the traj and create the effects for the new info and the new domains
    extendTraj(domains);
    createEffects(hits,exings,domains);
    // update all the effects for this new configuration
    for(auto& ieff : effects_ ) ieff->updateConfig(config());
    // now refit the track
    fit();
  }

  // replace domains when BField correction is added or changed. the traj must also be replaced, so that
  // the pieces correspond to the new domains
  template <class KTRAJ> void Track<KTRAJ>::replaceDomains(DOMAINCOL const& domains) {
    // if domains exist, clear them and remove all BField effects
    if(domains_.size() > 0){
      domains_.clear();
      // remove all existing BField effects
      auto ieff = effects_.begin();
      while(ieff != effects_.end()){
        const KKBFIELD* kkbf = dynamic_cast<const KKBFIELD*>(ieff->get());
        if(kkbf != 0){
          ieff = effects_.erase(ieff);
        } else {
          ++ieff;
        }
      }
    }
    // create new traj
    auto newtraj = std::make_unique<PTRAJ>();
    // loop over domains
    for(auto const& domain : domains) {
      double dtime = domain.begin();
      // Set the BField to the start of this domain
      auto bf = bfield_.fieldVect(fittraj_->position3(dtime));
      // loop until we're either out of this domain or the piece is out of this domain
      while(dtime < domain.end()){
        // find the nearest piece of the traj
        auto index = fittraj_->nearestIndex(dtime);
        auto const& oldpiece = *fittraj_->pieces()[index];
        // create a new piece
        KTRAJ newpiece(oldpiece,bf,dtime);
        // set the range as needed
        double endtime = (index < fittraj_->pieces().size()-1) ? std::min(domain.end(),oldpiece.range().end()) : domain.end();
        newpiece.range() = TimeRange(dtime,endtime);
        newtraj->append(newpiece);
        // update the time
        static double epsilon(1e-10);
        dtime = newpiece.range().end()+epsilon; // to avoid boundary
      }
    }
    // update all effects to reference this trajectory
    for (auto& eff : effects_) {
      auto const& ltrajptr = newtraj->nearestTraj(eff->time());
      eff->updateReference(ltrajptr);
    }
    // swap
    fittraj_.swap(newtraj);
  }

  template <class KTRAJ> void Track<KTRAJ>::extendTraj(DOMAINCOL const& domains ) {
    if(domains.size() > 0){
      TimeRange newrange(std::min(fittraj_->range().begin(),domains.begin()->begin()),
          std::max(fittraj_->range().end(),domains.rbegin()->end()));
      fittraj_->setRange(newrange);
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::createTraj(PTRAJ const& seedtraj , TimeRange const& range, DOMAINCOL const& domains ) {
    // if we're making local BField corrections, divide the trajectory into domain pieces.  Each will have equivalent parameters, but relative
    // to the local field
    if(config().bfcorr_ ) {
      if(fittraj_)throw std::invalid_argument("Initial reference trajectory must be empty");
      if(domains.size() == 0)throw std::invalid_argument("Empty domain collection");
      fittraj_ = std::make_unique<PTRAJ>();
      for(auto const& domain : domains) {
        // Set the BField to the start of this domain
        auto bf = bfield_.fieldVect(seedtraj.position3(domain.begin()));
        KTRAJ newpiece(seedtraj.nearestPiece(domain.begin()),bf,domain.begin());
        newpiece.range() = domain;
        fittraj_->append(newpiece);
      }
    } else {
      // use the middle of the range as the nominal BField for this fit:
      double tref = range.mid();
      VEC3 bf = bfield_.fieldVect(seedtraj.position3(tref));
      // create the first piece.  Note this constructor adjusts the parameters according to the local field
      KTRAJ firstpiece(seedtraj.nearestPiece(tref),bf,tref);
      firstpiece.range() = range;
      // create the piecewise trajectory from this
      fittraj_ = std::make_unique<PTRAJ>(firstpiece);
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::createEffects( HITCOL& hits, EXINGCOL& exings,DOMAINCOL const& domains ) {
    // pre-allocate as needed
//    effects_.reserve(effects_.size()+hits.size()+exings.size()+domains.size());
    // append the effects.  First, loop over the hits
    for(auto& hit : hits ) {
      // create the hit effects and insert them in the collection
      effects_.emplace_back(std::make_unique<KKMEAS>(hit));
      // update hit reference; this should be done on construction FIXME
      hit->updateReference(fittraj_->nearestTraj(hit->time()));
    }
    //add material effects
    for(auto& exing : exings) {
      effects_.emplace_back(std::make_unique<KKMAT>(exing,*fittraj_));
      // update xing reference; should be done on construction FIXME
      exing->updateReference(fittraj_->nearestTraj(exing->time()));
    }
    // add BField effects
    for( auto const& domain : domains) {
      // create the BField effect for integrated differences over this range
      effects_.emplace_back(std::make_unique<KKBFIELD>(config(),bfield_,domain));
    }
    // sort
//    std::sort(effects_.begin(),effects_.end(),KKEFFComp ());
    effects_.sort(KKEFFComp ());
    // store the inputs; these are just for and may not be in time order
    hits_.insert(hits_.end(),hits.begin(),hits.end());
    exings_.insert(exings_.end(),exings.begin(),exings.end());
    domains_.insert(domains.begin(),domains.end());
  }
  // fit the
  template <class KTRAJ> void Track<KTRAJ>::fit() {
    // execute the schedule of meta-iterations
    for(auto imiconfig=config().schedule().begin(); imiconfig != config().schedule().end(); imiconfig++){
      auto miconfig  = *imiconfig;
      // keep the meta-iteration count correct even if we extend the fit.
      unsigned nmeta = history_.size() == 0? 0 : fitStatus().miter_ + 1;
      unsigned niter(0);
      do{
        history_.push_back(Status(nmeta,niter++));
        // catch exceptions and record them in the status
        try {
          iterate(miconfig);
        } catch (std::exception const& error) {
          status().status_ = Status::failed;
          status().comment_ = error.what();
        }
      } while(canIterate());
      if(!status().usable())break;
    }
    // if the fit is usable, process the passive effects on either end
    if(config().ends_ && status().usable()){
      history_.push_back(fitStatus());
      status().comment_ = "EndProcessing";
      try {
        processEnds();
      } catch (std::exception const& error) {
        status().status_ = Status::failed;
        status().comment_ = status().comment_ + error.what();
      }
   }
    if(config().plevel_ > Config::none)print(std::cout, config().plevel_);
  }

  // single algebraic iteration
  template <class KTRAJ> void Track<KTRAJ>::iterate(MetaIterConfig const& miconfig) {
    if(config().plevel_ >= Config::basic)std::cout << "Processing fit iteration " << fitStatus().iter_ << std::endl;
    // update the effects for this configuration; this will sort the effects and find the iteration bounds
    bool first = status().iter_ == 0; // 1st iteration of a meta-iteration: update the effect internals
    for(auto& ieff : effects_ ) ieff->updateState(miconfig,first);
    // sort the sites, and set the iteration bounds
    effects_.sort(KKEFFComp ());
    KKEFFFWDBND fwdbnds;
    KKEFFREVBND revbnds;
    setBounds(fwdbnds,revbnds );
    FitStateArray states;
    TimeRange fitrange((*fwdbnds[0])->time(),(*revbnds[0])->time());
    initFitState(states, fitrange, config().dwt_/miconfig.varianceScale());
    // loop over relevant effects, adding their info to the fit state.  Also compute chisquared
    for(auto feff=fwdbnds[0];feff!=fwdbnds[1];++feff){
      auto effptr = feff->get();
      // update chisquared increment WRT the current state: only needed once
      Chisq dchisq = effptr->chisq(states[0].pData());
      status().chisq_ += dchisq;
      // process
      effptr->process(states[0],TimeDir::forwards);
      if(config().plevel_ >= Config::detailed && dchisq.nDOF() > 0){
        std::cout << "Chisq increment " << dchisq << " ";
        effptr->print(std::cout,config().plevel_-Config::detailed);
      }
    }
    if(status().chisq_.nDOF() >= (int)config().minndof_) { // I need a better way to define coverage as just having sufficien NDOF doesn't mean all parameters are constrained TODO
      double mintime(std::numeric_limits<double>::max());
      double maxtime(-std::numeric_limits<float>::max());
      for(auto beff = revbnds[0]; beff!=revbnds[1]; ++beff){
        auto effptr = beff->get();
        effptr->process(states[1],TimeDir::backwards);
        mintime = std::min(mintime,effptr->time());
        maxtime = std::max(maxtime,effptr->time());
      }
      // convert the fit result into a new trajectory
      // initialize the parameters to the backward processing end
      auto front = fittraj_->front();
      front.params() = states[1].pData();
      // extend range if needed
//      TimeRange maxrange(std::min(fittraj_->range().begin(),fwdbnds[0]->get()->time()),
//          std::max(fittraj_->range().end(),revbnds[0]->get()->time()));
      TimeRange maxrange(mintime-0.1,maxtime+0.1); //fixed time buffer shouldn't be needed FIXME
      front.setRange(maxrange);
      auto ptraj = std::make_unique<PTRAJ>(front);
      // process forwards, adding pieces as necessary.  This also sets the effects to reference the new trajectory
      for(auto& ieff=fwdbnds[0]; ieff != fwdbnds[1]; ++ieff) {
        ieff->get()->append(*ptraj,TimeDir::forwards);
      }
      setStatus(ptraj); // set the status for this iteration
      // prepare for the next iteration: first, update the references for effects outside the fit range
      // (the ones inside the range were updated above in 'append')
      if(status().usable()){
        for(auto feff=fwdbnds[1]; feff != effects_.end(); ++feff)
          feff->get()->updateReference(ptraj->nearestTraj(feff->get()->time()));
        for(auto beff=revbnds[1]; beff != effects_.rend(); ++beff)
          beff->get()->updateReference(ptraj->nearestTraj(beff->get()->time()));
      }
      // now all effects reference the new traj: we can swap it with the old
      fittraj_.swap(ptraj);
      if(config().plevel_ >= Config::complete)fittraj_->print(std::cout,1);
    } else {
      status().status_ = Status::lowNDOF;
    }
  }

  // initialize states used before iteration
  template <class KTRAJ> void Track<KTRAJ>::initFitState(FitStateArray& states, TimeRange const& fitrange, double dwt) {
    auto fwdtraj = fittraj_->nearestPiece(fitrange.begin());
    auto revtraj = fittraj_->nearestPiece(fitrange.end());
    // dweight the covariance, scaled by the temperature.
    fwdtraj.params().covariance() *= dwt;
    revtraj.params().covariance() *= dwt;
    auto fwdeff = Weights(fwdtraj.params());
    auto reveff = Weights(revtraj.params());
    states[0].append(fwdeff);
    states[1].append(reveff);
  }

  // finalize after iteration
  template <class KTRAJ> void Track<KTRAJ>::setStatus(PTRAJPTR& ptraj) {
    // compute parameter difference WRT previous iteration.  Compare at front and back ends
    auto const& ffront = ptraj->front();
    auto const& sfront = fittraj_->nearestPiece(ffront.range().mid());
    DVEC dpfront = ffront.params().parameters() - sfront.params().parameters();
    DMAT frontwt = sfront.params().covariance();
    if(! frontwt.Invert())throw std::runtime_error("Reference covariance uninvertible");
    double dpchisqfront = ROOT::Math::Similarity(dpfront,frontwt);
    // back
    auto const& fback = ptraj->back();
    auto const& sback = fittraj_->nearestPiece(fback.range().mid());
    DVEC dpback = fback.params().parameters() - sback.params().parameters();
    DMAT backwt = sback.params().covariance();
    if(! backwt.Invert())throw std::runtime_error("Reference covariance uninvertible");
    double dpchisqback = ROOT::Math::Similarity(dpback,backwt);
    // fit chisquared chang3
    double dchisq = config().convdchisq_ + 1e-4;  // initialize to insure 0th iteration doesn't converge
    if(fitStatus().iter_ > 0){
      auto prevstat = history_.rbegin();
      prevstat++;
      dchisq = fitStatus().chisq_.chisqPerNDOF() - prevstat->chisq_.chisqPerNDOF();
    }
    // check gap
    size_t igap;
    double maxgap,avggap;
    ptraj-> gaps(maxgap,igap,avggap);
    // test and update status
    if(avggap > config().divgap_ ) {
      status().status_ = Status::gapdiverged;
    } else if (dpchisqfront > config().pdchisq_ || dpchisqback > config().pdchisq_) {
      status().status_ = Status::paramsdiverged;
    } else if (dchisq > config().divdchisq_) {
      status().status_ = Status::chisqdiverged;
    } else if (status().chisq_.nDOF() < (int)config().minndof_){
      status().status_ = Status::lowNDOF;
    } else if(fabs(dchisq) < config().convdchisq_) {
      status().status_ = Status::converged;
    } else
      status().status_ = Status::unconverged;
  }

  // update between iterations
  template <class KTRAJ> void Track<KTRAJ>::setBounds(KKEFFFWDBND& fwdbnds, KKEFFREVBND& revbnds) {
// find bounds between first and last measurement
    for(auto ieff=effects_.begin();ieff!=effects_.end();++ieff){
      auto const* kkmeas = dynamic_cast<const KKMEAS*>(ieff->get());
      if(kkmeas != 0 && kkmeas->active()){
        fwdbnds[0] = ieff;
        revbnds[1] = KKEFFREV(ieff);
        break;
      }
    }
    for(auto ieff=effects_.rbegin();ieff!=effects_.rend();++ieff){
      auto const* kkmeas = dynamic_cast<const KKMEAS*>(ieff->get());
      if(kkmeas != 0 && kkmeas->active()){
        revbnds[0] = ieff;
        fwdbnds[1] = ieff.base();
        break;
      }
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::processEnds() {
    KKEFFFWDBND fwdbnds; // bounds for iterating
    KKEFFREVBND revbnds;
    setBounds(fwdbnds,revbnds);
    FitStateArray states;
    TimeRange fitrange((*fwdbnds[0])->time(),(*revbnds[0])->time());
    initFitState(states, fitrange, 1.0); // no deweighting
    for(auto feff=fwdbnds[1]; feff != effects_.end(); ++feff)
      feff->get()->process(states[1],TimeDir::forwards);
    for(auto reff=revbnds[1]; reff != effects_.rend(); ++reff)
      reff->get()->process(states[0],TimeDir::backwards);
    // finally, append the effects to the trajectory, using these states
    // skip any states that migrated to an unprocessed region
    for(auto feff=fwdbnds[1]; feff != effects_.end(); ++feff)
      if(feff->get()->time() > fittraj_->back().range().begin())feff->get()->append(*fittraj_,TimeDir::forwards);
    for(auto reff=revbnds[1]; reff != effects_.rend(); ++reff)
      if(reff->get()->time() < fittraj_->front().range().rbegin())reff->get()->append(*fittraj_,TimeDir::backwards);
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
      fitTraj().print(ost,detail-2);
    }
    if(detail > Config::complete) {
      ost << " Effects " << endl;
      for(auto const& eff : effects()) eff.get()->print(ost,detail-3);
    }
  }
  // divide a trajectory into magnetic 'domains' used to apply the BField corrections
  template<class KTRAJ> void Track<KTRAJ>::createDomains(PTRAJ const& ptraj, TimeRange const& range, DOMAINCOL& domains, double tol) const {
    double tstart = range.begin();
    do {
      // see how far we can go on the current traj before the BField change causes the momentum estimate to go out of tolerance
      // note this assumes the trajectory is accurate (geometric extrapolation only)
      auto const& ktraj = ptraj.nearestPiece(tstart);
      double tend = tstart + bfield_.rangeInTolerance(ktraj,tstart,tol);
      domains.insert(Domain(TimeRange(tstart,tend),tol));
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
      tmin = std::min(tmin,exing->time());
      tmax = std::max(tmax,exing->time());
    }
    return TimeRange(tmin,tmax);
  }

  template<class KTRAJ> template <class XTEST> bool Track<KTRAJ>::extrapolate(ExtraConfig const& xconfig, XTEST const& xtest) {
    bool retval(false);
    using enum TimeDir;
    if(this->fitStatus().usable()){
      if(config().bfcorr_){
        // iterate until the extrapolation condition is met
        double time = xconfig.xdir_ == forwards ? domains_.rbegin()->end() : domains_.begin()->begin();
        double tstart = time;
        while(xtest.needsExtrapolation(*this) && fabs(time-tstart) < xconfig.maxdt_){
          // create a domain for this extrapolation
          auto const& ktraj = fittraj_.nearestPiece(time);
          double dt = bfield_.rangeInTolerance(ktraj,time,xconfig.tol_); // always positive
          TimeRange range = xconfig.xdir_ == forwards ? TimeRange(time,time+dt) : TimeRange(time-dt,time);
          Domain domain(range,xconfig.tol_);
          addDomain(domain,xconfig.xdir_);
          time = xconfig.xdir_ == forwards ? domain.end() : domain.begin();
        }
        retval = true;
      }
    } else {
      retval = true;
    }
    return retval;
  }

  template<class KTRAJ> void Track<KTRAJ>::addDomain(Domain const& domain,TimeDir const& tdir) {
    using enum TimeDir;
    if(tdir == forwards){
      auto const& ktraj = fittraj_.nearestPiece(domains_.begin());
      FitState fstate(ktraj.params());
      effects_.emplace_back(std::make_unique<KKBFIELD>(config(),bfield_,domain));
      effects_.back()->process(fstate,tdir);
      effects_.back()->append(*fittraj_,tdir);
    } else {
      auto const& ktraj = fittraj_.nearestPiece(domains_.end());
      FitState fstate(ktraj.params());
      effects_.emplace_front(std::make_unique<KKBFIELD>(config(),bfield_,domain));
      effects_.front()->process(fstate,tdir);
      effects_.front()->append(*fittraj_,tdir);
    }
    domains_.insert(domain);
 }
}
#endif
