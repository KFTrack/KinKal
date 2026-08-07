#ifndef KinKal_Track_hh
#define KinKal_Track_hh
//
//  Primary class of the Kinematic Kalman fit.  This class owns the state describing
//  the fit inputs (measurements, material interactions, BField changes, etc), the result of the fit,
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
#include "KinKal/Fit/DomainWall.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/Status.hh"
#include "KinKal/Fit/Domain.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/CloneContext.hh"
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
#include <optional>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <ostream>

namespace KinKal {
  using DOMAINPTR = std::shared_ptr<Domain>;
  bool operator < (DOMAINPTR const& a, DOMAINPTR const&  b) { return a->begin() < b->begin(); }
  bool operator < (DOMAINPTR const& a, double val) { return a->begin() < val; }
  bool operator < (double val, DOMAINPTR const& b) { return val < b->begin(); }

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
      using KKDW = DomainWall<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PKTRAJPTR = std::unique_ptr<PKTRAJ>;
      using HIT = Hit<KTRAJ>;
      using HITPTR = std::shared_ptr<HIT>;
      using HITCOL = std::vector<HITPTR>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      using EXINGCOL = std::vector<EXINGPTR>;
      using DOMAINCOL = std::set<DOMAINPTR>;
      using CONFIGCOL = std::vector<Config>;
      using FitStateArray = std::array<FitState,2>;
      // construct from a set of hits and passive material crossings. This will compute the magnetic domains implicitly from the configuration
      Track(Config const& config, BFieldMap const& bfield, KTRAJ const& seedtraj, HITCOL& hits, EXINGCOL& exings );
      // reconstitute from a previous track
      Track(Config const& config, BFieldMap const& bfield, PKTRAJPTR& fittraj, HITCOL& hits, EXINGCOL& exings, DOMAINCOL& domains);
      // copy constructor
      Track(const Track& rhs, CloneContext& context);
      // extend an existing track with either new configuration, new hits, and/or new material xings
      void extend(Config const& config, HITCOL& hits, EXINGCOL& exings );
      // extrapolate the fit through the magnetic field with the given config until the given predicate is satisfied. This function requires
      // the fit be valid, otherwise the return code is false.  If successful the status, domains, and trajectory of the fit are updated
      // Note that the actual fit itself is unchanged
      template <class XTEST> bool extrapolate(TimeDir tdir, XTEST const& XTest);
      // extrapolate the fit through a material xing. Return value is true if the effect is added and fit updated.
      bool extrapolate(EXINGPTR const& exing,TimeDir const& tdir);
      // accessors
      std::vector<Status> const& history() const { return history_; }
      Status const& fitStatus() const { return history_.back(); } // most recent status
      PKTRAJ const& fitTraj() const { return *fittraj_; }
      bool hasTraj() const { return static_cast<bool>(fittraj_); } // false for a fit that failed before the traj was built
      KKEFFCOL const& effects() const { return effects_; }
      Config const& config() const { return config_.back(); }
      CONFIGCOL const& configs() const { return config_; }
      BFieldMap const& bfield() const { return bfield_; }
      HITCOL const& hits() const { return hits_; }
      EXINGCOL const& exings() const { return exings_; }
      DOMAINCOL const& domains() const { return domains_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
      TimeRange activeRange() const; // time range of active hits
      void extendTraj(TimeRange const& newrange);
      static TimeRange detectorRange(HITCOL& hits, EXINGCOL& exings,bool active=false);
    protected:
      Track(Config const& cfg, BFieldMap const& bfield);
      void fit(HITCOL& hits, EXINGCOL& exings, KTRAJ const& seedtraj);
      void fit(HITCOL& hits, EXINGCOL& exings, DOMAINCOL& domains, PKTRAJPTR& fittraj);
    private:
      void createEffects( HITCOL& hits, EXINGCOL& exings, DOMAINCOL const& domains);
      void convertSeed(KTRAJ const& seedtraj,TimeRange const& refrange, DOMAINCOL& domains);
      void fit(); // process the effects and create the trajectory.  This executes the current schedule
      bool createDomains(PKTRAJ const& ptraj, TimeRange const& range, DOMAINCOL& domains) const;
      // build one usable domain starting at tstart, or nothing if the field there can't support one
      std::optional<Domain> createDomain(PKTRAJ const& ptraj, double tstart, double tend) const;
      // input preconditions that don't depend on the BField. Checking them here keeps unusable input out
      // of the domain walk and out of the trajectory; records the reason and returns false on failure.
      bool validInput(TimeRange const& detrange);
      bool validInput(TimeRange const& detrange, KTRAJ const& seedtraj);
      bool setBounds(KKEFFFWDBND& fwdbnds,KKEFFREVBND& revbnds); // set the bounds.  Returns false if the bounds are empty
      bool extendDomains(TimeRange const& fitrange); // extend domains if the fit range changes.  Return value says if domains were added
      void updateDomains(PKTRAJ const& ptraj); // Update domains between iterations
      void iterate(MetaIterConfig const& miconfig);
      void setStatus(PKTRAJPTR& ptrajptr);
      void initFitState(FitStateArray& states, TimeRange const& fitrange, double dwt=1.0);
      PKTRAJPTR initTraj(FitState& state, TimeRange const& fitrange);
      bool canIterate() const;
      bool replaceDomains(DOMAINCOL const& domains);
      void extendTraj(DOMAINCOL const& domains);
      void processEnds();
      // add a single domain within the tolerance and extend the fit in the specified direction.
      void addDomain(Domain const& domain,TimeDir const& tdir,bool exact=false);
      auto& status() { return history_.back(); } // most recent status
      // divide a trajectory into magnetic 'domains' within which the BField can be considered constant (parameter change is within tolerance)
      // payload
      CONFIGCOL config_; // configuration
      BFieldMap const& bfield_; // magnetic field map
      std::vector<Status> history_; // fit status history; records the current iteration
      PKTRAJPTR fittraj_; // result of the current fit
      KKEFFCOL effects_; // effects used in this fit, sorted by time
      HITCOL hits_; // hits used in this fit
      EXINGCOL exings_; // material xings used in this fit
      DOMAINCOL domains_; // domains used in this fit, contiguous and sorted by time
  };
  // minimal constructor for subclasses. The resulting object has no fit.
  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield) :
    config_{cfg}, bfield_(bfield)
    {
      if(config().schedule().size() ==0)throw std::invalid_argument("Invalid configuration: no schedule");
      history_.push_back(Status(0,0,Status::unfit, "Construction"));
    }

  // construct from configuration, reference (seed) fit, hits,and materials specific to this fit. This will compute the domains according to the configuration before fitting.
  //
  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield, KTRAJ const& seedtraj, HITCOL& hits, EXINGCOL& exings)
    : Track(cfg,bfield)
  {
    fit(hits,exings, seedtraj);
  }

  template <class KTRAJ> Track<KTRAJ>::Track(Config const& cfg, BFieldMap const& bfield, PKTRAJPTR& fittraj, HITCOL& hits, EXINGCOL& exings, DOMAINCOL& domains)
    : Track(cfg,bfield)
  {
    fit(hits,exings,domains,fittraj);
  }

  template <class KTRAJ> void Track<KTRAJ>::fit(HITCOL& hits, EXINGCOL& exings, KTRAJ const& seedtraj) {
    auto detrange = detectorRange(hits,exings,true);
    if(!validInput(detrange,seedtraj)) return;
    // convert the seed traj to a piecewaise traj. This creates the domains
    DOMAINCOL domains;
    convertSeed(seedtraj,detrange,domains);
    if(status().needsFit()){
      // convert all the primary info to fit effects
      createEffects(hits,  exings, domains);
      // now fit
      fit();
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::fit(HITCOL& hits, EXINGCOL& exings, DOMAINCOL& domains, PKTRAJPTR& fittraj) {
    fittraj_ = std::move(fittraj); // steal the underlying object
    // truncate the domains and fit trajectory to be within the detector range
    auto detrange = detectorRange(hits,exings,true);
    if(!validInput(detrange)) return;
    if(domains.size() > 0){
      auto idom = domains.begin();
      // stop at the 1st domain overlaping the detector range, and erase all elements up to that point
      while(idom != domains.end() && !(detrange.overlaps((*idom)->range())))++idom;
      if(idom != domains.begin())domains.erase(domains.begin(),--idom);// leave the overlapping piece
      auto jdom= domains.rbegin();
      while(jdom != domains.rend() && !(detrange.overlaps((*jdom)->range())))++jdom;
      domains.erase(jdom.base(),domains.end()); // base points 1 past the reverse iterator
      // hit times can fall outside every saved domain when the domain span is shorter than the
      // trajectory piece; soft-fail rather than dereference an empty set
      if(domains.empty()){
        history_.emplace_back(0,0,Status::outsidemap, "Empty domains after detector-range trim");
        return;
      }
      // trim the trajectory to this range
      detrange.combine((*domains.begin())->range());
      detrange.combine((*domains.rbegin())->range());
    }
    fittraj_->setRange(detrange,true);
    createEffects(hits,exings,domains);
    fit();
  }

  // copy constructor
  template<class KTRAJ> Track<KTRAJ>::Track(const Track& rhs, CloneContext& context) :
    config_(rhs.configs()),
    bfield_(rhs.bfield()),
    history_(rhs.history())
  {
    fittraj_ = std::make_unique<PKTRAJ>(rhs.fitTraj());
    hits_.reserve(rhs.hits().size());
    for (const auto& ptr: rhs.hits()){
      auto hit = ptr->clone(context);
      hits_.push_back(hit);
    }
    exings_.reserve(rhs.exings().size());
    for (const auto& ptr: rhs.exings()){
      auto xng = ptr->clone(context);
      exings_.push_back(xng);
    }
    for (const auto& ptr: rhs.domains()){
      auto dmn = context.get(ptr);
      domains_.insert(dmn);
    }
    for (const auto& ptr: rhs.effects()){
      auto eff = ptr->clone(context);
      effects_.push_back(std::move(eff));
    }
  };

  // extend an existing track
  template <class KTRAJ> void Track<KTRAJ>::extend(Config const& cfg, HITCOL& hits, EXINGCOL& exings) {
    // require the existing fit to be usable
    if(!fitStatus().usable())return;
    // configuation check
    if(cfg.schedule().size() ==0)throw std::invalid_argument("Invalid configuration: no schedule");
    // save the new configuration
    config_.push_back(cfg);
    // remember the previous config
    auto oldciter = config_.rbegin(); oldciter++;
    auto const& oldconfig = *oldciter;
    // find the range of the added information, and extend as needed
    TimeRange exrange = fittraj_->range();
    if(hits.size() >0 || exings.size() > 0){
      TimeRange newrange = detectorRange(hits,exings);
      exrange.combine(newrange);
    }
    DOMAINCOL domains;
    bool dok(true);
    if(config().bfcorr_ ) {
      // if the previous configuration didn't have domains, then create them for the full reference range
      if(!oldconfig.bfcorr_ || oldconfig.tol_ != config().tol_){
        // create domains for the whole range
        dok &= createDomains(*fittraj_,exrange, domains);
        // replace previous  domains with these.  This replaces the trajectory and bfield-related effects
        if(bfield_.protecting() && dok && domains.empty()){
          // Map-edge stop before any domain: do not call replaceDomains on an empty set.
          dok = false;
        } else if(dok){
          // a tighter extension tolerance can flip omega near a collapsing field, which the
          // parameterization reads as a charge change; leave the track untouched and record why
          if(!replaceDomains(domains)){
            history_.push_back(Status(0));
            status().status_ = Status::incompatiblepiece;
            status().comment_ = std::string("Domain replacement: incompatible piece");
            return;
          }
        }
      } else {
        // create domains just for the extensions
        TimeRange exlow(exrange.begin(),fittraj_->range().begin());
        if(exlow.range()>0.0) {
          DOMAINCOL lowdomains;
          dok &= createDomains(*fittraj_,exlow, lowdomains);
          if(dok)domains.insert(lowdomains.begin(),lowdomains.end());
        }
        TimeRange exhigh(fittraj_->range().end(),exrange.end());
        if(exhigh.range()>0.0){
          DOMAINCOL highdomains;
          dok &= createDomains(*fittraj_,exhigh, highdomains);
          if(dok)domains.insert(highdomains.begin(),highdomains.end());
        }
      }
    }
    if(!dok){
      // domain calculation failed. Under low-field protection, keep a previously usable fit
      // (map-edge truncation is preferred inside createDomains; this is a safety net).
      if(bfield_.protecting() && fitStatus().usable()) return;
      history_.push_back(Status(0));
      status().status_ = Status::outsidemap;
      status().comment_ = std::string("Extension error");
      return;
    }
    // Extend the traj and create the effects for the new info and the new domains
    extendTraj(domains);
    createEffects(hits,exings,domains);
    // update all the effects for this new configuration
    for(auto& ieff : effects_ ) ieff->updateConfig(config());
    // now refit the track
    fit();
  }

  // replace domains when DomainWall correction is added or changed. the traj must also be replaced, so that
  // the pieces correspond to the new domains. The new traj is geometrically equivalent, but not parametrically equal.
  // Two-phase: build the replacement without touching any member, then commit. Returns false, leaving
  // the track untouched, if a transformed piece cannot describe the same particle as those before it.
  template <class KTRAJ> bool Track<KTRAJ>::replaceDomains(DOMAINCOL const& domains) {
    if(domains.size() == 0) return false;
    TimeRange drange(domains.begin()->get()->begin(),domains.rbegin()->get()->end());
    TimeRange oldrange = fittraj_->range();
    auto newtraj = std::make_unique<PKTRAJ>();
    // loop over domains, splitting the overlapping traj pieces at the domain walls, and transforming them to reference the domain's field
    // This increases the number of traj pieces.
    for(auto const& domain : domains) {
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using DKTRAJ = std::deque<KTRAJPTR>;
      using DKTRAJCITER = typename DKTRAJ::const_iterator;
      // clamp the lookup to the existing traj: domains may extend past it, and extending fittraj_ to
      // reach them is the pre-mutation that used to need rollback. drange is restored on newtraj below.
      double tlo = std::max(domain->begin(),oldrange.begin());
      double thi = std::min(domain->end(),oldrange.end());
      if(tlo < thi){
        // find the range of existing ptraj pieces that overlaps with this domain's range
        DKTRAJCITER first,last;
        fittraj_->pieceRange(TimeRange(tlo,thi),first,last);
        // loop over these pieces; first and last can be the same!
        auto olditer = first;
        do {
          auto const& oldpiece = **olditer;
          // copy this piece, translating bnom to this domain's field
          KTRAJ newpiece(oldpiece,domain->bnom(),domain->range().mid());
          // set the range for this piece, making sure it is non-zero
          double tstart = std::max(domain->begin(), oldpiece.range().begin());
          double tend = std::min(domain->end(),oldpiece.range().end());
          if(tstart < tend){
            // test only where the old code would actually have appended (and thrown)
            if(!newtraj->compatible(newpiece)) return false;
            newpiece.range() = TimeRange(tstart,tend);
            newtraj->append(newpiece);
          }
          if(olditer != last)++olditer;
        } while(olditer != last);
      } else {
        // domain lies entirely outside the current traj: transform the nearest piece to cover it
        auto const& oldpiece = fittraj_->nearestPiece(domain->range().mid());
        KTRAJ newpiece(oldpiece,domain->bnom(),domain->range().mid());
        if(!newtraj->compatible(newpiece)) return false;
        newpiece.range() = domain->range();
        newtraj->append(newpiece);
      }
    }
    if(newtraj->pieces().size() == 0) return false;
    // restore the full domain span on the replacement, matching the range the old code reached by
    // extending fittraj_ up front
    newtraj->setRange(drange);
    // commit: clear old domains / DomainWall effects, retarget remaining effects, swap traj
    if(domains_.size() > 0){
      domains_.clear();
      auto ieff = effects_.begin();
      while(ieff != effects_.end()){
        const KKDW* kkbf = dynamic_cast<const KKDW*>(ieff->get());
        if(kkbf != 0){
          ieff = effects_.erase(ieff);
        } else {
          ++ieff;
        }
      }
    }
    for (auto& eff : effects_) {
      eff->updateReference(*newtraj);
    }
    fittraj_.swap(newtraj);
    return true;
  }

  template <class KTRAJ> void Track<KTRAJ>::extendTraj(DOMAINCOL const& domains ) {
    if(domains.size() > 0){
      TimeRange newrange(std::min(fittraj_->range().begin(),domains.begin()->get()->begin()),
          std::max(fittraj_->range().end(),domains.rbegin()->get()->end()));
      fittraj_->setRange(newrange);
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::extendTraj(TimeRange const& newrange){
    TimeRange temprange(std::min(fittraj_->range().begin(),newrange.begin()),
        std::max(fittraj_->range().end(),newrange.end()));
    fittraj_->setRange(temprange);
  }

  // no active hits or material crossings: there is nothing to fit, and a null range would otherwise be
  // walked for domains and then set on the trajectory
  template <class KTRAJ> bool Track<KTRAJ>::validInput(TimeRange const& detrange) {
    if(detrange.null()){
      history_.emplace_back(0,0,Status::lowNDOF, "No active hits or material crossings");
      return false;
    }
    return true;
  }

  // as above, plus the seed itself must be finite: a NaN seed otherwise propagates into the
  // parameterization and is only caught much later, if at all
  template <class KTRAJ> bool Track<KTRAJ>::validInput(TimeRange const& detrange, KTRAJ const& seedtraj) {
    if(!validInput(detrange)) return false;
    double tref = detrange.mid();
    auto spos = seedtraj.position3(tref);
    if(!std::isfinite(seedtraj.momentum(tref)) || !std::isfinite(spos.X()) ||
        !std::isfinite(spos.Y()) || !std::isfinite(spos.Z())){
      history_.emplace_back(0,0,Status::failed, "Non-finite seed trajectory");
      return false;
    }
    return true;
  }

  template <class KTRAJ> void Track<KTRAJ>::convertSeed(KTRAJ const& seedtraj,TimeRange const& range, DOMAINCOL& domains) {
    // if we're making local DomainWall corrections, divide the trajectory into domain pieces.  Each will have equivalent parameters, but relative
    // to the local field
    if(config().bfcorr_ ) {
      // convert the seedtraj to a PKTRAJ
      PKTRAJ straj(seedtraj);
      bool dok = createDomains(straj, range, domains);
      if(!dok){
        // failed initialization
        history_.emplace_back(0,0,Status::outsidemap, "Domain initialization error");
        return;
      }
      // create the initial fittraj from the seed. Each piece will have 'the same' physical trajectory as the seed, but reference the local domain BField
      fittraj_ = std::make_unique<PKTRAJ>();
      for(auto const& domain : domains) {
        // Set the DomainWall to the start of this domain
        auto bf = bfield_.fieldVect(seedtraj.position3(domain->mid()));
        KTRAJ newpiece(seedtraj,bf,domain->mid());
        newpiece.range() = domain->range();
        // same routine incompatibility as replaceDomains; this append was previously unguarded, so a
        // CentralHelix omega sign flip here threw std::invalid_argument out of the Track constructor
        if(!fittraj_->compatible(newpiece)){
          history_.emplace_back(0,0,Status::incompatiblepiece, "Seed conversion: incompatible piece");
          return;
        }
        fittraj_->append(newpiece);
      }
      // zero domains leaves fittraj_ empty, which createEffects would report by throwing
      // std::length_error out of the constructor; soft-fail instead
      if(fittraj_->pieces().empty()){
        history_.emplace_back(0,0,Status::outsidemap, "Empty seed trajectory (no domains)");
        return;
      }
    } else {
      // use the middle of the range as the nominal BField for this fit:
      double tref = range.mid();
      VEC3 bf = bfield_.fieldVect(seedtraj.position3(tref));
      // create the first piece.  Note this constructor adjusts the parameters according to the local field
      KTRAJ firstpiece(seedtraj,bf,tref);
      firstpiece.range() = range;
      // create the piecewise trajectory from this
      fittraj_ = std::make_unique<PKTRAJ>(firstpiece);
    }
  }

  template <class KTRAJ> void Track<KTRAJ>::createEffects( HITCOL& hits, EXINGCOL& exings,DOMAINCOL const& domains ) {
    // create and append the effects.  First, loop over the hits
    for(auto& hit : hits ) {
      effects_.emplace_back(std::make_unique<KKMEAS>(hit,*fittraj_));
    }
    //add material effects
    for(auto& exing : exings) {
      effects_.emplace_back(std::make_unique<KKMAT>(exing,*fittraj_));
    }
    // add DomainWall effects
    if(config().bfcorr_ && domains.size() > 1){
      auto nextdom = domains.cbegin();
      auto prevdom = nextdom;
      ++nextdom;
      while( nextdom != domains.cend() ){
        // must be contiguous
        if(fabs(prevdom->get()->end()-nextdom->get()->begin())>1e-10)throw std::invalid_argument("Invalid domains");
        effects_.emplace_back(std::make_unique<KKDW>(*prevdom,*nextdom ,*fittraj_));
        prevdom = nextdom;
        ++nextdom;
      }
    }
    // store the inputs; these are just for and may not be in time order
    hits_.insert(hits_.end(),hits.begin(),hits.end());
    exings_.insert(exings_.end(),exings.begin(),exings.end());
    domains_.insert(domains.cbegin(),domains.cend());
  }

  // fit the track
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
    // if the fit is usable, extend the track as necessary
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

  // single algebraic iteration of the fit
  template <class KTRAJ> void Track<KTRAJ>::iterate(MetaIterConfig const& miconfig) {
    if(config().plevel_ >= Config::basic)std::cout << "Processing fit iteration " << fitStatus().iter_ << std::endl;
    // update the effects for this configuration; this will sort the effects and find the iteration bounds
    bool first = status().iter_ == 0; // 1st iteration of a meta-iteration: update the effect internals
    for(auto& ieff : effects_ ) ieff->updateState(miconfig,first);
    // sort the sites, and set the iteration bounds
    effects_.sort(KKEFFComp());
    KKEFFFWDBND fwdbnds; // bounding sites used for fitting
    KKEFFREVBND revbnds;
    if(!setBounds(fwdbnds,revbnds)){
      status().status_ = Status::lowNDOF;
      return;
    }
    // update the domains
    TimeRange fitrange(fwdbnds[0]->get()->time(),revbnds[0]->get()->time());
    if(fitrange.range() == 0){
      status().status_ = Status::lowNDOF;
      return;
    }
    if(config().bfcorr_){
      // update the limits if new DW effects were added
      if(extendDomains(fitrange))setBounds(fwdbnds,revbnds);
    }
    FitStateArray states;
    initFitState(states, fitrange, config().dwt_/miconfig.varianceScale());
    // process the relevant effects forwards in time, adding their info to the fit state.  Also compute chisquared
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
    if(status().chisq_.nDOF() < (int)config().minndof_) { // I need a better way to define coverage as this test doesn't guarantee all parameters are constrained TODO
      status().status_ = Status::lowNDOF;
      return;
    }
    // Now process effects backwards in time
    for(auto beff = revbnds[0]; beff!=revbnds[1]; ++beff){
      auto effptr = beff->get();
      effptr->process(states[1],TimeDir::backwards);
    }
    // convert the innermost state into a new trajectory
    auto ptraj = initTraj(states[1],fitrange);
    // append forwards in time
    for(auto& ieff=fwdbnds[0]; ieff != fwdbnds[1]; ++ieff) {
      ieff->get()->append(*ptraj,TimeDir::forwards);
    }
    setStatus(ptraj); // set the status for this iteration.  This compares the trajs, so it must be done before swapping out the old fit
    if(status().usable()){
      // prepare for the next iteration.  Update the domains
      updateDomains(*ptraj);
      // update the effects to reference the new trajectory
      for(auto& effptr : effects_) effptr->updateReference(*ptraj);
      fittraj_.swap(ptraj);
      // now all effects reference the new traj: we can swap the fit to that.
      if(config().plevel_ >= Config::complete)fittraj_->print(std::cout,1);
    }
  }

  template <class KTRAJ> std::unique_ptr<ParticleTrajectory<KTRAJ>> Track<KTRAJ>::initTraj(FitState& state, TimeRange const& fitrange) {
    // use the existing fit to define the parameterization
    auto front = fittraj_->front();
    // initialize the parameters to the earliest backward processing end
    front.params() = state.pData();
    // set bnom for the front piece parameters to the domain used
    if(config().bfcorr_){
      // find the relevant domain
      double ftime = fitrange.begin();
      for(auto const& domain : domains_) {
        if(domain->range().inRange(ftime)){
          front.resetBNom(domain->bnom());
          break;
        }
      }
    }
    // set the range
    front.setRange(fitrange);
    return std::make_unique<PKTRAJ>(front);
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
  template <class KTRAJ> void Track<KTRAJ>::setStatus(PKTRAJPTR& ptraj) {
    // compute parameter difference WRT previous iteration.  Compare at front and back ends
    auto const& ffront = ptraj->front();
    // test covariance
    auto fdiag = ffront.params().covariance().Diagonal();
    for(size_t ipar = 0; ipar < NParams(); ++ipar){
      if(fdiag[ipar] < 0.0 || std::isnan(fdiag[ipar])){
        status().status_ = Status::failed;
        return;
      }
    }
    auto const& sfront = fittraj_->nearestPiece(ffront.range().mid());
    DVEC dpfront = ffront.params().parameters() - sfront.params().parameters();
    DMAT frontwt = sfront.params().covariance();
    if(! frontwt.Invert())throw std::runtime_error("Reference covariance uninvertible");
    double dpchisqfront = ROOT::Math::Similarity(dpfront,frontwt);
    // back
    auto const& fback = ptraj->back();
    auto bdiag = fback.params().covariance().Diagonal();
    for(size_t ipar = 0; ipar < NParams(); ++ipar){
      if(bdiag[ipar] < 0.0 || std::isnan(bdiag[ipar])){
        status().status_ = Status::failed;
        return;
      }
    }
    auto const& sback = fittraj_->nearestPiece(fback.range().mid());
    DVEC dpback = fback.params().parameters() - sback.params().parameters();
    DMAT backwt = sback.params().covariance();
    if(! backwt.Invert())throw std::runtime_error("Reference covariance uninvertible");
    double dpchisqback = ROOT::Math::Similarity(dpback,backwt);
    // fit chisquared change
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

  template <class KTRAJ> bool Track<KTRAJ>::setBounds(KKEFFFWDBND& fwdbnds,KKEFFREVBND& revbnds) {
    // find the bounding sites for algebraic processing this iteraiton.  This excludes effects which can't change the fit
    // find first measurement
    bool retval(false);
    fwdbnds[0] = effects_.end();
    for(auto ieff=effects_.begin();ieff!=effects_.end();++ieff){
      auto const* kkmeas = dynamic_cast<const KKMEAS*>(ieff->get());
      if(kkmeas != 0 && kkmeas->active()){
        fwdbnds[0] = ieff;
        revbnds[1] = KKEFFREV(ieff);
        break;
      }
    }
    if(fwdbnds[0] != effects_.end()){
      retval = true;
      // now the last measurement
      for(auto ieff=effects_.rbegin();ieff!=effects_.rend();++ieff){
        auto const* kkmeas = dynamic_cast<const KKMEAS*>(ieff->get());
        if(kkmeas != 0 && kkmeas->active()){
          revbnds[0] = ieff;
          fwdbnds[1] = ieff.base();
          break;
        }
      }
    }
    return retval;
  }

  template <class KTRAJ> void Track<KTRAJ>::updateDomains(PKTRAJ const& ptraj) {
    for(auto& domain : domains_) {
      domain->updateBNom(bfield_.fieldVect(ptraj.position3(domain->mid())));
    }
  }

  template <class KTRAJ> bool Track<KTRAJ>::extendDomains(TimeRange const& fitrange) {
    bool retval(false);
    // then, check if the range has
    TimeRange drange(domains().begin()->get()->begin(),domains().rbegin()->get()->end());
    if(!drange.contains(fitrange)){
      retval = true;
      // we need to extend the domains.  First backwards
      if(drange.begin() > fitrange.begin()){
        double time = drange.begin();
        while(time > fitrange.begin()){
          auto const& ktraj = fittraj_->nearestPiece(time);
          double dt = std::max(bfield_.rangeInTolerance(ktraj,time,config().tol_),config().mindtstep_);
          // clamp the domain low bound to the active range minus domainmargin_ (max = unclamped)
          double dlo = std::max(time-dt, fitrange.begin() - config().domainmargin_);
          TimeRange range(dlo,time);
          // sample BNom at the domain-midpoint piece when confined (domainmargin_ set), else at the nearest piece
          auto const& straj = (config().domainmargin_ < std::numeric_limits<double>::max()) ? fittraj_->nearestPiece(range.mid()) : ktraj;
          Domain domain(range,bfield_.fieldVect(straj.position3(range.mid())));
          addDomain(domain,TimeDir::backwards);
          time = domain.begin();
        }
      }
      // then forwards
      if(drange.end() < fitrange.end()){
        double time = drange.end();
        while(time < fitrange.end()){
          auto const& ktraj = fittraj_->nearestPiece(time);
          double dt = std::max(bfield_.rangeInTolerance(ktraj,time,config().tol_),config().mindtstep_);
          // clamp the domain high bound to the active range plus domainmargin_ (max = unclamped)
          double dhi = std::min(time+dt, fitrange.end() + config().domainmargin_);
          TimeRange range(time,dhi);
          // sample BNom at the domain-midpoint piece when confined (domainmargin_ set), else at the nearest piece
          auto const& straj = (config().domainmargin_ < std::numeric_limits<double>::max()) ? fittraj_->nearestPiece(range.mid()) : ktraj;
          Domain domain(range,bfield_.fieldVect(straj.position3(range.mid())));
          addDomain(domain,TimeDir::forwards);
          time = domain.end();
        }
      }
    }
    return retval;
  }

  template <class KTRAJ> void Track<KTRAJ>::processEnds() {
    // sort effects in case ends have migrated
    effects_.sort(KKEFFComp());
    // update domains as needed to cover the end effects
    if(config().bfcorr_){
      TimeRange endrange(effects_.front()->time(),effects_.back()->time());
      extendDomains(endrange);
    }
    KKEFFFWDBND fwdbnds; // bounding sites used for fitting
    KKEFFREVBND revbnds;
    if(!setBounds(fwdbnds,revbnds)){
      status().status_ = Status::lowNDOF;
      return;
    }
    // initialize the fit state where we left off processing
    FitStateArray states;
    TimeRange fitrange(fwdbnds[0]->get()->time(),revbnds[0]->get()->time());
    initFitState(states, fitrange, 1.0); // no deweighting
    // process forwards and backwards
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
    // convertSeed returns before fittraj_ is built when domain initialization fails, so a soft-failed
    // fit has no trajectory to print; dereferencing it here segfaulted
    if(!hasTraj()){
      ost << "(no trajectory: " << fitStatus().comment_ << ")" << endl;
      return;
    }
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
  // build one domain starting at tstart, clipped to end no later than tend. Returns nothing when the
  // field can't support a domain here; the map decides that, this just reports it.
  template<class KTRAJ> std::optional<Domain> Track<KTRAJ>::createDomain(PKTRAJ const& ptraj, double tstart, double tend) const {
    auto const& ktraj = ptraj.nearestPiece(tstart);
    if(!bfield_.usable(ktraj.position3(tstart))) return std::nullopt;
    double trange = bfield_.domainStep(ktraj,tstart,config().tol_,config().mindtstep_);
    double dhi = std::min(tstart+trange,tend);
    if(dhi <= tstart) return std::nullopt;
    TimeRange drange(tstart,dhi);
    // the domain carries the field sampled at its midpoint. Requiring that sample to be usable is what
    // protection buys; unprotected, fieldVect just reports null outside the map, as it always did.
    // Sample the midpoint on the midpoint's own piece only when confined (domainmargin_ set); the
    // unconfined walk samples it on the start piece, as extendDomains does.
    bool confined = config().domainmargin_ < std::numeric_limits<double>::max();
    auto const& straj = confined ? ptraj.nearestPiece(drange.mid()) : ktraj;
    VEC3 midpos = straj.position3(drange.mid());
    if(bfield_.protecting() && !bfield_.usable(midpos)) return std::nullopt;
    return Domain(drange,bfield_.fieldVect(midpos));
  }

  // divide a trajectory into magnetic 'domains' used to apply the DomainWall corrections
  template<class KTRAJ> bool Track<KTRAJ>::createDomains(PKTRAJ const& ptraj, TimeRange const& range, DOMAINCOL& domains) const {
    if(!config().bfcorr_) return true;
    // No usable domain means: soft stop when protection is on (keep what we have, let the caller
    // proceed), otherwise the failure the map used to signal by throwing out of fieldDeriv.
    if(config().domainmargin_ < std::numeric_limits<double>::max()){
      // confined (domainmargin_ set): walk only within the active range +/- margin, clipping the last domain
      double const tlo = range.begin() - config().domainmargin_;
      double const thi = range.end() + config().domainmargin_;
      double tstart = tlo;
      while(tstart < thi){
        auto domain = createDomain(ptraj,tstart,thi);
        if(!domain) return bfield_.protecting();
        domains.emplace(std::make_shared<Domain>(*domain));
        tstart = domain->end();
      }
    } else {
      // Unconfined (default): half-domain overhang so the first/last effect sits mid-domain. The loop
      // bound tracks the current step, as upstream, so domains are unclipped and full length.
      auto const& ktraj0 = ptraj.nearestPiece(range.begin());
      if(!bfield_.usable(ktraj0.position3(range.begin()))) return bfield_.protecting();
      double trange = bfield_.domainStep(ktraj0,range.begin(),config().tol_,config().mindtstep_);
      double tstart = range.begin() - 0.5*trange;
      do {
        auto domain = createDomain(ptraj,tstart,std::numeric_limits<double>::max());
        if(!domain) return bfield_.protecting();
        trange = domain->range().range();
        domains.emplace(std::make_shared<Domain>(*domain));
        tstart = domain->end();
      } while(tstart < range.end() + 0.5*trange);
    }
    return true;
  }

  template<class KTRAJ> TimeRange Track<KTRAJ>::detectorRange(HITCOL& hits, EXINGCOL& exings,bool active) {
    double tmin = std::numeric_limits<double>::max();
    double tmax = -std::numeric_limits<double>::max();
    // can't assume effects are sorted
    for(auto const& hit : hits){
      if((!active) || hit->active()){
        tmin = std::min(tmin,hit->time());
        tmax = std::max(tmax,hit->time());
      }
    }
    for(auto const& exing : exings){
      if((!active) || exing->active() ){
        tmin = std::min(tmin,exing->time());
        tmax = std::max(tmax,exing->time());
      }
    }
    // no (active) effects leaves tmin>tmax; return a null range instead of an invalid one that would throw
    if(tmax < tmin) return TimeRange();
    return TimeRange(tmin,tmax);
  }

  template<class KTRAJ> template <class XTEST> bool Track<KTRAJ>::extrapolate(TimeDir tdir, XTEST const& xtest) {
    bool retval = fitStatus().usable();
    if(retval){
      if(config().bfcorr_){
        // opt-in low-field handoff, for extrapolation that must leave the field map; with minfield_ 0
        // the unprotected domain walk below is used unchanged
        if(bfield_.protecting()){
          auto geometricExtend = [&](double tmax_remaining) {
            if(tmax_remaining <= 0.0) return;
            auto& endpiece = tdir == TimeDir::forwards ? fittraj_->backPtr() : fittraj_->frontPtr();
            double time = tdir == TimeDir::forwards ? endpiece->range().end() : endpiece->range().begin();
            double tstart = time;
            bool needsext(true);
            do {
              TimeRange newrange = tdir == TimeDir::forwards ?
                TimeRange(endpiece->range().begin(),endpiece->range().end()+xtest.maxDtStep())
                :
                TimeRange(endpiece->range().begin()-xtest.maxDtStep(),endpiece->range().end());
              endpiece->setRange(newrange);
              time = tdir == TimeDir::forwards ? endpiece->range().end() : endpiece->range().begin();
              needsext = xtest.needsExtrapolation(*fittraj_,tdir);
            } while(needsext && fabs(time-tstart) < tmax_remaining);
          };

          bool handed_off = false;
          double time = tdir == TimeDir::forwards ? domains_.crbegin()->get()->end() : domains_.cbegin()->get()->begin();
          double tstart = time;
          try {
            while(fabs(time-tstart) < xtest.maxDt() && xtest.needsExtrapolation(*fittraj_,tdir) ){
              auto const& ktraj = fittraj_->nearestPiece(time);
              if( !std::isfinite(ktraj.momentum(time)) ) break;

              // the map decides whether this point can carry field-corrected transport; asking it first
              // also keeps rangeInTolerance from sampling fieldDeriv out of range
              if(!bfield_.usable(ktraj.position3(time))){
                // leave the bfcorr / DomainWall path; free-particle continuation is a geometric
                // range-extend of the current end piece, with no parameter rebuild at B ~ 0
                handed_off = true;
                break;
              }

              double dt = std::clamp(bfield_.rangeInTolerance(ktraj,time,xtest.dpTolerance()),config().mindtstep_,xtest.maxDtStep());
              TimeRange range = tdir == TimeDir::forwards ? TimeRange(time,time+dt) : TimeRange(time-dt,time);
              VEC3 midpos = ktraj.position3(range.mid());
              if(!bfield_.usable(midpos)){
                handed_off = true;
                break;
              }
              Domain domain(range,bfield_.fieldVect(midpos));
              addDomain(domain,tdir,true);
              time = tdir == TimeDir::forwards ? domain.end() : domain.begin();
            }
          } catch (std::exception const& error) {
            history_.push_back(Status(0));
            status().status_ = Status::outsidemap;
            status().comment_ = std::string("Extrapolation error");
            retval = false;
          }
          if(retval && handed_off && xtest.needsExtrapolation(*fittraj_,tdir)){
            geometricExtend(xtest.maxDt() - fabs(time-tstart));
          }
        } else {
          // no low-field protection: the unprotected bfcorr domain walk, unchanged
          try {
            double time = tdir == TimeDir::forwards ? domains_.crbegin()->get()->end() : domains_.cbegin()->get()->begin();
            double tstart = time;
            while(fabs(time-tstart) < xtest.maxDt() && xtest.needsExtrapolation(*fittraj_,tdir) ){
              auto const& ktraj = fittraj_->nearestPiece(time);
              double dt = std::clamp(bfield_.rangeInTolerance(ktraj,time,xtest.dpTolerance()),config().mindtstep_,xtest.maxDtStep());
              TimeRange range = tdir == TimeDir::forwards ? TimeRange(time,time+dt) : TimeRange(time-dt,time);
              auto domainfield = bfield_.fieldVect(ktraj.position3(range.mid()));
              Domain domain(range,domainfield);
              addDomain(domain,tdir,true);
              time = tdir == TimeDir::forwards ? domain.end() : domain.begin();
            }
          } catch (std::exception const& error) {
            history_.push_back(Status(0));
            status().status_ = Status::outsidemap;
            status().comment_ = std::string("Extrapolation error");
            retval = false;
          }
          retval = true;
        }
      } else {
        // geometric extrapolation of the end piece; no need to protect
        auto& endpiece = tdir == TimeDir::forwards ? fittraj_->backPtr() : fittraj_->frontPtr();
        double time = tdir == TimeDir::forwards ? endpiece->range().end() : endpiece->range().begin();
        double tstart = time;
        bool needsext(true);
        do {
          TimeRange newrange = tdir == TimeDir::forwards ?
            TimeRange(endpiece->range().begin(),endpiece->range().end()+xtest.maxDtStep())
            :
            TimeRange(endpiece->range().begin()-xtest.maxDtStep(),endpiece->range().end());
          endpiece->setRange(newrange);
          time = tdir == TimeDir::forwards ? endpiece->range().end() : endpiece->range().begin();
          needsext = xtest.needsExtrapolation(*fittraj_,tdir);
        } while(needsext && fabs(time-tstart) < xtest.maxDt());
      }
    }
    return retval;
  }

  template<class KTRAJ> void Track<KTRAJ>::addDomain(Domain const& domain,TimeDir const& tdir,bool exact) {
    auto dptr = std::make_shared<Domain>(domain);
    if(tdir == TimeDir::forwards){
      auto const& prevdom = *domains_.crbegin();
      auto const& ktraj = fittraj_->nearestPiece(prevdom->end());
      auto& dw = effects_.emplace_back(std::make_unique<KKDW>(prevdom,dptr,ktraj));
      if(exact){
        dw->extrapolate(*fittraj_,tdir);
      } else {
        FitState fstate(ktraj.params());
        effects_.back()->process(fstate,tdir);
        effects_.back()->append(*fittraj_,tdir);
      }
    } else {
      auto const& nextdom = *domains_.cbegin();
      auto const& ktraj = fittraj_->nearestPiece(nextdom->begin());
      auto& dw = effects_.emplace_front(std::make_unique<KKDW>(dptr,nextdom,ktraj));
      if(exact){
        dw->extrapolate(*fittraj_,tdir);
      } else {
        FitState fstate(ktraj.params());
        effects_.front()->process(fstate,tdir);
        effects_.front()->append(*fittraj_,tdir);
      }
    }
    domains_.insert(dptr);
  }

  template<class KTRAJ> bool Track<KTRAJ>::extrapolate(EXINGPTR const& exingptr,TimeDir const& tdir) {
    bool retval = fitStatus().usable();
    if(retval){
      if(tdir == TimeDir::forwards){
        // make sure the time is legal
        retval = fittraj_->back().range().inRange(exingptr->time());
        if(retval){
          auto& exmat = effects_.emplace_back(std::make_unique<KKMAT>(exingptr,*fittraj_));
          exmat->extrapolate(*fittraj_,tdir);
        }
      } else {
        retval = fittraj_->front().range().inRange(exingptr->time());
        if(retval){
          auto& exmat = effects_.emplace_front(std::make_unique<KKMAT>(exingptr,*fittraj_));
          exmat->extrapolate(*fittraj_,tdir);
        }
      }
    }
    return retval;
  }


  template<class KTRAJ> TimeRange Track<KTRAJ>::activeRange() const {
    double tmin = std::numeric_limits<double>::max();
    double tmax = -std::numeric_limits<double>::max();
    // can't assume effects are sorted
    for(auto const& hit : hits_){
      tmin = std::min(tmin,hit->time());
      tmax = std::max(tmax,hit->time());
    }
    return TimeRange(tmin,tmax);
  }
}

#endif
