#ifndef KinKal_KKTrk_hh
#define KinKal_KKTrk_hh
//
//  Primary class of the Kinematic Kalman fit.  This class owns the state describing
//  the fit inputs (hits, material interactions, BField corrections, etc), the result of the fit,
//  and the methods for computing it.  The fit result is expressed as a piecewise kinematic covariant
//  particle trajectory, providing position, momentum etc information about the particle with covariance
//  as a function of physical time.
//
//  KKTrk is templated on a simple kinematic trajectory class representing the 1-dimensional path and
//  momentum of a particle traveling through empty space in a constant magnetic field, as a function of time.
//  The piecewise kinematic trajectory fit result is expressed as a time sequence of these simple trajectory objects.
//  Material effects and spatial variation of magnetic fields are modeled through changes between adjacent simple
//  trajectories.
//  To instantiate KKTrk the particle trajectory class must satisfy the following interface:
//      void position(Vec4& pos) const;
//      void position(float time, Vec3& pos) const;
//      void velocity(float time, Vec3& vel) const;
//      double speed(float time) const;
//      void direction(float time, Vec3& dir) const;
//      void print(std::ostream& ost, int detail) const;
//      void momentum(double t,Mom4& mom) const; // momentum in MeV/c, mass in MeV/c^2 as a function of time
//      void momentum(Vec4 const& pos, Mom4& mom) const { return momentum(pos.T(),mom); }
//      double momentum(float time) const; // momentum and energy magnitude in MeV/
//      double momentumVar(float time) const; // variance on momentum value
//      double energy(float time) const; 
//      void rangeInTolerance(TRange& range, BField const& bfield, double tol);
//      PDATA const& params() const;
//
//  The PDATA object provides a minimal basis from which the geometric and kinematic properties of the particle as a function
//  of time can be computed.  For instance, a kinematic helix in space requires a PDATA instance with 6 parameters.  The physical
//  interpretation of the PDATA payload is made in the kinematic trajectory class.
//
//  KKTrk uses the root SVector and SMatrix classes for algebraic manipulation, and GenVector classes for geometric and
//  kinematic particle descriptions, both part of the root Math package.  These are described on the root website https://root.cern.ch/root/html608/namespaceROOT_1_1Math.html
//
//  The underlying processing model is a progressive BLUE fit first used in the geometric track fit implementation used by the BaBar
//  collaboration, described in "D.N. Brown, E.A. Charles, D.A. Roberts, The BABAR track fitting algorithm, Proceedings of CHEP 2000, Padova, Italy, 2000"
//
//  KKTrk is constructed from a configuration object which can be shared between many instances, and a unique set of hit and
//  material interactions.  The configuration object controls the fit iteration convergence testing, including simulated
//  annealing and interactions with the external environment such as the material model and the magnetic field map.
//  The fit is performed on construction.
//
//  The KinKal package is licensed under Adobe v2, and is hosted at 
//
#include "KinKal/PKTraj.hh"
#include "KinKal/KKEff.hh"
#include "KinKal/KKEnd.hh"
#include "KinKal/KKMHit.hh"
#include "KinKal/KKHit.hh"
#include "KinKal/KKMat.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/THit.hh"
#include "KinKal/KKConfig.hh"
#include "KinKal/FitStatus.hh"
#include "KinKal/BField.hh"
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
  template<class KTRAJ> class KKTrk {
    public:
      typedef KKEff<KTRAJ> KKEFF;
      typedef std::shared_ptr<KKConfig> KKCONFIGPTR;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef THit<KTRAJ> THIT;
      typedef DXing<KTRAJ> DXING;
      typedef std::shared_ptr<THIT> THITPTR;
      typedef std::shared_ptr<DXING> DXINGPTR;
      typedef std::vector<THITPTR> THITCOL;
      typedef std::vector<DXINGPTR> DXINGCOL;
      struct KKEFFComp { // comparator to sort effects by time
	bool operator()(std::unique_ptr<KKEFF> const& a, std::unique_ptr<KKEFF> const&  b) const {
	  if(a.get() != b.get())
	    return a->time() < b->time();
	  else
	    return false;
	}
      };
      typedef std::vector<std::unique_ptr<KKEFF> > KKEFFCOL; // container type for effects
      // construct from a set of hits and passive material crossings
      KKTrk(KKCONFIGPTR const& kkconfig, PKTRAJ const& reftraj, THITCOL& thits, DXINGCOL& dxings ); 
      void fit(); // process the effects.  This creates the fit
      // accessors
      std::vector<FitStatus> const& statusHistory() const { return history_; }
      FitStatus const& fitStatus() const { return history_.back(); } // most recent status
      PKTRAJ const& refTraj() const { return reftraj_; }
      PKTRAJ const& fitTraj() const { return fittraj_; }
      KKEFFCOL const& effects() const { return effects_; }
      KKConfig const& config() const { return *kkconfig_; }
      THITCOL const& timeHits() const { return thits_; } 
      DXINGCOL const& detMatXings() const { return dxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
    private:
      // helper functions
      void update(FitStatus const& fstat, MConfig const& mconfig);
      void fitIteration(FitStatus& status);
      bool canIterate() const;
      bool oscillating(FitStatus const& status) const;
      void createFieldDomains();
      // payload
      KKCONFIGPTR kkconfig_; // shared configuration
      std::vector<FitStatus> history_; // fit status history; records the current iteration
      PKTRAJ reftraj_; // reference against which the derivatives were evaluated and the current fit performed
      PKTRAJ fittraj_; // result of the current fit, becomes the reference when the fit is algebraically iterated
      KKEFFCOL effects_; // effects used in this fit, sorted by time
      THITCOL thits_; // shared collection of hits
      DXINGCOL dxings_; // shared collection of material crossings/interactions
  };

// construct from configuration, reference (seed) fit, hits,and materials specific to this fit.  Note that hits
// can contain associated materials.
  template <class KTRAJ> KKTrk<KTRAJ>::KKTrk(KKCONFIGPTR const& kkconfig, PKTRAJ const& reftraj,  THITCOL& thits, DXINGCOL& dxings) : 
    kkconfig_(kkconfig), reftraj_(reftraj), fittraj_(reftraj.range(),reftraj.mass(),reftraj.charge()), 
    thits_(thits), dxings_(dxings) {
    // create the effects.  First, loop over the hits
      for(auto& thit : thits_ ) {
	// create the hit effects and insert them in the set
	// if there's associated material, create a combined material and hit effect, otherwise just a hit effect
	if(thit->hasMaterial()){
	  dxings_.push_back(thit->detCrossing());
	  effects_.emplace_back(std::make_unique<KKMHit<KTRAJ> >(thit,reftraj));
	} else{ 
	  effects_.emplace_back(std::make_unique<KKHit<KTRAJ> >(thit,reftraj));
	}
      }
      //add pure material effects
      for(auto& dxing : dxings) {
	effects_.emplace_back(std::make_unique<KKMat<KTRAJ> >(dxing,reftraj));
      }
      // reset the range 
      reftraj_.setRange(TRange(std::min(reftraj_.range().low(),effects_.begin()->get()->time() - config().tbuff_),
	    std::max(reftraj_.range().high(),effects_.rbegin()->get()->time() + config().tbuff_)));
      // add BField inhomogeneity effects; not yet implemented FIXME!
      // createFieldDomains();
      // create end effects; this should be last to avoid confusing the BField correction
      effects_.emplace_back(std::make_unique<KKEnd<KTRAJ>>(reftraj,TDir::forwards,config().dwt_));
      effects_.emplace_back(std::make_unique<KKEnd<KTRAJ>>(reftraj,TDir::backwards,config().dwt_));
      // now fit the track
      fit();
    }

  // fit iteration management 
  template <class KTRAJ> void KKTrk<KTRAJ>::fit() {
   // execute the schedule of meta-iterations
    for(auto imconfig=config().schedule().begin(); imconfig != config().schedule().end(); imconfig++){
      auto mconfig  = *imconfig;
      mconfig.miter_  = std::distance(config().schedule().begin(),imconfig);
      // algebraic convergence iteration
      FitStatus fstat(mconfig.miter_);
      history_.push_back(fstat);
      while(canIterate()) {
	update(fstat,mconfig);
	fitIteration(fstat);
      }
    }
  }

  // single algebraic iteration 
  template <class KTRAJ> void KKTrk<KTRAJ>::fitIteration(FitStatus& fstat) {
  // catch exceptions and record them in the status
    try {
    // reset counters
      fstat.chisq_ = 0.0;
      fstat.ndof_ = -(int)KTRAJ::NParams();
      fstat.iter_++;
      // fit in both directions (order doesn't matter)
      auto feff = effects_.begin();
      // start with empty fit information; each effect will modify this as necessary, and cache what it needs for later processing
      KKData<KTRAJ::NParams()> ffitdata;
      while(feff != effects_.end()){
	auto ieff = feff->get();
	// update chisquared; only needed forwards
	fstat.ndof_ += ieff->nDOF();
	fstat.chisq_ += ieff->chisq(ffitdata.pData());
	// process
	ieff->process(ffitdata,TDir::forwards);
	feff++;
      }
      fstat.prob_ = TMath::Prob(fstat.chisq_,fstat.ndof_);
      // reset the fit information and process backwards (the order does not matter)
      KKData<KTRAJ::NParams()> bfitdata;
      auto beff = effects_.rbegin();
      while(beff != effects_.rend()){
	auto ieff = beff->get();
	ieff->process(bfitdata,TDir::backwards);
	beff++;
      }
      // convert the fit result into a new trajectory; start with an empty ptraj
      fittraj_ = PKTRAJ(reftraj_.range(),reftraj_.mass(),reftraj_.charge());
      // process forwards, adding pieces as necessary
      for(auto& ieff : effects_) {
	ieff->append(fittraj_);
      }
      // trim the range to the physical elements (past the end sites)
      feff = effects_.begin(); feff++;
      beff = effects_.rbegin(); beff++;
      fittraj_.range().low() = (*feff)->time() - config().tbuff_;
      fittraj_.range().high() = (*beff)->time() + config().tbuff_;
      // update status
      if(fabs(fstat.chisq_ -fitStatus().chisq_) < config().convdchisq_) {
	fstat.status_ = FitStatus::converged;
      } else if (fstat.chisq_-fitStatus().chisq_ > config().divdchisq_) {
	fstat.status_ = FitStatus::diverged;
      } else if (fstat.ndof_ < config().minndof_){
	fstat.status_ = FitStatus::lowNDOF;
      } else if(oscillating(fstat)){
	fstat.status_ = FitStatus::oscillating;
      } else
	fstat.status_ = FitStatus::unconverged;
    } catch (std::exception const& error) {
      fstat.status_ = FitStatus::failed;
      fstat.comment_ = error.what();
    }
    // record this status in the history
    history_.push_back(fstat);
  }

  // update between iterations 
  template <class KTRAJ> void KKTrk<KTRAJ>::update(FitStatus const& fstat, MConfig const& mconfig) {
    if(fstat.iter_ < 0) { // 1st iteration of a meta-iteration: update the state
      if(mconfig.miter_ > 0)// if this isn't the 1st meta-iteration, swap the fit trajectory to the reference
	reftraj_ = fittraj_;
      for(auto& ieff : effects_ ) ieff->update(reftraj_,mconfig);
    } else {
    //swap the fit trajectory to the reference
      reftraj_ = fittraj_;
    // update the effects to use the new reference
      for(auto& ieff : effects_) ieff->update(reftraj_);
    }
// sort the effects by time
    std::sort(effects_.begin(),effects_.end(),KKEFFComp ());
  }

  template<class KTRAJ> bool KKTrk<KTRAJ>::canIterate() const {
    auto const& fstat = fitStatus();
    if( (fstat.status_ == FitStatus::needsfit ||
	  fstat.status_ == FitStatus::unconverged) &&
	fstat.iter_ < config().maxniter_) {
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
    //  double tstep = reftraj_.timeStep(reftraj_.range().mid(),config().dtol_);
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
    //    float time = domain.range_.mid();
    //    Vec3 tpos,fvec;
    //    reftraj_.position(time,tpos);
    //    bfield_.fieldVect(fvec,tpos);
    //
    //    // integrate 
    //  }
  }

  template <class KTRAJ> void KKTrk<KTRAJ>::print(std::ostream& ost, int detail) const {
    using std::endl;
    if(detail == 0) 
      ost <<  "Fit Status " << fitStatus();
    else {
      ost <<  "Fit History " << endl;
      for(auto const& stat : statusHistory()) ost << stat << endl;
    }
    ost << "Fit Result ";
    fitTraj().print(ost,detail);
    if(detail > 1) {
      ost << "Reference ";
      refTraj().print(ost,detail);
    }
    if(detail > 2) {
      ost << "Effects " << endl;
      for(auto const& eff : effects()) eff.get()->print(ost,detail);
    }
  }

}
#endif
