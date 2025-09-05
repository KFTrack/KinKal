#ifndef KinKal_DomainWall_hh
#define KinKal_DomainWall_hh
//
// Effect describing the change in fit parameters for the change in BField crossing between 2 domains
// This effect adds no information content or noise (presently), just transports the parameters
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Fit/Domain.hh"
#include "KinKal/Fit/Config.hh"
#include "KinKal/Fit/FitState.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/BFieldMap.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class DomainWall : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using DOMAINPTR = std::shared_ptr<Domain>;
      double time() const override { return prev_->end(); }
      bool active() const override { return true; } // always active
      void process(FitState& kkdata,TimeDir tdir) override;
      void updateState(MetaIterConfig const& miconfig,bool first) override;
      void updateConfig(Config const& config) override {}
      void updateReference(PTRAJ const& ptraj) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void append(PTRAJ& fit,TimeDir tdir) override;
      void extrapolate(PTRAJ& fit,TimeDir tdir) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();} // no information added
      auto const& parameterChange() const { return dpfwd_; }
      virtual ~DomainWall(){}
      // disallow copy and equivalence
//    DomainWall(DomainWall const& ) = delete;
      DomainWall& operator =(DomainWall const& ) = delete;
      // clone op for reinstantiation
      DomainWall(DomainWall const&);
      std::unique_ptr< Effect<KTRAJ> > clone(CloneContext&) const override;
      // specific DomainWall interface
      DomainWall(DOMAINPTR const& prevdomain,DOMAINPTR const& nextdomain, PTRAJ const& ptraj);
      // previous and next domains
      auto const& prevDomain() const { return *prev_; }
      auto const& nextDomain() const { return *next_; }
      // other accessors
      auto const& prevPtr() const { return prev_; }
      auto const& nextPtr() const { return next_; }
      auto const& fwdChange() const { return dpfwd_; }
      auto const& prevWeight() const { return prevwt_; }
      auto const& nextWeight() const { return nextwt_; }
      auto const& fwdCovarianceRotation() const { return dpdpdb_; }
      void setPrevPtr(DOMAINPTR const& ptr){ prev_ = ptr; }
      void setNextPtr(DOMAINPTR const& ptr){ next_ = ptr; }

    private:
      DOMAINPTR prev_, next_; // pointers to previous and next domains
      DVEC dpfwd_; // parameter change across this domain wall in the forwards time direction
      Weights prevwt_, nextwt_; // cache of weights
      PSMAT dpdpdb_; // forward rotation of covariance matrix going in the forwards direction

  };

  template<class KTRAJ> DomainWall<KTRAJ>::DomainWall(
      DOMAINPTR const& prevdomain, DOMAINPTR const& nextdomain, PTRAJ const& ptraj) :
    prev_(prevdomain), next_(nextdomain) {
      updateReference(ptraj);
    }

  template<class KTRAJ> void DomainWall<KTRAJ>::process(FitState& fstate,TimeDir tdir) {
    if(tdir == TimeDir::forwards) {
      prevwt_ += fstate.wData();
      fstate.append(dpfwd_,tdir);
      // fstate.pData().covariance() = ROOT::Math::Similarity(dpdpdb_,fstate.pData().covariance());  Not tested TODO
      nextwt_ += fstate.wData();
    } else {
      nextwt_ += fstate.wData();
      fstate.append(dpfwd_,tdir);
      // fstate.pData().covariance() = ROOT::Math::SimilarityT(dpdpdb_,fstate.pData().covariance());
      prevwt_ += fstate.wData();
    }
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // reset the cached weights
    prevwt_ = nextwt_ = Weights();
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // update change in parameters for passage through this DW
    auto const& refpiece = ptraj.nearestPiece(time()-1e-5); // disambiguate derivativates
    auto db = next_->bnom() - prev_->bnom();
    dpfwd_ = refpiece.dPardB(this->time(),db);
    dpdpdb_ = refpiece.dPardPardB(this->time(),db);
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
    // make sure the piece is appendable
    if((tdir == TimeDir::forwards && ptraj.back().range().begin() > time()) ||
        (tdir == TimeDir::backwards && ptraj.front().range().end() < time()) )
      throw std::invalid_argument("DomainWall: Can't append piece");
    auto const& oldpiece = (tdir == TimeDir::forwards) ? ptraj.back() : ptraj.front();
    KTRAJ newpiece(oldpiece);
    if( tdir == TimeDir::forwards){
      newpiece.range() = TimeRange(next_->begin(),std::max(ptraj.range().end(),next_->end()));
      newpiece.params() = Parameters(nextwt_);
      newpiece.resetBNom(next_->bnom());
      ptraj.append(newpiece);
    } else {
      newpiece.range() = TimeRange(std::min(ptraj.range().begin(),prev_->begin()),prev_->end());
      newpiece.params() = Parameters(prevwt_);
      newpiece.resetBNom(prev_->bnom());
      ptraj.prepend(newpiece);
    }
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::extrapolate(PTRAJ& ptraj,TimeDir tdir) {
  // make sure the piece is appendable
    if((tdir == TimeDir::forwards && ptraj.back().range().begin() > time()) ||
        (tdir == TimeDir::backwards && ptraj.front().range().end() < time()) )
      throw std::invalid_argument("DomainWall: Can't append piece");
    // sample the particle state at this domain wall
    auto pstate = ptraj.stateEstimate(time());
    if( tdir == TimeDir::forwards){
      // re-intepret that as a trajectory, using the next domain's magnetic field. This gives exact continuity in position and momentum
      KTRAJ newpiece(pstate, nextDomain().bnom(),nextDomain().range());
      ptraj.append(newpiece);
    } else {
      // same as above, in the opposite direction
      KTRAJ newpiece(pstate, prevDomain().bnom(),prevDomain().range());
      ptraj.prepend(newpiece);
    }
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "DomainWall " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " previous domain " << *prev_ << " next domain " << *next_;
    ost << " effect " << dpfwd_ << std::endl;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, DomainWall<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

  // clone op for reinstantiation
  template <class KTRAJ>
  DomainWall<KTRAJ>::DomainWall(DomainWall const& rhs):
      dpfwd_(rhs.fwdChange()),
      prevwt_(rhs.prevWeight()),
      nextwt_(rhs.nextWeight()),
      dpdpdb_(rhs.fwdCovarianceRotation()){
    /**/
  }

  template <class KTRAJ>
  std::unique_ptr< Effect<KTRAJ> > DomainWall<KTRAJ>::clone(CloneContext& context) const{
    auto casted = std::make_unique< DomainWall<KTRAJ> >(*this);
    DOMAINPTR prev = context.get(prev_);
    casted->setPrevPtr(prev);
    DOMAINPTR next = context.get(next_);
    casted->setNextPtr(next);
    //auto rv = std::make_unique< Effect<KTRAJ> >(casted);
    auto rv = std::move(casted);
    return rv;
  }
}
#endif
