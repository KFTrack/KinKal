#ifndef KinKal_DomainWall_hh
#define KinKal_DomainWall_hh
//
// Effect describing the change in fit parameters for the change in BField crossing between 2 domains
// This effect adds no information content or noise (presently), just transports the parameters
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/Fit/Domain.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/Fit/Config.hh"
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
      void updateState(MetaIterConfig const& miconfig,bool first) override {}
      void updateConfig(Config const& config) override {}
      void updateReference(PTRAJ const& ptraj) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void append(PTRAJ& fit,TimeDir tdir) override;
      Chisq chisq(Parameters const& pdata) const override { return Chisq();} // no information added
      auto const& parameterChange() const { return dpfwd_; }
      virtual ~DomainWall(){}
      // disallow copy and equivalence
      DomainWall(DomainWall const& ) = delete;
      DomainWall& operator =(DomainWall const& ) = delete;
      // specific DomainWall interface
      // create from the domain and BField
      DomainWall(BFieldMap const& bfield,DOMAINPTR const& prevdomain,DOMAINPTR const& nextdomain, PTRAJ const& ptraj);
      // previous and next domains
      auto const& prevDomain() const { return *prev_; }
      auto const& nextDomain() const { return *next_; }

    private:
      BFieldMap const& bfield_; // bfield
      DOMAINPTR prev_, next_; // pointers to previous and next domains
      DVEC dpfwd_; // parameter change across this domain wall in the forwards time direction
  };

  template<class KTRAJ> DomainWall<KTRAJ>::DomainWall(BFieldMap const& bfield,
      DOMAINPTR const& prevdomain, DOMAINPTR const& nextdomain, PTRAJ const& ptraj) :
    bfield_(bfield), prev_(prevdomain), next_(nextdomain) {
      updateReference(ptraj);
    }

  template<class KTRAJ> void DomainWall<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    kkdata.append(dpfwd_,tdir);
    // rotate the covariance matrix for the change in reference BField.  TODO
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // update BField at the domain centers given the new trajectory
    prev_->updateBNom(bfield_.fieldVect(ptraj.position3(prev_->mid())));
    next_->updateBNom(bfield_.fieldVect(ptraj.position3(next_->mid())));
    // update change in parameters for passage through this DW
    auto const& refpiece = ptraj.nearestPiece(prev_->mid()); // by convention, use previous domains parameters to define the derivative
    dpfwd_ = refpiece.dPardB(this->time(),next_->bnom());
  }

  template<class KTRAJ> void DomainWall<KTRAJ>::append(PTRAJ& ptraj,TimeDir tdir) {
    double etime = time();
    // make sure the piece is appendable
    if((tdir == TimeDir::forwards && ptraj.back().range().begin() > etime) ||
        (tdir == TimeDir::backwards && ptraj.front().range().end() < etime) )
      throw std::invalid_argument("DomainWall: Can't append piece");
    auto const& oldpiece = (tdir == TimeDir::forwards) ? ptraj.back() : ptraj.front();
    KTRAJ newpiece(oldpiece);
    newpiece.range() = (tdir == TimeDir::forwards) ? TimeRange(next_->begin(),std::max(ptraj.range().end(),next_->end())) :
      TimeRange(std::min(ptraj.range().begin(),prev_->begin()),prev_->end());
    if( tdir == TimeDir::forwards){
      newpiece.params() += dpfwd_;
      newpiece.bnom() = next_->bnom(); // set the parameters according to what was used in processing
      newpiece.setBNom(next_->mid(),bfield_.fieldVect(ptraj.position3(next_->mid()))); // update to reference the BField at the next
      ptraj.append(newpiece);
    } else {
      newpiece.params().parameters() -= dpfwd_; // reverse effect going backwards.  Should also rotate covariance TODO
      newpiece.bnom() = prev_->bnom();
      newpiece.setBNom(prev_->mid(),bfield_.fieldVect(ptraj.position3(prev_->mid())));
      ptraj.prepend(newpiece);
    }
    // update for the next iteration
    prev_->updateBNom(bfield_.fieldVect(oldpiece.position3(prev_->mid())));
    next_->updateBNom(bfield_.fieldVect(newpiece.position3(next_->mid())));
    dpfwd_ = oldpiece.dPardB(etime,next_->bnom());
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

}
#endif
